%% 
clear all
close all
clc

% test parameters
K_Val           = [1,2,4,8]; 	% the number of stages for the sequential test. Test performance was evaluated for different K
Alpha           = 0.025;      	% the type-I error rate for the full K-staged sequential test
Beta            = 0.1;          % the type-II error rate for the full K-staged sequential test
Gamma         	= 0.9;       	% the statistical power for the full K-staged sequential test
Gamma_max       = 0.9;          % the required statistical power for removing the assumed effect sizes 
ES_True         = 0.25;         % the true effect size for Simulation II. For Simulations I, ES_True = 1
ES_Assumed      = [0.25, 0.5, 1];   % the assumed effect sizes for Simulation II. For Simulations I, ES_Assumed = 1
Na_Val          = [10, 20, 50];     % the number of accrued samples before data monitoring is initiated. Test performance was evaluated for different "Na".
Sigma2_Val      = 1:50;             % the Sample variance for the simulated data. Test performance was evaluated for different "Sigma2"
Q               = 1;                % the number of features for the simulated data set. The obne-sample t-test is used for Q=1
Res_parameter   = 10000;            % a parameter for controlling the resolution of the numerical distributions 

% the following simulations were run 25,000 times using the Iridis supercomputer at Southampton University. 
for i0=1:length(Na_Val)
    for i1=1:length(Sigma2_Val)
        for i2=1:length(K_Val)
            
            % unpack
            Na      = Na_Val(i0);           % the number of accrued samples before data monitoring is initiated
            Sigma2	= Sigma2_Val(i1);       % the variance of the data to be simulated
            K       = K_Val(i2);            % the number of stages for the sequential analysis
            Alpha_k = ones(1,K) * (Alpha / K);                      % the stage-wise type-I error rates
            Beta_k  = ones(1,K) * (Beta / K);                       % the stage-wise type-II error rates
            Gamma_k	= Get_Gamma_k( K_Val(i2), Q, Alpha, Gamma);     % the stage-wise statistical powers (cumulative across stages)

            % initialise 
            k               = 0;        % the curernt stage of the analysis
            AllData         = [];       % all accrued data, up to stage k
            T_k             = [];       % the stage-wise test statistics
            StagePowers     = [];       % the estimated stage-wise statistical powers
            
            % reset
            clear FAxis_Values_Alt Alt AltTrunc     % delete distributions and parameters from previously simulated trials
            All_ES          = ES_Assumed;           % this keeps track of the assumed effect sizes throughout the analysis (the assumed effect sizes can be modified throughout the trial)
            Resolution_Alt  = [];                   % the resolution or smoothness of the alternative distributions, to be determined below. 
            ContinueTesting	= true;         % testing continues untill H0 has been accepted or rejected, or untill the final stage of the analysis has not yet been reached.
            while ContinueTesting            

                k           = k + 1;            % increment stage index
                Data        = [];               % no data yet for stage k
                N_k         = 0;                % sample size for stage k is 0
                AddN        = 1;                % monitor the sample variance after every "AddN" new samples, starting with 1 (this will be increased as the trial progresses)
                Monitoring  = true;         
                while Monitoring            % sample variance is monitored untill the desired statistical power for stage k has been obtained

                    if N_k > 50
                        AddN = ceil( N_k * 0.05 );  % as the trial progresses and N_k increases, data is monitored less often (AddN is increased) to reduce computation time
                    end
                    if N_k == 0
                        Data = ES_True + (sqrt(Sigma2) * randn(Na,Q));          	% generate "Na" new observations
                    else
                        Data = [Data; ES_True + (sqrt(Sigma2) * randn(AddN,Q)) ];	% generate "AddN" new observations
                    end
                    N_k = size(Data,1);                                             % current sample size for stage k

                    % generate the null distribution 
                    if k==1                                                     
                        FRange_Null         = finv(1-0.000001, Q, N_k-Q);               	% the upper boundary for the axis along which the null distribution is defined
                        FAxis_Values_Null	= 0:(FRange_Null/Res_parameter):FRange_Null;	% the axis values along which the distribution is defined
                        Resolution_Null     = 1 / (FRange_Null/Res_parameter);              % the resolution of the axis 
                    end
                    Null = fpdf(FAxis_Values_Null, Q, N_k-Q);               % the null distribution for stage k 
                    if Q==1
                        Null(1) = Null(2) + (Null(2)-Null(3));              % avoid infinities
                    end
                    Null = Null/sum(Null);                                  % normalise so the area equals 1
                    if k>1                                                  
                        Null = conv(Null, NullTrunc);                       % for stages 2 and later, convolve with the truncated null distribution from the previous stage
                        Null = (Null ./ sum(Null)) * sum(NullTrunc);        % normalise the area so that it equals 1 minus the area of the trunacted regions
                    end
                    
                    % generate the efficacy threshold
                    Percentile          = sum(Null) - Alpha_k(k);                                 	% the percentile associated with the efficacy threshold
                    Efficacy_Ind_Null	= Get_Truncation_Index(Null, Percentile, Resolution_Null);  % find the location of the efficacy threshold along the null distribution
                    A_k                 = Efficacy_Ind_Null(1)/Resolution_Null;                     % the efficacy threshold

                    % generate alternative distributions and estimate the statistical powers
                    CovMatrix       = cov( [AllData;Data] );        % the covariance matrix to monitor
                    InvCov          = inv( CovMatrix );             % and its inverse 
                    for ESi = 1:length(All_ES)
                        NC_Estimate	= N_k * All_ES(ESi) * InvCov * All_ES(ESi);                     % the estimated non-centrality parameter for effect size ESi
                        if k==1         
                            FRange_Alt              = ncfinv(1-0.000001, Q, N_k-Q, NC_Estimate );   % the upper boundary for the axis values along which the alternative distribution will be defined
                            FAxis_Values_Alt{ESi}	= 0:(FRange_Alt/Res_parameter):FRange_Alt;      % the axis values along which the alternative distribution is defined
                            Resolution_Alt(ESi)     = 1 / (FRange_Alt/Res_parameter);               % and the axis resolution
                        end
                        Alt{ESi} = ncfpdf(FAxis_Values_Alt{ESi}, Q, N_k-Q, NC_Estimate );           % the alternative distribution, given by a non-central F-distribution
                        if Alt{ESi}(1)==inf      
                            Alt{ESi}(1) = Alt{ESi}(2) + (Alt{ESi}(2)-Alt{ESi}(3));                  % avoid infinities
                        end
                        Alt{ESi} = Alt{ESi} / sum(Alt{ESi});                                        % normalise so the area is equal to 1
                        if k>1                                                                      % from stage 2 onwards, convolve with the truncated distribution from the previous stage
                            Alt{ESi} = conv( Alt{ESi}, AltTrunc{ESi} );                             % carry out the convolution
                            Alt{ESi} = ( Alt{ESi} ./ sum( Alt{ESi}) ) * sum( AltTrunc{ESi} );       % and normalise so the area equals 1 minus the areas of the truncated regions
                        end

                        % estimate power
                        Efficacy_Ind_Alt    = round( A_k * Resolution_Alt(ESi) );           % the location along the alternative distribution associated with the efficacy threshold A_k
                        StagePowers(k, ESi) = sum( Alt{ESi}(Efficacy_Ind_Alt(1):end) );     % statistical power: the area under the alternative to the right of the A_k threshold

                    end

                    % evalaute stop condition
                    if mean( sum(StagePowers,1) ) >= Gamma_k(k)     % if sufficient power has been obtained
                        Monitoring = false;                         % stop data monitoring
                    end
                    
                    % display the estimated powers under the assumed effect sizes
                    StagePowers
                end
                
                % generate the B_k futility threshold
                ESi  = 1;	% the index associated with the smallest anticipated effect size
               	if Alt{ESi}(1) > Beta_k(k)      % the area under the alternative distribution at index 1 already exceeds the specified type-II error rate for this stage
                    B_k         = 0;            % futility threshold set to zero
                    try                         % and the type-II error rate for this stage is bumped up to the next stage
                        Beta_k(k+1) = Beta_k(k+1) + Beta_k(k);  
                    end
                    Beta_k(k)   = 0;            % 
                else
                    Futility_Ind_Null	= Get_Truncation_Index( Alt{ESi}, Beta_k(k), Resolution_Alt(ESi) );     % the index associated with the location of the futility threshold
                    B_k                 = Futility_Ind_Null(1)/Resolution_Alt(ESi);                          	% and the actual futility threshold
                    try         % update the stage-wise type-II error rates in accordance with numerical rounding errors 
                        Beta_k(k+1) = Beta_k(k+1) + (Beta_k(k) - sum( Alt{ESi}(1:Futility_Ind_Null(1)) ) );
                    end
                end
                
                % analyse the data
                T_k(k) = Calculate_F(Data);     % the stage k test statistic
                S_k(k)  = sum(T_k);             % and the stage k summary statistic
                
                % update all accrued data 
                AllData = [AllData;Data];

                % check efficacy stopping
                if S_k(k) >= A_k              	% reject H0
                    ContinueTesting = false;
                else
                    % check futility stopping
                    if S_k(k) <= B_k            % accept H0
                        ContinueTesting = false;
                    end
                end
                % check if final stage reached
                if k == K                           
                    ContinueTesting = false;             
                end

                if ContinueTesting  % prepare for the next stage

                    % truncations for the null distribution
                    if B_k == 0     
                        Futility_Ind_Null = 0;                           	% index associated with the B_k threshold
                    else
                        Futility_Ind_Null = round( B_k * Resolution_Null );	% index associated with the B_k threshold
                    end
                    NullTrunc                           = Null;             % prepare to truncate
                    NullTrunc(Efficacy_Ind_Null(1):end)	= zeros( 1, length(NullTrunc(Efficacy_Ind_Null(1):end)) );	% all values larger than the efficacy threshold are set to zero
                    NullTrunc(1:Futility_Ind_Null(1))	= zeros( 1, length(NullTrunc(1:Futility_Ind_Null(1))) );    % and all values smaller than futility threshold are set to zero

                    % truncations for the alternative distributions and the effect size 
                    clear Alt_New AltTrunc_New      % reset
                    ES_New          = [];   % temporary array to store the effect sizes that are to be carried foward to the next stage
                    StagePowers_New = [];   % and a temporary array to store the estimated statistical powers for the effect sizes
                    C               = 0;    % counter to keep track of how many effect sizes to be taken foward to the next stage
                    for ESi = 1:length(All_ES)
                        if sum(StagePowers(:,ESi)) < Gamma_max  % if the estimated power under this effect size has not exceeded Gamma_max...
                            C                       = C + 1;                % increment C
                            ES_New(C)               = All_ES(ESi);          % store the effect size
                            Alt_New{C}              = Alt{ESi};             % store the alternative distribution
                            StagePowers_New(:,C)	= StagePowers(:,ESi);   % and store the estimated powers under this effect size
                            Efficacy_Ind_Alt        = round( A_k * Resolution_Alt(ESi) ); 	% find the location of the futility threshold along the alternative distribution for this effect size
                            if B_k == 0
                                Futility_Ind_Alt    = 0;                                    % the location of the futility threshold
                            else
                                Futility_Ind_Alt	= round( B_k * Resolution_Alt(ESi) );   % the location of the futility threshold
                            end
                            % truncate the distribution
                            AltTrunc_New{C}                             = Alt{ESi};     % prepare for truncation
                            AltTrunc_New{C}(Efficacy_Ind_Alt(1):end)	= zeros( 1, length(AltTrunc_New{C}(Efficacy_Ind_Alt(1):end)) );     % all values larger than the efficacy threshold are set to zero
                            AltTrunc_New{C}(1:Futility_Ind_Alt(1))  	= zeros( 1, length(AltTrunc_New{C}(1:Futility_Ind_Alt(1))) );    	% and all values smaller than futility threshold are set to zero
                        end
                    end
                    % replace 
                    Alt         = Alt_New;
                    AltTrunc	= AltTrunc_New;
                    All_ES      = ES_New;
                    StagePowers = StagePowers_New;
                    
                end
            end
            
            % store results
            % ...
            
        end
    end
end

