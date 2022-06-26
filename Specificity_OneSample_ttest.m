%%
clear all
close all
clc

% test parameters
K_Val           = [1,2,4,8]; 	% the number of stages for the sequential test. Test performance was evaluated for different K
Alpha           = 0.025;      	% the type-I error rate for the full K-staged sequential test
Gamma         	= 0.9;       	% the statistical power for the full K-staged sequential test
ES_Assumed      = 1;                % the assumed simulated effect size
Na_Val          = [10, 20, 50];     % the number of accrued samples before data monitoring is initiated. Test performance was evaluated for different "Na".
Sigma2_Val      = 1:50;             % the sample variance for the simulated data. Test performance was evaluated for different "Sigma2"
Q               = 1;                % the number of features for the simulated data set. The obne-sample t-test is used for Q=1
Resolution      = (1/0.0005);        	% the resolution for the axis along which the distributions are generated
XAxis_Values	= 0:1/Resolution:60;	% x-axis values along which the distributions are generated

% the following simulations were run 200,000 times using the Iridis supercomputer at Southampton University. 
for i0=1:length(Na_Val)
    for i1=1:length(Sigma2_Val)
        for i2=1:length(K_Val)
            
            % unpack
            Na          = Na_Val(i0);           % the number of accrued samples before data monitoring is initiated
            Sigma2      = Sigma2_Val(i1);       % the variance of the data to be simulated
            K           = K_Val(i2);            % the number of stages for the sequential analysis
            Alpha_k      = ones(1,K) * (Alpha / K);                     % the stage-wise type-I error rates
            Gamma_k     = Get_Gamma_k( K_Val(i2), Q, Alpha, Gamma );    % The stage-wise statistical powers (cumulative across stages)
            
            % initialise 
            k               = 0;        % the curernt stage of the analysis
            AllData         = [];       % all accrued data, up to stage k
            T_k             = [];       % the stage-wise test statistics
            StagePowers     = [];       % the estimated stage-wise statistical powers
            ContinueTesting	= true;
            while ContinueTesting       % the sequential test starts here

                k           = k + 1; 	% increment the stage index
                N_k         = 0;        % sample size for stage k is 0
                Data        = [];       % no data collected at this stage yet
                AddN        = 1;        % monitor the sample variance after every "AddN" new samples, starting with 1 (this will be increased as the trial progresses)
                Monitoring  = true;     
                while Monitoring      	% continuously monitor sample variance 

                    % generate data
                    if N_k == 0
                        Data	= sqrt(Sigma2) * randn(Na,Q);          	% add Na observations
                    else
                        Data    = [Data; sqrt(Sigma2) * randn(AddN,Q)]; % add "AddN" observation
                    end
                    N_k = size(Data,1);

                    % increment AddN
                    if N_k > 50
                        AddN = ceil( N_k * 0.05 );              % as the trial progresses, data is monitored less often to reduce computation time
                    end
                    
                    % generate the null distribution
                    Null = fpdf(XAxis_Values, Q, N_k-Q);        % the null distribution for stage k 
                    if Q==1
                        Null(1) = Null(2) + (Null(2)-Null(3));  % avoid infinities
                    end
                    Null        = Null/sum(Null);               % normalise so the area equals 1
                    if k>1                                                  % for stage two onwards... 
                        Null	= conv(Null, NullTrunc);                    % convolve with the truncated null distribution from the previous stage
                        Null	= (Null ./ sum(Null)) * sum(NullTrunc); 	% normalise the area so that it equals 1 minus the area of the trunacted regions
                    end

                    % generate the efficacy threshold
                    Percentile      = sum(Null) - Alpha_k(k);                               % the percentile associated with the efficacy threshold
                    Efficacy_Ind	= Get_Truncation_Index(Null, Percentile, Resolution);   % find the location of the efficacy threshold along the null distribution
                    A_k             = Efficacy_Ind(1)/Resolution;                           % the efficacy threshold

                    % generate the alternative distribution
                    CovMatrix       = cov( [AllData;Data] );    % the covariance matrix to monitor
                    InvCov          = inv( CovMatrix );         % and its inverse
                    NonCen_Estimate	= N_k * ES_Assumed * InvCov * ES_Assumed';  % the non-centrality parameter for this effect size           
                    Alt             = ncfpdf( XAxis_Values, Q, N_k - Q, NonCen_Estimate );
                    if Alt(1)==inf
                        Alt(1) = Alt(2) + (Alt(2)-Alt(3));  % avoid infinities
                    end
                    Alt             = Alt / sum(Alt);       % normalise so the area equals 1
                    if k>1                                              % for stage two onwards... 
                        Alt     = conv( Alt, AltTrunc );                % convolve with the truncated null distribution from the previous stage
                        Alt     = (Alt ./ sum(Alt)) * sum( AltTrunc );  % and normalise the area so that it equals 1 minus the area of the trunacted regions
                    end

                    % estimate statistical power
                    StagePowers(k) = sum( Alt(Efficacy_Ind(1):end) );   % estimated power (cumulative across stages)
                    if sum( StagePowers ) >= Gamma_k(k)                 % if the power exceeds the pre-specified power level
                        Monitoring = false;                             % stop data monitoring
                    end
                    
                    % display the estimated powers
                    StagePowers
                end

                % analyse data
                T_k(k)  = Calculate_F(Data);    % the stage k test statistic
                S_k(k)  = sum(T_k);             % the stage k summary statistic
                
                % update all accrued data
                AllData = [AllData;Data];

                % evaluate stop conditions
                if S_k(k) >= A_k              	% reject H0
                    ContinueTesting	= false;
                end
                if k == K                       % final stage reached 
                    ContinueTesting = false;             
                end

                % prepare for the next stage
                if ContinueTesting
                    NullTrunc                       = Null;
                    NullTrunc(Efficacy_Ind(1):end)	= zeros(1, length(NullTrunc(Efficacy_Ind(1):end)));	% all values larger than the efficacy threshold are set to zero
                    AltTrunc                        = Alt;
                    AltTrunc(Efficacy_Ind(1):end)	= zeros(1, length(AltTrunc(Efficacy_Ind(1):end)));	% all values larger than the efficacy threshold are set to zero
                end

            end

            % store results here
            % ...

        end
    end
end



