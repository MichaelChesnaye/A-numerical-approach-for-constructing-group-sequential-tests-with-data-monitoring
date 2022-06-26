%%
clear all
close all
clc

% test parameters
K_Val               = [1,2,4,8]; 	% the number of stages for the sequential test. Test performance was evaluated for different K
Alpha               = 0.025;      	% the type-I error rate for the full K-staged sequential test
Gamma               = 0.9;       	% the statistical power for the full K-staged sequential test
ES_Assumed          = 1;         	% the assumed simulated effect size
Na_Val              = [10, 20, 50];     % the number of accrued samples before data monitoring is initiated. Test performance was evaluated for different "Na".
Sigma2_Val          = 1:50;             % the sample variance for the simulated data. Test performance was evaluated for different "Sigma2"
Q                   = 1;                % the number of features for the simulated data set. The obne-sample t-test is used for Q=1
Resolution          = (1/0.0005);        	% the resolution of the axis along which the distributions are generated
XAxis_Values        = 0:1/Resolution:60;	% x-axis values along which the distributions are generated

% the following simulations were run 200,000 times using the Iridis supercomputer at Southampton University. 
for i0=1:length(Na_Val)
	for i1=1:length(Sigma2_Val)
        for i2=1:length(K_Val)
        
            % unpack
            Na          = Na_Val(i0);           % the number of accrued samples before data monitoring is initiated
            Sigma2      = Sigma2_Val(i1);       % the variance of the data to be simulated
            K           = K_Val(i2);            % the number of stages for the sequential analysis
            Alpha_k    	= ones(1,K) * (Alpha / K);                      % the stage-wise type-I error rates
            Gamma_k     = Get_Gamma_k( K_Val(i2), Q, Alpha, Gamma );    % The stage-wise statistical powers (cumulative across stages)
            
         	% initialise 
            k               = 0;        % the curernt stage of the analysis
            AllData1        = [];       % all accrued data for group 1, up to stage k
            AllData2        = [];       % all accrued data for group 2, up to stage k
            T_k             = [];       % the stage-wise test statistics
            StagePowers     = [];       % the estimated stage-wise statistical powers
            ContinueTesting	= true;
            while ContinueTesting       % sequential test starts here

                k           = k + 1; 	% increment stage index    
                Data_1      = [];       % group 1 data for stage k
                Data_2      = [];       % group 2 data for stage k
                N_k         = 0;        % sample size for stage k is 0
                Monitoring  = true;     % monitor sample variance boolean
                AddN        = 1;        % monitor the sample variance after every "AddN" new samples, starting with 1 (this will be increased as the trial progresses)
                    
                while Monitoring    	% continuously monitor sample variance 

                    % generate data
                    if N_k == 0
                        Data_1	= sqrt(Sigma2) * randn(Na,Q);	% add Na observations for group 1
                        Data_2	= sqrt(Sigma2) * randn(Na,Q);   % add Na observations for group 2
                    else
                        Data_1 	= [Data_1; sqrt(Sigma2) * randn(AddN,Q)];	% add AddN observations for group 1
                        Data_2 	= [Data_2; sqrt(Sigma2) * randn(AddN,Q)];   % add AddN observations for group 2
                    end
                    N_k = 2*size(Data_1,1);

                    % increment AddN
                    if N_k > 50
                        AddN = ceil( N_k * 0.05 );	% as the trial progresses, data is monitored less often to reduce computation time
                    end

                    % generate null distribution
                    Null = fpdf(XAxis_Values, Q, N_k-Q-1);      % the null distribution for stage k 
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
                    Cov1        = cov( [AllData1;Data_1] );     % group 1 covariance matrix
                    Cov2        = cov( [AllData2;Data_2] );     % group 2 covariance matrix
                    Cov_Pooled  = (Cov1 + Cov2) / 2;          	% pooled covariance matrix under the assumption of equal group sizes, else Cov_Pooled  = ((n1-1)*Cov1 + (n2-1)*Cov2) / (n1+n2-2);
                    InvCov    	= inv( Cov_Pooled );            % the inverse
                    NC_estimate = N_k * ES_Assumed * InvCov * ES_Assumed';          % the estimated non-centrality parameter     
                    Alt       	= ncfpdf(XAxis_Values, Q, N_k-Q-1, NC_estimate );   % and the alternative distribution under the assumed effect size
                    if Alt(1)==inf
                        Alt(1) = Alt(2) + (Alt(2)-Alt(3));  % avoid infinities
                    end
                    Alt = Alt / sum(Alt);   % normalise so the area equals 1
                    if k>1                                          % for stage two onwards... 
                        Alt = conv( Alt, AltTrunc );                % convolve with the truncated null distribution from the previous stage
                        Alt = (Alt ./ sum(Alt)) * sum( AltTrunc );  % and normalise the area so that it equals 1 minus the area of the trunacted regions
                    end

                    % estimate statistical power
                    StagePowers(k) = sum( Alt(Efficacy_Ind(1):end) );   % estimated power (cumulative across stages)
                    if sum( StagePowers ) >= Gamma_k(k)                 % if the power exceeds the pre-specified power level
                        Monitoring = false;                             % stop data monitoring
                    end
                    
                    % display estimated powers
                    StagePowers
                end
                
                % analyse data using the two-sample t-test
                n1      = size(Data_1,1);   % sample size group 1
                n2      = size(Data_2,1);   % sample size group 2
                M1      = mean(Data_1);     % group 1 mean
                M2      = mean(Data_2);     % group 2 mean
                Cov1    = cov(Data_1);      % group 1 variance
                Cov2    = cov(Data_2);      % group 2 variance
                Cov_P   = ((n1-1)*Cov1 + (n2-1)*Cov2) / (n1+n2-2);              % pooled variance
                T2      = ((n1*n2)/(n1+n2)) * (M1-M2) * inv(Cov_P) * (M1-M2)';  % squared t statistic (i.e. the T2 statistic of the Hotelling's T2 test). 
                T_k(k)	= ( (n1+n2-Q-1) / ((n1+n2-2)*Q) ) * T2;  	% the stage k test statistic
                S_k(k)  = sum(T_k);                                 % the stage k summary statistic
                
                % update all accrued data
                AllData1 = [AllData1;Data_1];
                AllData2 = [AllData2;Data_2];
                    
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



