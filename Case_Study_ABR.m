%% Section 5. A case-study for ABR detection
clear all
close all
clc

load CaseStudy_Data	% load the data and the assumed effect sizes
% Data from: Lv J., Simpson D.M., & Bell S.L. “Objective detection of evoked potentials using a bootstrap technique”  Med. Eng. Phys., vol. 29, no. 2, pp. 191–198, 2007. doi: 10.1016/j.medengphy.2006.03.001

% The effect sizes
All_ES{1}   = SmallTemplate;        % small effect size
All_ES{2}   = MediumTemplate;       % medium effect size
All_ES{3}   = BigTemplate;          % large effect size

% Sequential test parameters
K           = 3;                	% number of stages for the analysis
Alpha       = 0.05;                 % Type-I error rate for the full K-staged sequential test  
Beta        = 0.05;                 % Type-II error rate for the full K-staged sequential test	
Gamma       = 0.95;                 % Statistical power for the full K-staged sequential test   
Gamma_max	= 0.95;                 % The required statistical power for removing one of the assumed effect sizes 
Alpha_k  	= ones(1,K) * (Alpha / K);      % The stage-wise type-I error rates
Beta_k      = ones(1,K) * (Beta / K);       % The stage-wise type-II error rates
Gamma_k   	= [0.5, 0.85, 0.95];            % The desired power levels following each stage of the analysis
Nb          = 50;                           % The data-monitoring interval in samples
Resolution  = 1000;                         % the F-axis resolution for the null and alternative distributions
FValues     = 0:(1/Resolution):50;          % the F-axis along which the null and alternative distributions are defined
Q           = 25;                           % the number of features of the multivariate data set to-be-analysed by the Hotelling's T2 test. 

% initialise "simulation"
All_Collected_Data  = [];           % all accured data up to this point
T_k                 = [];           % stage-wise test statistics 
N_k              	= zeros(1, K);	% stage-wise sample sizes 
A_k                 = [];           % the stage-wise critical thresholds for rejecting the null hypothesis
B_k                 = [];           % the stage-wise critical thresholds for rejecting the alternative hypothesis

% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % STAGE ONE % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % %
Estimated_Powers    = zeros(1, 3);	% stage-wise estimated statistical powers for each assumed effect size (3 total)
k                   = 1;            % current stage of the sequential analysis
Data_k              = [];           % Data collected at this stage
Monitoring          = true;         % continue collecting data until the estimated statistical power exceeds gamma_k
while Monitoring    	

    N_k(k)	= N_k(k) + Nb;          % increment the sample size
    Data_k	= Data(1:N_k(k),:);     % and increase the sample 

    % 1) generate the null distribution 
    Null = fpdf(FValues, Q, N_k(k)-Q);      % For the Hotelling's T2 test, the null distribution is given by a central F distribution
    Null = Null/sum(Null);                  % normalise so that the area equals 1
    
    % 2) generate the efficacy threshold for rejecting the null hypothesis
    Percentile      = sum(Null) - Alpha_k(k);                               % The percentile associated with the efficacy threshold for stage 1
    TruncInd_R(k)	= Get_Truncation_Index(Null, Percentile, Resolution);	% The location of the efficacy threshold along the null distribution
    A_k(k)          = TruncInd_R(k)/Resolution;                             % The efficacy threshold for stage 1

    % 3) Monitor sample variance and estimate statistical power for each of the asumed effect sizes
    CovMatrix   = cov( [All_Collected_Data;Data_k] );   % Sample covariance matrix
    InvCov      = inv(CovMatrix);                       % and it's inverse 
    for ESi = 1:length(All_ES)                          
        NonCen_Estimate     = N_k(k) * All_ES{ESi} * InvCov * All_ES{ESi}';     % The non-centrality parameter for effect size ESi
        Alt{ESi}            = ncfpdf(FValues, Q, N_k(k) - Q, NonCen_Estimate );	% The alternative distribution for this effect size, given by a non-central F-distribution
        Alt{ESi}            = Alt{ESi} / sum(Alt{ESi});                         % normalise so the area equals 1
        StagePowers(ESi)	= sum( Alt{ESi}(TruncInd_R(k):end) );               % the estimated statistical power for effect size ESi
    end

    % Stop data collection: the mean estimated power exceeds the pre-specified power 
    if mean( StagePowers ) >= Gamma_k(k)
        Monitoring = false;
    end

end

% 4) Generate the B_k futility threshold 
Distr_delta_m  	= Alt{1};                                                       % the alternative distribution under the smallest anticipated effect size
TruncInd_L(k)	= Get_Truncation_Index( Distr_delta_m, Beta_k(k), Resolution );	% The location of the B_k threshold along the alternative distribution
B_k(k)          = TruncInd_L(k)/Resolution;                                     % the B_k threshold for futility stopping

% 5) Analyse data
T_k(k)	= Calculate_F(Data_k);              % calculate the stage k test statistic, here given by an F statistic generated by the Hotelling's T2 test
S_k(k)	= sum(T_k);
    
% stop criteria
if ( S_k(k) >= A_k(k) ) || ( S_k(k) <= B_k(k) )  || ( k == K )
    ContinueTesting = false;        % stop testing if H0 is rejected or accepted, or if we've reached the final stage
else
    ContinueTesting = true;         % else go to stage two
end

% print stage one results
sprintf('%s%s%s%s%s\n%s%s%s%s%s%s%s\n%s%s\n%s%s\n%s%s', ...
    'Stage one critical thresholds: [', num2str(B_k(k)), ', ', num2str(A_k(k)), ']', ...
    'Stage one estimated powers: [', num2str(StagePowers(1)), ', ', num2str(StagePowers(2)),', ', num2str(StagePowers(3)), ']', ...
    'Stage one test statistic: ', num2str(T_k(k)) , ...
    'Stage one summary statistic: ', num2str(S_k(k)),...
    'Stage one sample size: ', num2str(N_k(k)))


% % % % % % % % % % % % % % % % %
% % % prepare for stage two % % % 
% % % % % % % % % % % % % % % % %
All_Collected_Data  = [All_Collected_Data;Data_k];      % Keep track of collected data 
Data                = Data( (N_k(k))+1:end,:);          % Update the available data 

% truncations for the null distribution
NullTrunc                       = Null;         
NullTrunc(TruncInd_R(k):end)	= zeros( 1, length(NullTrunc(TruncInd_R(k):end)) );	% Truncate the null distribution: anything larger than TruncInd_R(k) is set to zero
NullTrunc(1:TruncInd_L(k))      = zeros( 1, length(NullTrunc(1:TruncInd_L(k))) );   % and anything smaller than TruncInd_L(k) is set to zero

% truncations for the alternative distribution and the effect size
ES_Stage2           = {};           % the assumed effect sizes, to be take foward to stage two
StagePowers_Stage1	= [];           % and the corresponding statistical powers, achieved at stage 1
C                   = 0;            % counts the number of effect sizes to be taken foward to stage two
for ESi = 1:length(All_ES)
    if StagePowers(ESi) < Gamma_max	% if the power for this is less than Gamma_max, the effect size is taken fowards to stage two
        C                                       = C + 1;                % increment the counter
        ES_Stage2{C}                            = All_ES{ESi};          % store the effect size 
        StagePowers_Stage1(C)                   = StagePowers(ESi);     % store the corresponding power, achieved at stage 1
        AltTrunc_Stage1{C}                      = Alt{ESi};             % make a copy of the alternative distribution and...
        AltTrunc_Stage1{C}(TruncInd_R(k):end)	= zeros( 1, length(Alt{ESi}(TruncInd_R(k):end)) );	% truncate it: anything larger than TruncInd_R(k) is set to zero
        AltTrunc_Stage1{C}(1:TruncInd_L(k))  	= zeros( 1, length(Alt{ESi}(1:TruncInd_L(k))) );    % and anything smaller than TruncInd_L(k) is set to zero
    end
end


% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % STAGE TWO % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % %
clear Alt
k           = k + 1;    % increment stage index
Data_k      = [];       % Data for stage k
StagePowers = [];       % estimated statistical powers for stage 2
Monitoring  = true;
while Monitoring            

    % simulate data collection
    N_k(k)  = N_k(k) + Nb;
    Data_k	= Data(1:N_k(k),:);

    % 1) generate the null distribution 
    Null = fpdf(FValues, Q, N_k(k)-Q);               % For the Hotelling's T2 test, the null distribution is given by a central F distribution
    Null = Null/sum(Null);                           % normalise so that the area equals 1
    Null = conv(Null, NullTrunc);                    % convolvte with the truncated null distribution from stage 1
    Null = (Null ./ sum(Null)) * sum(NullTrunc); 	% normalise so the area equals 1 minus the previously trunacted percentiles
    
    % 2) generate the efficacy threshold for rejecting the null hypothesis
    Percentile      = sum(Null) - Alpha_k(k);                           	% The percentile associated with the efficacy threshold
    TruncInd_R(k)	= Get_Truncation_Index(Null, Percentile, Resolution);   % find the location of the efficacy threshold along the null distribution
    A_k(k)         	= TruncInd_R(k)/Resolution;                             % efficacy threshold

    % 3) Monitor sample variance and estimate statistical power for each of the asumed effect sizes
    CovMatrix   = cov( [All_Collected_Data;Data_k] );   % Sample covariance matrix
    InvCov      = inv(CovMatrix);                       % and it's inverse 
    for ESi = 1:length(ES_Stage2)                          
        NonCen_Estimate     = N_k(k) * ES_Stage2{ESi} * InvCov * ES_Stage2{ESi}';	% The non-centrality parameter for effect size ESi
        Alt{ESi}            = ncfpdf(FValues, Q, N_k(k) - Q, NonCen_Estimate );     % The alternative distribution for this effect size, given by a non-central F-distribution
        Alt{ESi}            = conv( Alt{ESi}, AltTrunc_Stage1{ESi} );                       % convolve with the stage 1 truncated alternative distribution
        Alt{ESi}            = ( Alt{ESi} ./ sum( Alt{ESi}) ) * sum( AltTrunc_Stage1{ESi} ); % normalise so the area equals 1 minus the area of the truncated regions
        StagePowers(ESi)	= sum( Alt{ESi}(TruncInd_R(k):end) );                   % the estimated statistical power for effect size ESi
    end
    
    % Stop data collection: the mean estimated power exceeds the pre-specified power 
    if mean(  StagePowers_Stage1 + StagePowers  ) >= Gamma_k(k)
        Monitoring = false;
    end
    
end

% 4) Generate the B_k futility threshold 
Distr_delta_m  	= Alt{1};                                                       % the alternative distribution under the smallest anticipated effect size
TruncInd_L(k)	= Get_Truncation_Index( Distr_delta_m, Beta_k(k), Resolution );	% The location of the B_k threshold along the alternative distribution
B_k(k)          = TruncInd_L(k)/Resolution;                                     % the B_k threshold for futility stopping

% 5) Analyse data
T_k(k)	= Calculate_F(Data_k);              % calculate the stage k test statistic, here given by an F statistic generated by the Hotelling's T2 test
S_k(k)	= sum(T_k);
    
% stop criteria
if ( S_k(k) >= A_k(k) ) || ( S_k(k) <= B_k(k) )  || ( k == K )
    ContinueTesting = false;        % stop testing if H0 is rejected or accepted, or if we've reached the final stage
else
    ContinueTesting = true;         % else go to stage two
end

% print stage two results
sprintf('%s%s%s%s%s\n%s%s%s%s\n%s%s\n%s%s\n%s%s', ...
    'Stage two critical thresholds: [', num2str(B_k(k)), ', ', num2str(A_k(k)), ']', ...
    'Stage two estimated powers: [', num2str(StagePowers(1)), ', ', num2str(StagePowers(2)), ...
    'Stage two test statistic: ', num2str(T_k(k)) , ...
    'Stage two summary statistic: ', num2str(S_k(k)),...
    'Stage two sample size: ', num2str(N_k(k)))



% % % % % % % % % % % % % % % % % %
% % % prepare for stage three % % % 
% % % % % % % % % % % % % % % % % %
All_Collected_Data  = [All_Collected_Data;Data_k];      % Keep track of collected data 
Data                = Data( (N_k(k))+1:end,:);          % Update the available data 

% truncations for the null distribution
NullTrunc                       = Null;         
NullTrunc(TruncInd_R(k):end)	= zeros( 1, length(NullTrunc(TruncInd_R(k):end)) );	% Truncate the null distribution: anything larger than TruncInd_R(k) is set to zero
NullTrunc(1:TruncInd_L(k))      = zeros( 1, length(NullTrunc(1:TruncInd_L(k))) );   % and anything smaller than TruncInd_L(k) is set to zero

% truncations for the alternative distribution and the effect size
ES_Stage3                   = {};           % the assumed effect sizes, to be take foward to stage three
StagePowers_Stage2          = [];           % and the corresponding statistical powers, achieved at stage 2
StagePowers_Stage1_updated	= [];           % and the corresponding statistical powers, achieved at stage 1
C                           = 0;            % count the number of effect sizes to be taken foward to stage three
for ESi = 1:length(ES_Stage2) 	% loop through all remaining effect sizes
    if ( StagePowers_Stage1(ESi) + StagePowers(ESi) ) < Gamma_max           % if the accumulated power (across stages one and two) is less than Gamma_max, then the effect size is taken fowards to stage three
        C                                       = C + 1;                    % increment the counter
        ES_Stage3{C}                            = ES_Stage2{ESi};           % store the effect size 
        StagePowers_Stage2(C)                   = StagePowers(ESi);         % store the corresponding power, achieved at stage 2
        StagePowers_Stage1_updated(C)           = StagePowers_Stage1(ESi);  % store the corresponding power, achieved at stage 1
        AltTrunc_Stage2{C}                      = Alt{ESi};                 % make a copy of the alternative distribution and...
        AltTrunc_Stage2{C}(TruncInd_R(k):end)	= zeros( 1, length(Alt{ESi}(TruncInd_R(k):end)) );	% truncate it: anything larger than TruncInd_R(k) is set to zero
        AltTrunc_Stage2{C}(1:TruncInd_L(k))  	= zeros( 1, length(Alt{ESi}(1:TruncInd_L(k))) );    % and anything smaller than TruncInd_L(k) is set to zero
    end
end
StagePowers_Stage1 = StagePowers_Stage1_updated;


% % % % % % % % % % % % % % % % % % % % % 
% % % % % % % STAGE THREE % % % % % % % %
% % % % % % % % % % % % % % % % % % % % %
clear Alt
k           = k + 1;        % increment stage index
Data_k      = [];           % Data for stage k
StagePowers = [];           % estimated statistical powers for stage three
Monitoring  = true;         
while Monitoring            

    % simulate data collection
    N_k(k)  = N_k(k) + Nb;
    Data_k	= Data(1:N_k(k),:);

    % 1) generate the null distribution 
    Null = fpdf(FValues, Q, N_k(k)-Q);              % For the Hotelling's T2 test, the null distribution is given by a central F distribution
    Null = Null/sum(Null);                          % normalise so that the area equals 1
    Null = conv(Null, NullTrunc);                   % convolve with the truncated distribution from stage two
    Null = (Null ./ sum(Null)) * sum(NullTrunc);	% normalise so the area equals 1 minus the previously trunacted percentiles
    
    % 2) generate the efficacy threshold for rejecting the null hypothesis
    Percentile      = sum(Null) - Alpha_k(k);                           	% The percentile associated with the efficacy threshold
    TruncInd_R(k)	= Get_Truncation_Index(Null, Percentile, Resolution);   % find the location of the efficacy threshold along the null distribution
    A_k(k)         	= TruncInd_R(k)/Resolution;                             % efficacy threshold

    % 3) Monitor sample variance and estimate statistical power for each of the asumed effect sizes
    CovMatrix   = cov( [All_Collected_Data;Data_k] );   % Sample covariance matrix
    InvCov      = inv(CovMatrix);                       % and it's inverse 
    for ESi = 1:length(ES_Stage3)                          
        NonCen_Estimate     = N_k(k) * ES_Stage3{ESi} * InvCov * ES_Stage3{ESi}';	% The non-centrality parameter for effect size ESi
        Alt{ESi}            = ncfpdf(FValues, Q, N_k(k) - Q, NonCen_Estimate );     % The alternative distribution for this effect size, given by a non-central F-distribution
        Alt{ESi}            = conv( Alt{ESi}, AltTrunc_Stage2{ESi} );                       % convolve with the stage 1 truncated alternative distribution
        Alt{ESi}            = ( Alt{ESi} ./ sum( Alt{ESi}) ) * sum( AltTrunc_Stage2{ESi} ); % normalise so the area equals 1 minus the area of the truncated regions
        StagePowers(ESi)	= sum( Alt{ESi}(TruncInd_R(k):end) );                   % the estimated statistical power for effect size ESi
    end
    
    % Stop data collection: the mean estimated power exceeds the pre-specified power 
    if mean(  StagePowers_Stage1 + StagePowers_Stage2 + StagePowers  ) >= Gamma_k(k)
        Monitoring = false;
    end
    
end
    
% 4) Generate the B_k futility threshold 
Distr_delta_m  	= Alt{1};                                                       % the alternative distribution under the smallest anticipated effect size
TruncInd_L(k)	= Get_Truncation_Index( Distr_delta_m, Beta_k(k), Resolution );	% The location of the B_k threshold along the alternative distribution
B_k(k)          = TruncInd_L(k)/Resolution;                                     % the B_k threshold for futility stopping

% 5) Analyse data
T_k(k)	= Calculate_F(Data_k);              % calculate the stage k test statistic, here given by an F statistic generated by the Hotelling's T2 test
S_k(k)	= sum(T_k);
    
% stop criteria
if ( S_k(k) >= A_k(k) ) || ( S_k(k) <= B_k(k) )  || ( k == K )
    ContinueTesting = false;        % stop testing if H0 is rejected or accepted, or if we've reached the final stage
else
    ContinueTesting = true;         % else go to stage two
end

% print stage three results
sprintf('%s%s%s%s%s\n%s%s\n%s%s\n%s%s\n%s%s', ...
    'Stage three critical thresholds: [', num2str(B_k(k)), ', ', num2str(A_k(k)), ']', ...
    'Stage three estimated powers: [', num2str(StagePowers(1)), ...
    'Stage three test statistic: ', num2str(T_k(k)) , ...
    'Stage three summary statistic: ', num2str(S_k(k)),...
    'Stage three sample size: ', num2str(N_k(k)))   

