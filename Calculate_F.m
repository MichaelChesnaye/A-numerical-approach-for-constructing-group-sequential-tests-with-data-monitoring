%%
function F = Calculate_F(Data)
[N,Q]   = size(Data);                       % degrees of freedom
InvCov  = pinv(cov(Data));                  % the inverse of the covariance matrix
MeanEst	= mean(Data);                       % the Q-dimensional vector of feature means
T2      = N * MeanEst * InvCov * MeanEst';  % the T2 statistic
F       = (N-Q)/(Q*(N-1)) * T2;             % and the F statistic