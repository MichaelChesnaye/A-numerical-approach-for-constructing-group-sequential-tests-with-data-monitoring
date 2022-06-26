%%
function Gamma_k = Get_Gamma_k(K, Q, Alpha, Gamma)

% generate the power function
ES          = 0.001;                                                % a small effect size to obtain smooth contour of the power function
N           = (1+Q):60000;                                      	% sample size values: the axis along which the power function is generated
NonCen      = ES*N;                                                 % all corresponding non-centrality parameters
Pow         = 1 - ncfcdf(finv(1-Alpha, Q, N-Q), Q, N-Q, NonCen);    % HT2 power as a function of sample size

% locate Gamma
AMin        = abs( Pow - Gamma );               
MinIndex	= find( AMin==min(AMin) );          % Gamma location along the power function
n_k         = floor( N(MinIndex(1)) / K );    	% split the axis in K parts 
Gamma_k     = Pow( floor( n_k:n_k:(n_k*K) ) );  % and find the correponding power levels 




% 
% % generate the power function
% ES          = 0.001;                                                % a small effect size to obtain smooth contour of the power function
% N           = (1+Q):60000;                                      	% sample size values: the axis along which the power function is generated
% NonCen      = ES*N;                                                 % all corresponding non-centrality parameters
% Pow         = 1 - ncfcdf(finv(1-Alpha, Q, N-Q), Q, N-Q, NonCen);    % HT2 power as a function of sample size
% 
% % locate Gamma
% AMin        = abs( Pow - Gamma );               % 
% MinIndex	= find( AMin==min(AMin) );          % Gamma location along the power function
% Pow_Trunc	= Pow( 1:MinIndex(1) );             % The power function up to Gamma
% N_Trunc     = N( 1:MinIndex(1) );               % and the corresponding axis along which the (now trunacated) power function is defined
% N_k         = N_Trunc(end) / K;                 % take ~equidistance steps along the axis
% for k=1:K
%     AMin        = abs( N_Trunc - (N_k*k) );     % 
%     MinIndex	= find( AMin==min(AMin) );      % 
%     TPRs(k)     = Pow_Trunc( MinIndex(1) );     % 
% end
% TPRs(end) = Gamma;                              % 
% 
% 
% 
% 
% % generate the power function
% ES          = 0.001;                                                % a small effect size to obtain smooth contour of the power function
% N           = (1+Q):60000;                                      	% sample size values: the axis along which the power function is generated
% NonCen      = ES*N;                                                 % all corresponding non-centrality parameters
% Pow         = 1 - ncfcdf(finv(1-Alpha, Q, N-Q), Q, N-Q, NonCen);    % HT2 power as a function of sample size
% 
% % locate Gamma
% AMin        = abs( Pow - Gamma );               % 
% MinIndex	= find( AMin==min(AMin) );          % Gamma location along the power function
% n_k         = floor( N(MinIndex(1)) / K );               % split the axis in K parts 
% Gamma_k     = Pow( floor( n_k:n_k:(n_k*K) ) );  % and find the correponding power levels 
% 
% 
% 
% Pow_Trunc	= Pow( 1:MinIndex(1) );             % The power function up to Gamma
% N_Trunc     = N( 1:MinIndex(1) );               % and the corresponding axis along which the (now trunacated) power function is defined
% N_k         = N_Trunc(end) / K;                 % take ~equidistance steps along the axis
% for k=1:K
%     AMin        = abs( N_Trunc - (N_k*k) );     % 
%     MinIndex	= find( AMin==min(AMin) );      % 
%     TPRs(k)     = Pow_Trunc( MinIndex(1) );     % 
% end
% TPRs(end) = Gamma;                              % 

