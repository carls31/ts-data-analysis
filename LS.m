function a_hat = LS(F,f_obs)
% 
% Least Squares for Linear Interpolation
%
%  Input: 
%  -> F     = design matrix
%  -> f_obs = target values
%
%  Output:
%  -> a_hat = optimal LS coefficents vector for the given F
% 

N = (F'*F) + (2^(-50))*eye(size(F,2));

% [U,S,V] = svd(N);
% N_new = U*(S + (2^(-50))*eye(size(F,2)))*V';

u = F'*f_obs;
a_hat = N\u;

end