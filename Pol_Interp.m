function F = Pol_Interp(x_obs, params )
% 
% Pol_Interp build the design matrix required to comupute a polynomial
% interpolation of params-1 degree
% 
% Input:
% * x_obs  = observation abscissa
% * params = numbers of parameters of the found polynomial
% Output: 
% * F = design matrix

n_obs = length(x_obs); % numbers of observation


F = zeros(n_obs,params); % initialize design matrix
for i = 0:params-1
    F(:,i+1) = (x_obs).^i; % compute design matrix
end
% smooth_y = smoothdata(y)
end