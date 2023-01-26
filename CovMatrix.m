function tau = CovMatrix(x_obs,x_synth)
% 
% CovMatrix build the covariance matrix with the distances between the
% observed abscissa and the interpolating abscissa
% 
% Input:
% * x_obs   = observation abscissa
% * x_synth = interpolating values
% 
% Output:
% * tau = matrix of the distances

if nargin < 2
    x_synth = x_obs;
end

x_obs_std = STD(x_obs);

x_p_std = STD(x_synth);

[t1_p,t2_p] = meshgrid(x_p_std,x_obs_std);
tau = abs(t1_p-t2_p);

end