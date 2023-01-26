function [params, err] = Find_Pol(x_obs, f_obs, trsh)
% 
% Find_Pol find the best polynomial interpolation for the given observation
% 
% Input:
% * x_obs = observation abscissa
% * f_obs = observed values
% * trsh = treshold
%
% Output:
% * params = numbers of parameters of the found polynomial
% * err    = residual error between the polynomials of k-1 over k params

if nargin < 3
   trsh = 1;      % proportion of rows to select for training
end

err = sum((mean(f_obs)-f_obs).^2); 

for params = 2:40 % iterate on the numbers of parameters                 

    F = Pol_Interp(x_obs, params); % Polynomial Interpolation of Degree: params-1

    a_hat = LS(F,f_obs);           % compute LS opt solution
    
    err(params) = sum(((F*a_hat)-f_obs).^2);

    if err(params-1)/err(params) - 1 < trsh
        err = err(2:end);
        break
    end
end

end