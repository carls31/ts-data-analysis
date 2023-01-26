function cov = CovFunction(tau,sigma,len,model)
% 
% CovFunction computes the model matrix with the chosen hyperparameters of
% the covariance function
% 
% Input:
% * tau   = matrix of the distances
% * sigma = signal variance parameter
% * len   = lengthscale parameter 
% * model = 'Normal' or 'Exponential'
% 
% Output:
% * cov = covariance function

if nargin < 4
    model = "Normal";
    
    if nargin < 2
        sigma = 1;
        len = 1/sqrt(2);
    end
end

if ( strcmp(model,"Normal") ) 
    cov = sigma*exp(-(tau.^2)/(2*len^2));
elseif( strcmp(model,"Exponential") )
    cov = sigma*exp(-(tau)*len);
else
    disp('choose between Normal and Exponential')
end




end