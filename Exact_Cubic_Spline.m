function [S, x] = Exact_Cubic_Spline(x_obs,f_obs,x_synth)
% 
% Exact_Cubic_Spline compute an exact cubic spline interpolation 
% 
% Input:
% * x_obs   = observation abscissa
% * f_obs   = observed values
% * x_synth = interpolating values
% 
% Output:
% * S = interpolated observed values
% * x = interpolated observation abscissa

h = diff(x_obs);
A = diag(h(2:end-1),-1) + diag(2*(h(1:end-1)+h(2:end)),0) + diag(h(2:end-1),1);

d = 6*(f_obs(3:end) - f_obs(2:end-1))./h(2:end) -...
    6*(f_obs(2:end-1) - f_obs(1:end-2))./h(1:end-1);

y = d;
M = A\y;
M = [0; M; 0];

S = [];
x = [];
for j = 2:length(x_obs)
   
    x_j = x_synth(x_synth > x_obs(j-1) & x_synth <= x_obs(j)); 
    x = [x; x_j];                                                
    S = [S; f_obs(j-1) + ( (f_obs(j) - f_obs(j-1) )./h(j-1) - 2.*h(j-1).*M(j-1)/6 - h(j-1).*M(j)/6).*( x_j - x_obs(j-1) ) +...
                     ( M(j-1)/2 ).*(( x_j - x_obs(j-1) ).^2) +...
                     (( M(j) - M(j-1) )./( 6*h(j-1) ) ).*(( x_j - x_obs(j-1) ).^3)] ;

end

end