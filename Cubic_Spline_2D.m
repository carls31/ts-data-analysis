function F = Cubic_Spline_2D(x_obs,x_spl,y_obs,y_spl,h_x,h_y)
%
% Bidimensional Cubic Spline


F_x = Cubic_Spline(x_obs,x_spl,h_x);
F_y = Cubic_Spline(y_obs,y_spl,h_y);

F = F_x.*F_y;

end