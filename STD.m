function x_std = STD(x_obs)

mean_obs = mean(x_obs);
std_obs = std(x_obs);
x_std = (x_obs - mean_obs)/std_obs;

end