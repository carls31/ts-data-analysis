# Geospatial Data Analysis - lab
Graduate Project Assignment of the course "Geospatial Data Analysis" at Politecnico di Milano

The main source code is avaiable [here](https://github.com/carls31/GDA-Lab/blob/main/LAB_Assignment.m) and, in order to compare the results, the same analysis was conducted on two additional boreholds: [LAB_515](https://github.com/carls31/GDA-Lab/blob/main/LAB_515.m) and [LAB_518](https://github.com/carls31/GDA-Lab/blob/main/LAB_518.m)

The provided time series are real data made available by ARPA Veneto, the regional agency for environment protection. The whole dataset
is composed by measurements from several boreholes each identified by a numeric ID
The missing values are removed from the data table and they will be replaced if needed; the features used are:
 - epochs 'DATA' converted upstream in integer numbers;
 - piezometer measurements 'LIVELLOSTATICO [mslm]' taken between 1999 and 2021 (roughly);
 - boreholds ID planar cartographic coordinates in the district of Padua;
 - altitude of the boreholds ID 'QUOTA P.R. [mslm]';
 - interpolated equidistant epochs;
 
 Find the best fitting polynomial with a low degree using a treshold on the ratio between consecutive polynomials of dregree k and k+1 until
the improvement is too low.
Once obtained the best number of parameters for the Smoothing Least Square Interpolation, the polynomial model is build computing the
optimum solution that minimize the residuals.
The polynomial describes the overall overview of the samples behaviour and, assuming this as the deterministic component of the given
samples, it will be deployed later on the further analysis.

Each timeseries covers a different period and although there are 4 times per year the dates of acquisition are different. In order to get an
equally distant interval of sampling and fill the holes provided by the missing values, Exact Cubib Spline is a suitable interpolation
method based on a polynomial of third order in the Newton form.

Assuming the provided piezometer measurements to be samples of a signal that can be decomposed in a finite sum of sin and cosine
functions.
The resulting equidistant sampling got by the Exact Spline is required to build the power spectrum of the signal.
The power spectrum shows that there is no clear evidence to say if the frequencies are noise or harmonic component for example at 6
months.

The developped script can be more accurate in some boreholes than other, due to the variety in the samples distribution between different
time series. Considering clustering them could help Including more information about the data may be helpful, for example adding some
layers could help find a stronger pattern on the behaviour.
A proposal for further analysis may be the implementation of Linearized LS on the hyperparameters and also the implementation of Kriging
method. Furthermore, a comparison of the performance of different technique or different boreholes with evaluation metrics (eg R-squared,
R-squared adj, AIC, BIC)
