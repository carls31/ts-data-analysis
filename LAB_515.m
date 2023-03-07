%% GDA - Lab Assignment
% Student: LORENZO CARLASSARA 101724
% emma version 3.0
%% Data Preprocessing
close all, clear all, load PZ.mat 
% help Find_Pol
% help LS
% help Pol_Interp
% help Exact_Cubic_Spline
% help Cubic_Spline
% help Cubic_Spline_2D
% help CovMatrix
% help CovFunction
id = 25; 
%for id = 14:18%length(an_)
x_obs = pz.dataEpoch(pz.N_POZZO == an_(id));
f_obs = pz.LIVELLOSTATICO_MS_l_m_(pz.N_POZZO == an_(id));
x_syn = (x_obs(1):x_obs(end))';
h_obs = (x_obs(1)+1:91.25:x_obs(end))';
x_dt = datetime(86400*(x_obs-1),'ConvertFrom','epochtime','Epoch','1999-01-01','Format','dd/MMM/yyyy');
h_dt = datetime(86400*(h_obs-1),'ConvertFrom','epochtime','Epoch','1999-01-01','Format','dd/MMM/yyyy');
s_dt = datetime(86400*(x_syn-1),'ConvertFrom','epochtime','Epoch','1999-01-01','Format','dd/MMM/yyyy');
n_pz = sprintf(' ID %d',an_(id));
%% Smoothing Least Square Interpolation

trsh = 0.0001; 
[best, err] = Find_Pol(x_obs, f_obs, trsh);

params = best; 
F = Pol_Interp(x_obs, params); 
a_hat = LS(F,f_obs);  
F_interp = Pol_Interp(x_syn,params);
f_LS = F_interp*a_hat;

figure
hold on
deg = sprintf('polynomial degree: %d',params-1);
title("LS Polynomial Interpolation" + n_pz)
plot(s_dt,f_LS,'Color',[0.8500 0.3250 0.0980]);
plot(x_dt,f_obs,'o','Color',[0 0.4470 0.7410]); grid on;
xlabel('years'), ylabel('[m]'), legend(deg,"measurements",'Location','best')
hold off
%% Exact Cubic Spline

[f_h, x_h] = Exact_Cubic_Spline(x_obs,f_obs,h_obs);

figure, plot(h_dt,f_h,'-*','Color',[0.8500 0.3250 0.0980]) 
title("Exact Interpolation" + n_pz)
hold on, plot(x_dt,f_obs,'o','Color',[0 0.4470 0.7410]); grid on;
xlabel('years'), ylabel('[m]'),legend('equidistant cubic spline',"measurements",'Location','best')

%% Discrete Fourier Transformation
% Since the results provided by the Exact Cubic Spline are misleading 
% (overfitted)
% the DFT can not be taken in consideration

N = length(h_obs);
if (mod(N,2)==0) %even   
    f = 1/N * (-N/2 : N/2-1)'; 				% intervall [-fs/2 , fs/2]
else             %odd
	f = 1/N * (-(N-1)/2 : (N-1)/2)'; 		% intervall [-fs/2 , fs/2]
end
Fobs1 = 1/N * fftshift(fft(f_h));	% Fourier transformation

Fobs = (abs(Fobs1) - min(abs(Fobs1)))/(max(abs(Fobs1))-min(abs(Fobs1)));
figure,
plot(f,abs(Fobs1),'x-','Color',[0 0.4470 0.7410])
title('Fourier transform (power spectra)')
xlabel('frequency'),grid on
ylim([0,1])
deg = sprintf('obs signal - a0 = %g',max(abs(Fobs1)));
legend(deg,'Location','best')

%% Collocation Technique 

% Compute obs lag matrix
[t1,t2] = meshgrid(x_obs,x_obs);
TAO = abs(t1-t2);
TAO = triu(TAO,1);

% Assuming that the polynomial interpolation catches all the deterministic 
% component that characterizes the Random Process
f_LS = Pol_Interp(x_obs, params)*a_hat;
f_obs_nomean = f_obs - f_LS;
figure, plot(x_dt,f_obs_nomean), title("Residuals" + n_pz)
xlabel('years'), ylabel('[m]')

% Compute correspondent obs cov matrix
[p1,p2] = meshgrid(f_obs_nomean,f_obs_nomean);
C = p1.*p2;
C = triu(C,1);

%Reshape lag and cov into arrays 
tao = TAO(:);
c = C(:);

% Consider single combinations only
c(tao==0)=[];
tao(tao==0)=[];

% Empirical covariance (ECF) computation
delta_tao_ecf =365;		                       %ECF sampling interval: 1 year
max_tao_ecf = max(x_obs)/2;			               %farther ECF sample: half of the time span
tao_ecf = 0:delta_tao_ecf:max_tao_ecf;			    %ECF tao bin edges
tao_inbin = discretize(tao,tao_ecf);					%ECF tao binning
ECF=[];
for bin = 1:length(tao_ecf)-1							%ECF tao binning
	c_inbin = c(tao_inbin==bin);
	ECF(bin,1) = tao_ecf(bin)+ delta_tao_ecf/2;
	ECF(bin,2) = mean(c_inbin);
end
ECF = [[0 var(f_obs_nomean)]; ECF];						%add C(0) to ECF

figure
plot(ECF(:,1),ECF(:,2),'*b'), grid on, hold on
xlabel('days')
delta_tao_ecf =365/2;		                       %ECF sampling interval: 1 year
max_tao_ecf = max(x_obs)/2;			               %farther ECF sample: half of the time span
tao_ecf = 0:delta_tao_ecf:max_tao_ecf;			    %ECF tao bin edges
tao_inbin = discretize(tao,tao_ecf);					%ECF tao binning
ECF=[];
for bin = 1:length(tao_ecf)-1							%ECF tao binning
	c_inbin = c(tao_inbin==bin);
	ECF(bin,1) = tao_ecf(bin)+ delta_tao_ecf/2;
	ECF(bin,2) = mean(c_inbin);
end
ECF = [[0 var(f_obs_nomean)]; ECF];		

plot(ECF(:,1),ECF(:,2),'*r'),title("Empirical Covariance Function" + n_pz)
xlabel('distance [days]')

% Exact Polynomial interpolation
x_emp = ECF(2:4,1);
f_emp = ECF(2:4,2);
K = zeros(length(x_emp),length(x_emp)); % initialize design matrix
for i = 0:length(x_emp)-1
    K(:,i+1) = (x_emp).^i; % compute design matrix
end
u = K'*f_emp;
k_hat = inv(K'*K)*u;

A = k_hat(1); 
alpha = - log(ECF(2,2)/A)/2;
%end

%% Kernel method for linear regression interpolation
% Compare the result provided by the ECF with a priori assumption on the
% model distribution 
% Consider to merge the results provided by the Exact Cubic Spline and the
% GP Reggression

% find the hyperparameters of the kernel function using mle
model = "Normal";
phat = mle(f_obs, 'Distribution', model); 
sigma = phat(1); len = phat(2); 

 % Matrix of the distances
tau_ss = CovMatrix(x_obs);
tau_ps = CovMatrix(x_obs,h_obs); 
tau_pp = CovMatrix(h_obs);

% Covariance matrix
Css = CovFunction(tau_ss,sigma,len,model); 
Csp = CovFunction(tau_ps,sigma,len,model);
Cpp = CovFunction(tau_pp,sigma,len,model);

eps = 2^(-50); % prevent dividing by 0
noise = 0.05;   % noise
C = (Css + (eps + noise)*eye(length(x_obs))); 

L = chol(C)'; % Cholesky Decomposition
alfa = L'\(L\f_obs); %  L*L' = C
y_L = Csp'*alfa; 

y_pred = Csp' * (C\f_obs);

figure, plot(h_dt,y_pred,'*','Color',[0.8500 0.3250 0.0980])
hold on, plot(x_dt,f_obs,'o','Color',[0 0.4470 0.7410]) 
title("Random Process" + n_pz)
legend("prior: "+ model+" model",'measurements','Location','best')