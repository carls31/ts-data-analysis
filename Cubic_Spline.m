function F = Cubic_Spline(obs,spl,h)
% 
% Cubic_Spline compute 
% 
% 
% 
% 
% 
% 
% 
% 

n_obs = length(obs);
n_spl = length(spl);


[SPL,OBS] = meshgrid(spl,obs);

F = zeros(n_obs,n_spl);
D = (OBS-SPL)./h;
for i=1:n_obs
    for j=1:n_spl
        if (abs(D(i,j)) <= 2 )
            F(i,j) = (2-abs(D(i,j)))^3/6;
        
            if ( abs(D(i,j)) <= 1)
                F(i,j) = F(i,j) - 4*(1-abs(D(i,j)))^3/6;
            end
        end
    end
end

% [SPL,OBS] = meshgrid(x_spl,x_obs);
% D = (OBS-SPL)/h;
% A1 = 1/6*(D+2).^3;
% A2 = 1/6*((D+2).^3-4*(D+1).^3);
% A3 = 1/6*((2-D).^3-4*(1-D).^3);
% A4 = 1/6*(2-D).^3;
% A1(D<=-2 | D>=-1) = 0;
% A2(D<=-1 | D>=0) = 0;
% A3(D<=0 | D>=1) = 0;
% A4(D<=1 | D>=2) = 0;
% A = A1+A2+A3+A4;
% A(D==0) = 2/3;
% A(D==-1) = 1/6;
% A(D==+1) = 1/6;
% F = A;

end