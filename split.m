function [train,test] = split(dataA, p)

if nargin < 2
   p = .9;      % proportion of rows to select for training
end
dataA = dataA';
N = size(dataA,1);  % total number of rows 
tf = false(N,1);    % create logical index vector
tf(1:round(p*N)) = true;     
tf = tf(randperm(N));   % randomise order
train = dataA(tf,:)';
test = dataA(~tf,:)';

end