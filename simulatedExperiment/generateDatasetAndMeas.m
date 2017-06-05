%This function generates a random density matrix 'rho' with a specified
%dimensionality and purity. It generates a measurement matrix with many
%options (see commented text below).
%input
%param: parameters of the gererated density matrix and measurement matrix.
%outputs:
%rho: randomly generated density matrix
%A: measurement matrix (see options below)
%data: data generated from rho and A with multinomial or poisson noise 
%entovs: (optional) entropy of measurement matrix overlap

function [rho, A, data] = generateDatasetAndMeas(param)

d = param.d;
purity = param.purity;
counts = param.counts;
theta = param.theta;

x = 1:d; %eigenvalue index
lambda=0;
purityTemp=0;
%Generate exponialy decreasing eigenvalues with the specified purity
while purityTemp<purity
    lambda = lambda + 0.001; %increase std until reach correct purity
    lam = exp(-lambda*x); %exponential distribution of eigenvalues
    lamb = lam/sum(lam);
    purityTemp = sum(lamb.^2);
end
lambi(1:d) = lamb; %eigenvalues

%Generate completely random density matrix with predefined eigenvalues
rho = makeRandomDensityMatrix(lambi);

%Generate measurement matrix. 
%Uncomment/comment to chose system and measurement type
[A] = makeQubitILLMUB(log2(d),theta); 


% calculate exact probabilities
temp0 = (conj(A).*((A)*rho));
dataExact = real(sum(temp0,2))*counts;

% if the total number of detector clicks is small enough, add multinomial
% noise, else add Gaussian-approximated Poisson noise. 
%Multinomial noise is too computationally intensive beyond a million clicks
if ceil(counts/d*length(dataExact)) < 2E5 && license('test', 'Statistics_Toolbox')
    data = mnrnd(ceil(counts/d*length(dataExact)),dataExact/sum(dataExact))';
else
    data = dataExact + randn(length(dataExact),1).*sqrt(dataExact);
    data = ceil(data);
    data(data<0)=0;
end

%calculate log likelihood of Gaussian-approximated Poisson noise. This
%should be equal to unity on average. It's very much like a chi^2 figure of
%merit. Uncomment to see the value
% log_L_GP = sum( ((dataExact-data)./sqrt(dataExact+eps))  .^2)/length(data)

if ~exist('entovs'), entovs=0; end

