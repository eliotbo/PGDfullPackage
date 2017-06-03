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

function [rho, A, data, entropy] = generateDataset(param)

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
if mod(log2(d),1) == 0
    %      A = makeQubitMeasurementsSIC(log2(d),theta);
    if nargout == 4
        [A entropy] = makeQubitILLMUB(log2(d),theta);
    else
        [A] = makeQubitILLMUB(log2(d),theta); entropy = 0;
    end
    disp(['Simulation of a ' num2str(log2(d)) '-qubit system. The measurement matrix is made of MUB tensor products.'])
elseif d==25
    A = MUB5by5(1,5); 
    disp('Simulation of a two-qudit system. The measurement matrix is made of five-dimensional MUB tensor products.')
elseif d==50
    A = sic50(); %one 50-dimensional qudit and a 25-dim SIC-POVM
    disp('Simulation of a one-qudit system. The measurement matrix is made of a 50-dimensional SIC-POVM.')
elseif d==125
    A = MUB5by5by5GD(5); %three 5-dimensional qudits and 5-dim MUBs
    disp('Simulation of a three-qudit system. The measurement matrix is made of five-dimensional MUB tensor products.')
else
    A = makeRandomMeasurementsOneSystem(d,d); %d-dimensional qudit and d random bases 
    disp(['Simulation of a '  num2str(d) '-qudit system. The measurement matrix is made of ' num2str(d) ' random bases.'])
end

conditionNumber = cond(A);
disp(['The condition number of the measurement matrix is ' num2str(conditionNumber)])


%calculate exact probabilities
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

