%main function for running the full tomography algorithms. CVX is not
%included with this zip file; download it at http://cvxr.com/cvx/download/
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%r (optional), multiplication factor turning probabilities into clicks
%output
%rhoEstimate: maximum likelihood density matrix
%timeTaken: time taken by each method
%costs: log-likelihood as a function of iteration number

 function [rhoEstimates, timeTaken, costs] = fullTomography(data,A,r)

warning('off','all');
%check format of measurement matrix and data vector;
sD = size(data);
sA = size(A);
if sD(1)~=sA(1), error(['Number of measurements is not equal to number'...
        ' of outcomes. Should the measurement matrix and/or data vector'...
        ' be transposed?']); 
end
if sA(1)<sA(2)*3/2 && sA(1)>sA(2)*2/3, 
    error('Measurement matrix must be made of kets'); 
end

%fix dimensionality of system
d = sA(2);

%Specify method(s): comment/uncomment to enable/disable the method.
j=0;legend222{1} = '';
if exist('r','var')
    [rhoEstimates.PGDM, timeTaken.PGDM, costs.PGDM]   = runPGDM(data,A,r);
    [rhoEstimates.FISTA, timeTaken.FISTA, costs.FISTA] = runFISTA(data,A,r);
     [rhoEstimates.SPGD, timeTaken.SPGD, costs.SPGD]   = runSPGD(data,A,r);
     if d<50, [rhoEstimates.DIA timeTaken.DIA, costs.DIA] = runDIA(data,A,r); end
    if d<50 && exist('cvx','file'),
        [rhoEstimates.CVX, timeTaken.CVX] = runCVX(data,A,r);
    end;
else
    %Only PGDMfreeTrace and CVX are written to accept data without 'r'.
    [rhoEstimates.PGDMfreeTrace, timeTaken.PGDMfreeTrace,...
        costs.PGDMfreeTrace] = runPGDMfreeTrace(data,A);
    if d<50 && exist('cvx','file'),
        [rhoEstimates.CVX, timeTaken.CVX] = runCVX(data,A);
    end;
end





