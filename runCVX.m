%this function uses CVX to solve a full tomography problem or a compressive
%sensing tomography problem. Change the cost function below if the density
%matrix is known to be quasi pure.
%
%make sure CVX is installed on your computer
%
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%r (optional): multiplication factor that turns probabilities into estimated clicks
%output
%rhoCVX: maximum likelihood density matrix


function [rhoCVX, timeCVX]=runCVX(data,A,r)

sA = size(A);
d = sA(2);
N = sA(1);
conditionNumber = cond(A);

%The likelihood L_{PG} can't have a data vector with zeros, either add a
%small constant to the the zeros or change the likelihood used below
dataContainsZeros=0;
for i=1:N, if data(i)==0, dataContainsZeros=1; end; end

tic

cvx_begin quiet
    cvx_solver sdpt3
    cvx_precision low
    variable rho_estimate(d,d) complex semidefinite

    if ~exist('r')
        %r is included in rho_estimate
        r=1;
    else
        %enforce unit trace (r is not included in rho_estimate)
        trace(rho_estimate)==1;
    end
    
    %calculate estimated detector clicks
    temp1 = conj(A).*((A)*rho_estimate);
    Ax = sum(temp1,2)*r;

    %%%% if the density matrix is known to be quasi pure %%%%%%
    % minimize 0.008*norm_nuc(rho_estimate)+(norm((Ax  - data),'fro')/N) ;
    %%%% if the density matrix is known to be quasi pure %%%%%%
    
    
    if dataContainsZeros==1 
        minimize norm((Ax  - data),'fro')/sqrt(N);
    else
        minimize norm((Ax  - data),'fro')/sqrt(N);
%         minimize norm((Ax  - data)./sqrt(data+eps),'fro')/sqrt(N);
    end
cvx_end

timeCVX  = toc;
rho_estimate = rho_estimate/trace(rho_estimate);
rhoCVX = rho_estimate;