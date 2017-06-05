%This function uses the diluted iterative algorithm for full tomography.
%
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%counts: multiplication factor r that turns probabilities into estimated clicks
%output
%rhoDIA: maximum likelihood density matrix

function [rhoDIA, timeDIA, costFunction] = runDIA(data,A,counts)

tic
sA = size(A);
d = sA(2);
N = sA(1);

%maximum number of iterations
maxit = 5000;
%Preicision for the exit criterion
exitPrecision = 0.75e-4;

data2 = data;

%H is identity (times a constant) if A forms a POVM
H = ((A)'*(A))^(-1/2)*sqrt(d);

%fix epsilon if there is no epsilon search in the algorithm
bestEpsilon = 0.5;

%initial estimate
s = size(A);
rhok = eye(s(2))/s(2);

it=0; fid=[0 0];fidPast=[0 0];kkk=0;
%%%%% begin algorithm %%%%%%
condition=0;
while condition==0
    it = it +1;
    
    %calculate probabilities
    temp1 = conj(A).*((A)*rhok);
    p = abs(sum(temp1,2))+eps; %estimated probabilities
    p2 = p*counts; %estimated clicks
    
    %compute gradient of log L_M and calculate R
    temp2 = (-data./p);
    grad = (bsxfun(@times,A,temp2)'*A);
    R=-grad;
    R = H*R*H/counts;
    
    %search of the best epsilon by trying four different values
    k=0;
    EPSILON = [0.05 0.1 0.5 100];
    for epsilon = EPSILON
        k=k+1;
        
        %local estimates 
        rhok2 = (eye(d)+epsilon*R)/(1+epsilon)*rhok*(eye(d)+epsilon*R)/(1+epsilon);
        rhok2 = rhok2/trace(rhok2);
        
        %calculate probabilities as a function of epsilon
        MM = abs((A)*sqrtm(rhok2)).^2;
        p = sum((MM),2);
        p2 = p*counts;
        
        %simple sum of squares as the cost function
        cost(k) =sum( ( p2-data2)  .^2);
    end
    
    %keep the epsilon that minimises the cost function
    pos = find(abs(cost) == min(abs(cost)));
    bestEpsilon = EPSILON(pos(1));
    maxEp(it) = bestEpsilon;
    %search finished
    
    rhok = (eye(d)+bestEpsilon*R)/(1+bestEpsilon)*rhok*(eye(d)+bestEpsilon*R)/(1+bestEpsilon);
    rhok = rhok/trace(rhok);
    
    %calculate probabilities
    temp1 = conj(A).*((A)*rhok);
    p = abs(sum(temp1,2))+eps;
    p2 = p*counts;
    
    %log L_{PG}
    costFunction(it) =sum( (( p2-data2)./sqrt( p2+eps))  .^2)/length(p);
    
    %exit criterion depends on the mean gradient of the cost function for
    %the last 20 iterations
    if it >40 && costFunction(it) < 2 && mean(abs(diff(costFunction(end-20:end))))...
            <exitPrecision,condition=1; end
    if it>maxit,condition=1; end
    
    if it==1000
        disp('DIA took more than 1000 iterations')
    end
end

costFunction = abs(costFunction);
rhok = rhok/trace(rhok);
rhoDIA = rhok;
timeDIA = toc;






