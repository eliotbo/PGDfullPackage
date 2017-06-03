%This function uses the diluted iterative algorithm for full tomography.
%This method is often the fastest for ill-conditioned measurement matrices,
%but requires optimisation of the step size 'gamma' and the momentum 
%paramter 'alpha'. In addition to this, the freeing the trace makes the
%algorithm more likely to fall into a local minimum. 
%If the step size is too high, numerical instabilities ensue.
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%output
%rhoPGDM: maximum likelihood density matrix

function [rhoPGDM, timePGDM, costfunction] = runPGDMfreeTrace(data,A)
initialTic = tic;

sA = size(A);
d = sA(2);
N = sA(1);

%initial hyperparameters
condA = cond(A);

gam = 1;
alpha = 0.98;

%initial guessis identity
rhok = eye(d)*mean(data);

%initialisation
it=0;
vr = rhok*0;
orderMagnitudeCost = 1e100;

quit = 0;
while quit == 0
    it=it+1;
    
    %projection on the simplex
    [rhok] = positiveProjection(rhok);
    
    %Calculate probabilities
    temp1 = conj(A).*((A)*rhok);
    Ax = sum(temp1,2);
    clear('temp1');
    
    Axb = Ax-data;
    
    %log L_{GP}
    costFunction(it) = sum( ((Ax-data)./sqrt(data+eps))  .^2)/N;
    
    %increase the value of alpha with k
    if ceil(log10(costFunction(it)))< orderMagnitudeCost && condA>2
        orderMagnitudeCost = ceil(log10(costFunction(it)));
        beta = 1-condA/200; if beta<0.5, beta = 0.5; end
        alpha = (1-(1-alpha)*beta);
    end
    
    %Change estimate with momentum
    grad = (bsxfun(@times,A,Axb./sqrt(data))'*A);
    vr = vr*alpha - gam*grad ;
    rhok = rhok + vr;
    
    %exit criterion
    if it>100 && costFunction(end)<1.2 && mean(abs(diff(costFunction(end-20:end))))<1e-4,
        quit=1; end
    if it > 2000, quit = 1; end
    
    %detect numerical instabilities
    if costFunction(it)*1E-3 > costFunction(1)
        rhok = eye(d)/d; disp('Numerical instability with PGDM');
        costFunction=[]; costFunction=0; break;
    end
end

%plot log L_{GP}
figure(222);semilogy(abs(costFunction),'b--');drawnow
xlabel('Iteration number')
ylabel('log L_{PG}')

timePGDM = toc(initialTic);
[rhok] = positiveProjection(rhok);
rhok = rhok/trace(rhok);
rhoPGDM = rhok;
