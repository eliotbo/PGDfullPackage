 %This function uses the PGDM for full tomography. This method is often the
%fastest for ill-conditioned measurement matrices, but requires
%optimisation of the step size 'gamma' and the momentum paramter 'alpha'.
%If the step size is too high, numerical instabilities ensue.
%
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%r : multiplication factor that turns probabilities into estimated clicks
%output
%rhoPGDM: maximum likelihood density matrix

function [rhoPGDM, timePGDM, costFunctionReduced] = runPGDM(data,A,r)
% disp('PGDM')
initialTic = tic;

sA = size(A);
d = sA(2);
N = sA(1);

%trick to be able to solve high-purity tomographic problems: set pp to 0.9
pp=1;
pDataFake = data*pp + (1-pp)*ones(length(data),1)*r/d;
data = pDataFake;

%initial hyperparameters
alpha = 0.95;
gam = 1/d/r/2;
% gam = .1/d^2;

%initial guessis identity
rhok = eye(d)/d;

%initialisation
it=0;
vr = rhok*0;

costFunction = zeros(1,10000);
quit = 0;
while quit == 0
    it=it+1;
%     disp(['iteration: ' num2str(it)])
    %projection on the simplex
    [rhok] = simplexProj(rhok);
    
    %Calculate probabilities
    temp1 = conj(A).*((A)*rhok);
    Ax = sum(temp1,2)*r;
    clear('temp1');
    
    Axb = Ax-data;
    
    %log L_{GP}
    costFunction(it) = sum( ((Ax-data)./sqrt(data+eps))  .^2)/N;

%%%%%%%%%%%%%%% Experimenting with the momentum value %%%%%%%%%%%%%%%%%
%     %increase the value of alpha with k
%     if ceil(log10(costFunction(it)))< orderMagnitudeCost 
%         orderMagnitudeCost = ceil(log10(costFunction(it)));
%         alpha = (1-(1-alpha)*0.95);
%     end
%%%%%%%%%%%%%%% Experimenting with the momentum value %%%%%%%%%%%%%%%%%   

%%%%%%%%%%%%%%% Alternative cost funcion and gradient %%%%%%%%%%%%%%%%%
    %   % non-convex version of the Gaussian-Poisson approximation gradient
    %     tempg = ((Axb./Ax).*(2-(Axb./Ax)));
    %     grad = (bsxfun(@times,A,tempg)'*A);
%%%%%%%%%%%%%%% Alternative cost funcion and gradient %%%%%%%%%%%%%%%%%

    %Change estimate with momentum
    grad = (bsxfun(@times,A,Axb)'*A);
    vr = vr*alpha - gam*grad ;
    rhok = rhok + vr;
    
    %exit criterion
    if it>100 && costFunction(end)<2 && mean(abs(diff(costFunction(it-20:it))))<1e-5
        quit=1; end
    if it > 2000, quit = 1; end
    
    %detect numerical instabilities
    if costFunction(it)*1E-2 > costFunction(1)
        rhok = eye(d)/d; disp('Numerical instability with PGDM');
        costFunction=[]; costFunction=0; break;
    end
    
end

%high-purity trick finish
rhok  = (rhok - (1-pp)*eye(d)/d)/pp;

costFunction = abs(costFunction);
costFunctionReduced = costFunction(1:it);
timePGDM = toc(initialTic);
[rhok] = simplexProj(rhok);
rhok = rhok/trace(rhok);

rhoPGDM = rhok;
