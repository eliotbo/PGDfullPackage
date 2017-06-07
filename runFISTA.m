%This function uses the FISTA for full tomography.
%
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%r : multiplication factor that turns probabilities into estimated clicks
%output
%rhoFISTA: maximum likelihood density matrix
%timeFISTA: time taken by the algorithm
%costFunction: cost function value as a function of iteration number

function  [rhoFISTA, timeFISTA, costFunction]=runFISTA(data,A,r)
tic
% disp('FISTA')
sA = size(A);
d = sA(2);
N = sA(1);

pp=1;
dataFake = data * pp + (1-pp)*ones(length(data),1)*r/d;

data2 = dataFake;
data = dataFake/r;

gamAcc = 0.1/d;

guess_rho = eye(d)/d;
past_guess_rho = guess_rho;
stoppingCondition=0;
costFunction = [0 0];

k=0;
tic
while stoppingCondition == 0
    k=k+1;
    y = guess_rho+(k-2)/(k+1)*(guess_rho-past_guess_rho);
    past_guess_rho = guess_rho;
    temp = conj(A).*((A)*y);
    Av = real(sum(temp,2));
    Avb =  Av - data;
    
    costFunction(k) = sum( ((Av*r-data2)./sqrt(Av*r+eps)).^2)/N;
    gradient = (bsxfun(@times,A,Avb)'*A);
  
    guess_rho = y - gamAcc*gradient;
    [guess_rho] = simplexProj(guess_rho);

    if k>100 && costFunction(end)<2 && mean(abs(diff(costFunction(end-50:end))))<1e-4
        stoppingCondition = 1; end
    if k > 2000
        stoppingCondition = 1;
        
    end
end

 guess_rho  = (guess_rho - (1-pp)*eye(d)/d)/pp;

timeFISTA = toc;
costFunction = abs(costFunction);
guess_rho = simplexProj(guess_rho);
guess_rho = guess_rho/trace(guess_rho);
rhoFISTA = guess_rho;



    