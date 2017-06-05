%This function uses the FISTA for full tomography.
%
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%r : multiplication factor that turns probabilities into estimated clicks
%output
%rhoFISTA: maximum likelihood density matrix

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

% %set step size
% condA = cond(A);
% if condA > 28
%      gamAcc = 0.08+1.22*(exp(-log10(condA)/0.6));
% elseif condA > 5
%      gamAcc = 1.-0.5601*log10(condA);
% elseif condA > 1.66
%     gamAcc = 1.45-1.1533*log10(condA);
% else 
%    gamAcc =  1.3 -  0.4543*log10(condA);
%  end
% gamAcc = gamAcc/d;

gamAcc = 0.1/d;

guess_rho = eye(d)/d;
gradient = guess_rho*0;
past_guess_rho = guess_rho;
stoppingCondition=0;
k=0;
tic
while stoppingCondition == 0
    k=k+1;
%     disp(['iteration: ' num2str(k)])
    y = guess_rho+(k-2)/(k+1)*(guess_rho-past_guess_rho);
    past_guess_rho = guess_rho;
    temp = conj(A).*((A)*y);
    Av = real(sum(temp,2));
    Avb =  Av - data;
    
    
    costFunction(k) = sum( ((Av*r-data2)./sqrt(Av*r+eps)).^2)/N;
    gradient = (bsxfun(@times,A,Avb)'*A);
  
    guess_rho = y - gamAcc*gradient;
    [guess_rho] = simplexProj(guess_rho);

    if k>100 && costFunction(end)<2 && mean(abs(diff(costFunction(end-50:end))))<1e-4, 
        stoppingCondition = 1; end
    if k > 2000, 
        stoppingCondition = 1;
        
    end
    
    if k==1000
        disp('FISTA took more than 1000 iterations')
    end
end

 guess_rho  = (guess_rho - (1-pp)*eye(d)/d)/pp;

timeFISTA = toc;
costFunction = abs(costFunction);
guess_rho = simplexProj(guess_rho);
guess_rho = guess_rho/trace(guess_rho);
rhoFISTA = guess_rho;



    