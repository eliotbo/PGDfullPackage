function [rho,flag,iter,fe,costfunction] = pgqs(fun,rho0,tol,maxit,opt)
%
%   [rho,flag,iter] = pgqs(fun,rho0,tol,maxit,opt)
%
%   Projected Gradient method for minimization over density matrices
%
%   INPUT
%   fun         function to be minimized
%               the function should have the syntax 
%               [f,g] = myfun(rho,...) where f is the function value at rho
%               and g is the gradient
%
%   rho0        initial point
%   tol         tolerance in the stopping criteria
%   maxit       maximum number of iterations
%   opt 	opt=1: employs spectral scaling
%
%   OUTPUT
%   rho         solution found
%   flag        > 0 means convergence
%   iter        number of iterations
%   fe          counter of function evaluations
%
%   Reference:
%   D. S. Gon\c{c}alves, C. Lavor, M. A. Gomes Ruggiero,
%   A projected gradient method for optimization over density matrices
%   (submitted)
%
%   Copyright © 2015, D.S. Gon\c{c}alves, C. Lavor, M. A. Gomes Ruggiero
%   Campinas, 01/2015
%

% initialization
done=0;
iter=0;
rho=rho0;
gamma=1e-4;
%mu=1.0;
fe=0; % counter of function evaluations
M=5;
fmem = 1e+10*ones(M,1);

it=0;% main loop
while ~done
% for iterations=1:4000

    it = it+1;
%     disp(['iteration: ' num2str(it)])
    if (iter>=1) 
        g0 = g;
    end
    
    [f,g,p,cost] = feval(fun,rho);

    costfunction(it) = cost;
    normg = norm(g,'fro');
    
    if (iter==0)
        fm = f;
    end
    fe=fe+1;
    
    % classical scale factor
    mu=1.0;
    
     % spectral scale factor
     if (iter>=1 && opt==1)
        mu = min( max(trace(D*(g - g0))/(t*normD^2),1e-4), 1e+4 );
     end
    
     % projection of (rho - (1/mu)*g) onto S
     A = rho - (1/mu)*g;
    
%    % eigendecomposition
     [V,Di] = eig(A); 
    
%    % projection onto the unit simplex
     x = simplex_proj(diag(Di)); 
    
%    % resulting projection
     Pa = V*diag(x)*V';
     Pa = (Pa+Pa')/2;
    
    % search direction
    D = Pa - rho;
    normD = norm(D,'fro');

    gtd = real(trace(g*D));

    % potential reduction is too small
    if (gtd > -tol*1e-11)
        if ( (normD < 1e-8) )
          done=1;
          flag=2;
          continue
	else
% 	  fprintf('gtd: %f, normD: %f , mu: %f, normg: %f \n', gtd, normD, mu, normg);
	  done=1;
	  flag=-3;
% 	  break;
        end
    end
    
    if normD < tol
        done=1;
        flag=1;
        continue
    end
    
    t=1.0;
    
    rhon = rho + t*D;
    fn = feval(fun,rhon);
    
    % line-search
    while fn > fm + gamma*t*gtd
        t=0.5*t;
        
        % step-size is too small
        if t<1e-8
            done=1;
            flag=-2;
% 	    fprintf('gtd: %f, normD: %f normg: %f mu: %f \n', gtd, normD, normg,mu)
            break
        end
        
    if it>100 && costfunction(end)<2 && mean(abs(diff(costfunction(end-20:end))))<1e-4
        done=1; 
        disp('new exit')
    end
        
        
        rhon = rho + t*D;
        fn = feval(fun,rhon);
        fe=fe+1;
    end
    
    
    rho=rhon;
    f = fn;
    
    iter=iter+1;
    
     % classical line-search
     fm = fn;

%     % nonmonotone linesearch
%     fmem(1:M-1) = fmem(2:M);
%     fmem(M) = fn;
%     if (iter>M)
%     	fm = max(fmem);	
%     else
% 	fm = fn;
%     end
    
    % checking maximum number of iterations
    if iter>maxit
        done=1;
        flag=-1;
%         continue
    end

%      if costfunction(it)<0.4, done = 1; flag=2; end
    
%       figure(222);semilogy(abs(costfunction),'cyan');hold on;drawnow
end

