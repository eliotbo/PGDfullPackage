%This function uses the PGDB for full tomography.
%
%inputs
%data: Nx1 vector of integers corresponding to numbers of detector clicks
%A: Nxd matrix of the measurement kets
%counts : multiplication factor that turns probabilities into estimated clicks
%output
%rhoPGDB: maximum likelihood density matrix 
%timePGDB: time taken by the algorithm
%costFunction: cost function value as a function of iteration number


function  [rhoPGDB, timePGDB, costFunction] = runPGDB(data,A,counts)
%
% Example of usage of 'pgqs' routine
% --------------------------------------
% Maximum Likelihood tomography of W state

% disp('PGDB')

t1 = tic;
sA = size(A);
d = sA(2);

n = d; % order of the density matrix
normalisedData = data/counts;


% other options for grrhor
tol = 0.1e-4;     % tolerance used in convergence criteria
rho0 = (1/n)*eye(n); % initial guess for the density matrix
maxit = 2000; % maximum number of iterations

% this is the function to be minimized
% (For this example, the function is the least squares function) 
fun = @(rho) leastsquares(rho,A,normalisedData ,counts); % function handle

% calling grrhor
[rho,flag,iter,fe,costFunction] = pgqs(fun,rho0,tol,maxit,1);

if (flag>=0)
%     fprintf('Succesfully solved with PGDB. Flag: %d, Iter: %d, nfeval: %d \n', flag, iter, fe);
%     % printing the recovered density matrix
%     rho;
else
%     fprintf('PGDB did not converge to the maximum likelihood state. \n');
end


timePGDB = toc(t1);
costFunction  = abs(costFunction);
rhoPGDB = rho;




