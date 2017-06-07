%This example imports the data and feeds it to the tomography algorithm.
%Here, the measurement matrix A forms a POVM, and we can calculate the 
%multiplicative factor 'r' that turns probabilities into estimated 
%detector clicks. If 'r' is unknown, uncomment runPGDMfreeTrace(data,A).

% function rhoEstimates = example2()

clc;

%add subfolders as paths
currentFolderContents = dir(pwd);     
currentFolderContents (~[currentFolderContents.isdir]) = [];
for i = 3:length(currentFolderContents)           
   addpath(['./' currentFolderContents(i).name]) ;
end

%load dummy density matrix
load('rho')
%Measurement matrix: rows correspond to rank-1 normalised projectors
load('A');sA = size(A);
%load Nx1 raw-data vector of detector clicks
load('data')
%if A forms a POVM, r = d*sum(data)/N
r = sum(data)*sA(2)/sA(1);

%run a tomography algorithm: comment/uncomment to change algorithm
[rhoEstimates.PGDM, timeTaken.PGDM, costs.PGDM]   = runPGDM(data,A,r);
% [rhoEstimates.FISTA, timeTaken.FISTA, costs.FISTA] = runFISTA(data,A,r);
% [rhoEstimates.SPGD, timeTaken.SPGD, costs.SPGD]   = runSPGD(data,A,r);
% [rhoEstimates.DIA timeTaken.DIA, costs.DIA] = runDIA(data,A,r); 
% [rhoEstimates.CVX, timeTaken.CVX] = runCVX(data,A,r);
% 
%%% Only PGDMfreeTrace and CVX are written to accept data without r %%%%
% [rhoEstimates.PGDMfT, timeTaken.PGDMfT, costs.PGDMfT] = runPGDMfreeTrace(data,A);
% [rhoEstimates.CVX, timeTaken.CVX] = runCVX(data,A);

%detect which method ran
fields = fieldnames(rhoEstimates);

%plot log-likelihood log L_{GP} for each method
figure(222);hold off
semilogy(costs.(fields{1}))
grid on
legend(fields{1})

%caculate fidelity between dummy density matrix and recovered one
fid = abs(fidelityRho(rho,rhoEstimates.(fields{1})));

disp('The recovered density matrix is:')
disp(rhoEstimates.PGDM)
disp('The dummy density matrix is:')
disp(rho)
disp(['The fidelity between the dummy density matrix and the recovered '...
    'one is ' num2str(fid)])