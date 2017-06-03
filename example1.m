%This example simulates measurements made on a randomly generated density 
%matrix and applies the PGDM method (the unkown-r version). One can change
%the parameters of the simulations (dimensionality, purity, number of 
%clicks, type of measurement), but the current version uses
%a qubit SIC-POVM on a four-qubit 0.5 purity system with 10000 counts per
%outcome on average.

%function example1()

clc;clear;%%

%add subfolders as paths
currentFolderContents = dir(pwd);     
currentFolderContents (~[currentFolderContents.isdir]) = [];
for i = 3:length(currentFolderContents)           
   addpath(['./' currentFolderContents(i).name]) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Change parameters at will! 
%Notes: 
%1) Unit purity states render PGD algorithms unstable, but a trick can be
%performed inside to bypass the unstabilities (see runPGDM.m)
%2) The total number of detector clicks is N*param.counts/param.d, where N
%is the total number of outcomes (or projections in A)
%3) Try special systems and measurement matrices by setting d=25, 50 or 125
%3-b) if d=125, CVX 
%4-a) The modified SIC-POVM for qubits is generated if d=2^n.
%4-b) Theta is the angle between the project |0> and any of the three others 
%(see Figure 1 in the paper). We used theta = 60 degrees in the paper
param.d = 8;                   %dimensionality of the density matrix
param.purity = 0.5;             %purity of the density matrix
param.counts = param.d*1E4;     %multication factor r
param.theta = 90*pi/180;       %only used if d=2^n (many qubits). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [rho, A, data] = generateDataset(param);
%%% if the entropy is needed uncomment the following line (comp. intensive)
[rho, A, data] = generateDataset(param);

figure(222);semilogy(1,1);hold on

r = param.counts; %add as third input to fullTomography if known

%%%%%% if you have data of your own, the following lines are useful %%%%%%%
[rhoEstimates, timeTaken, costs] = fullTomography(data,A,r); %note that 'r' is optional here.

% close all
figure(222); 
%calculate fidelities and plot log likelihoods as a function of iteration
fields = fieldnames(costs);
for i = 1:numel(fields)
  fid.(fields{i}) = fidelityRho(rho,rhoEstimates.(fields{i}));
  %if ~strcmp(fields{i},'CVX')
      legend222{i} = fields{i};
%       if strcmp(fields{i},'PGDM')
%          semilogy(costs.(fields{i}),'--'); hold on
%       elseif strcmp(fields{i},'FISTA')
%          semilogy(costs.(fields{i}),'-.'); hold on
%       else
%          semilogy(costs.(fields{i})); hold on
%       end
  %end
end
legend(legend222)
grid on

disp(' ')
disp('The time taken by each method in seconds is')
disp(timeTaken)
disp('The fidelity between the actual state ''rho'' and the maximum likelihood state ''rhoEstimate'' is:')
disp(fid)

if exist('entropyA')
    disp(['The normalised measurement matrix overlap entropy is ' num2str(entropyA)])
end
