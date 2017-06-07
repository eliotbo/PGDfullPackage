%This example simulates measurements made on a randomly generated density 
%matrix and applies the PGDM method (the unkown-r version). One can change
%the parameters of the simulations (dimensionality, purity, number of 
%clicks, type of measurement), but the current version uses
%a qubit SIC-POVM on a four-qubit 0.5 purity system with 10000 counts per
%outcome on average.

function example1()

clc;

%add subfolders as paths
currentFolderContents = dir(pwd);     
currentFolderContents (~[currentFolderContents.isdir]) = [];
for i = 3:length(currentFolderContents)           
   addpath(['./' currentFolderContents(i).name]) ;
end

%%%%%%%%%%%%%%%%%%%%%  simulation parameters  %%%%%%%%%%%%%%%%%%%%%%%
param.d =8;                   %dimensionality of the density matrix
param.purity = 0.5;            %purity of the density matrix
param.counts = param.d*1E5;    %this is the multication factor r
r = param.counts;              
param.theta =60*pi/180;
%Notes: 
%1) Unit purity states render PGD algorithms unstable, but a trick can be
%performed inside to bypass the unstabilities (see runPGDM.m)
%2) The total number of detector clicks is N*param.counts/param.d, where N
%is the total number of outcomes (or projections in A)
%3) Theta is the angle between the project |0> and any of the three others
%in the Bloch sphere
%%%%%%%%%%%%%%%%%%%%%  simulation parameters  %%%%%%%%%%%%%%%%%%%%%%%%

%The density matrix rho has a fixed purity but is generated randomly 
%the data is computer generated with Poissonian noise
[rho, A, data] = generateDatasetAndMeas(param);

figure(222);semilogy(1,1);hold off

%call main function to process the simulated data with PGDM, FISTA, PGDB
%and DIA. CVX runs only if it is installed on the used machine
[rhoEstimates, timeTaken, costs] = fullTomography(data,A,r); 

figure(222); 
%calculate fidelities and plot log likelihoods as a function of iteration
fields = fieldnames(rhoEstimates);
sizeFields = numel(fields);
for i = 1:sizeFields
  fid.(fields{i}) = fidelityRho(rho,rhoEstimates.(fields{i}));
  if ~strcmp(fields{i},'CVX')
      legend222{i} = fields{i};
      if strcmp(fields{i},'PGDM')
         semilogy(costs.(fields{i}),'--'); hold on
      elseif strcmp(fields{i},'FISTA')
         semilogy(costs.(fields{i}),'-.'); hold on
      else
         semilogy(costs.(fields{i})); hold on
      end
  end
end
legend(legend222)
grid on

disp(' ')
disp('The time taken by each method in seconds is')
disp(timeTaken)
disp('The fidelity between the actual state ''rho'' and the maximum likelihood state ''rhoEstimate'' is:')
disp(fid)

