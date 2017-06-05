%This example simulates measurements made on a randomly generated density 
%matrix and applies the PGDM method (the unkown-r version). One can change
%the parameters of the simulations (dimensionality, purity, number of 
%clicks, type of measurement), but the current version uses
%a qubit SIC-POVM on a four-qubit 0.5 purity system with 10000 counts per
%outcome on average.

function [fid,timeTaken] = example1(param,A)


% [rho, A, data] = generateDataset(param);
%%% if the entropy is needed uncomment the following line (comp. intensive)


% figure(222);semilogy(1,1);hold on

r = param.counts; %add as third input to fullTomography if known

%%%%%% if you have data of your own, the following lines are useful %%%%%%%
[rho, data] = generateDataset(param,A);
[rhoEstimates, timeTaken, costs] = fullTomography(data,A,r); %note that 'r' is optional here.

% close all
% figure(222); hold off
%calculate fidelities and plot log likelihoods as a function of iteration
fields = fieldnames(costs);
fields2 = fieldnames(timeTaken);

for i = 1:numel(fields2)
    fid.(fields2{i}) = fidelityRho(rho,rhoEstimates.(fields2{i}));
end

for i = 1:numel(fields)
     semilogy(costs.(fields{i}),'-'); hold on
%     legend222{i} = fields{i};
end
% legend(legend222)
% grid on

% disp(' ')
% disp('The time taken by each method in seconds is')
% disp(timeTaken)
% disp('The fidelity between the actual state ''rho'' and the maximum likelihood state ''rhoEstimate'' is:')
% disp(fid)

% if exist('entropyA')
%     disp(['The normalised measurement matrix overlap entropy is ' num2str(entropyA)])
% end
