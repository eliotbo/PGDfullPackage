function [f,g,p] = negloglikemulti(rho,E,r)
%
% function [f,g,p] = negloglikemulti(rho,E,r)
%
% Negative log-likelihood function for ML Quantum State Tomography
% (Multinomial distribution)
%
% INPUT
% rho           density matrix
% E             cell array containing the POVM matrices E{i}
%
% OUTPUT
% f             function value at rho
% g             gradient at rho
% p             auxiliar vector; p(i) = trace(rho*E{i})
%

% initializing
m = length(r);
p = zeros(m,1);
n = length(rho);

% computing probabilities
for i=1:m
    p(i) = max( real(trace(rho*E{i})), 0.0 )  ;
end

%fprintf('qmin: %f \n', min(p) );

f = sum(r.*log(p));
f = -f;

% gradient evaluation
if nargout>1
    g=zeros(n,n);
    for i=1:m
        g=g + (r(i)/(p(i) ))*E{i};
    end
    g=-g;
end

