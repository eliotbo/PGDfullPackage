% function [f,g,p,costfunction] = leastsquares(rho,E,r)
function [f,g,p,costfunction] = leastsquares(rho,A1,r,counts)
%
% function [f,g,p] = leastsquares(rho,E,r)
%
% Least squares function for ML Quantum State Tomography
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
% m = length(r);
% p = zeros(m,1);
% n = length(rho);

% % computing probabilities
% for i=1:m
%     p(i) = max( real(trace(rho*E{i})), 0.0);
% end

    temp1 = conj(A1).*((A1)*rho);
    Ax = sum(temp1,2)+1E-10;

p = Ax;

f = 0.5*sum((Ax - r).^2);

% costfunction = sum((Ax/sum(Ax) - r/sum(r)).^2);
costfunction = (counts)*sum( ((Ax-r)./sqrt(Ax+eps))  .^2)/length(Ax);

Axb = Ax-r;
g = (bsxfun(@times,A1,Axb)'*A1);

% gradient evaluation
% if nargout>1
%     g=zeros(n,n);
%     for i=1:m
%         g = g + (p(i) - r(i))*E{i};
%     end
% end

