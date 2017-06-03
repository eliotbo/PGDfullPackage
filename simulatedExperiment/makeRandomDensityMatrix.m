function [rho] = makeRandomDensityMatrix(lambda)

  %generate random density matrix
    d = length(lambda);
    rho = zeros(d);
    randM = rand(d).*exp(1i*2*pi*rand(d));

    [Q, R] = qr(randM);
    U = Q*diag(sign(diag(R)));
    
    for ii=1:d
        psi = U(:,ii);
        rho = rho + psi*psi'*lambda(ii);
    end
    
% rho = zeros(d);
% rho(1) = 0.99;
% rho(2,2) =1-rho(1);