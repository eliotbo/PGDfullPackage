function A1 = makeRandomMeasurementsOneSystem(D,N)


    d=D;
    A1 = [];
    
for i= 1:N
    randM = rand(d).*exp(1i*2*pi*rand(d));

    [Q, R] = qr(randM);
    U = Q*diag(sign(diag(R)));
    
    
    A1 = [A1; U'];
    
%     for j=1:d
%         v= U(:,j);
%         for k = 1:d
%            w = U(:,k);
%            V = kron(v,w);
%            A1 = [A1; V'];
%         end
%          
%     end
    
   
end
    