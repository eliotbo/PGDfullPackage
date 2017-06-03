function result = cost(A1,rhor,pData2)

   temp1 = A1.*(conj(A1)*rhor);
    Ax = sum(temp1,2);
    clear('temp1');
Ax = real(Ax);

% Ax = Ax/sum(Ax);
% pData2 = pData2 / sum(pData2);

   Axb = Ax-pData2;
    result = sum(Axb.^2);
    

%    result = sum(-pData2.*log(Ax+eps));