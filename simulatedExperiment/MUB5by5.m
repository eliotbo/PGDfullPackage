function [A ]= MUB5by5(q,mubNum)

M{1} = eye(5);

M2e =    [0 0 0 0 0;
    0 2 4 -4 -2;
    0 4 -2 2 -4;
    0 -4 2 -2 4;
    0 -2 -4 4 2];
M{2} = exp(M2e*1i*pi/5)/sqrt(5);

M3e = [0 2 -2 -2 2;
    0 4 2 4 0;
    0 -4 -4 0 -2;
    0 -2 0 -4 -4;
    0 0 4 2 4];
M{3} = exp(M3e*1i*pi/5)/sqrt(5);

M4e =[0 4 -4 -4 4;
    0 -4 0 2 2;
    0 -2 4 -2 0;
    0 0 -2 4 -2;
    0 2 2 0 -4];
M{4} = exp(M4e*1i*pi/5)/sqrt(5);

M5e = [0 -4 4 4 -4;
    0 -2 -2 0 4;
    0 0 2 -4 2;
    0 2 -4 2 0;
    0 4 0 -2 -2];
M{5} = exp(M5e*1i*pi/5)/sqrt(5);

M6e = [0 -2 2 2 -2;
    0 0 -4 -2 -4;
    0 2 0 4 4;
    0 4 4 0 2;
    0 -4 -2 -4 0];
M{6} = exp(M6e*1i*pi/5)/sqrt(5);

A5 = [M{1} ;M{2}; M{3}; M{4}; M{5}; M{6}];

a=0;

for aa=1:q
    if aa==1
        U=eye(5);
    else
        %compute a random unitary matrix
        X = (randn(5) + 1i*randn(5))/sqrt(2);
        [Q,R] = qr(X);
        R = diag(diag(R)./abs(diag(R)));
        U = Q*R;
    end

    %rotate the MUB is the random direction of the unitary matrix
    m{1} = U'*M{1}*U;
    m{2} = U'*M{2}*U;
    m{3} = U'*M{3}*U;
    m{4} = U'*M{4}*U;
    m{5} = U'*M{5}*U;
    m{6} = U'*M{6}*U;
    
    for i = 1:mubNum
        for k=1:5
            for ii = 1:mubNum
                for kk=1:5
                    a=a+1;
                    proj = kron(M{i}(k,:),(M{ii}(kk,:)));
                    A(a,1:25) = proj;
                end
            end
        end
    end
end





