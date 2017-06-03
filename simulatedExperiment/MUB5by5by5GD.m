function [A ]= MUB5by5by5GD(mubNum)

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


A = zeros(125^2,125);

    
    for i = 1:mubNum
        for k=1:5
            for ii = 1:mubNum
                for kk=1:5
                    for iii = 1:mubNum
                        for kkk=1:5
                            a=a+1;
                            proj = kron(kron(M{i}(k,:),(M{ii}(kk,:))),M{iii}(kkk,:));
                            A(a,1:125) = proj;
                        end
                    end
                end
            end
        end
    end






