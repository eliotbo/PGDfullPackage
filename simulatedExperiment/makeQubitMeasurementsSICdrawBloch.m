%produces the 'modified' tensor product SIC-POVM projectors and stacks them in the
%measurement matrix. 
%inputs
%n: number of qubits
%theta: angle between the project |0> and any of the three others (see
%Figure 1 in the paper)

function [A] = makeQubitMeasurementsSICdrawBloch(n,theta)

n=1;
theta = 15*pi/180;

d=2^n;
[M] = permn([1 2 3 4], n);
ss = size(M);

lama = exp(1i*2*pi/3);

W = [1 0;
 cos(theta/2) sin(theta/2);
 cos(theta/2) lama*sin(theta/2);
 cos(theta/2) lama'*sin(theta/2)];

for i=1:4
    v = W(i,:);
    v = v/sqrt(v*v');
    W(i,:) = v;
end

sizeW = size(W);

for i=1:sizeW(1)
    v= W(i,:);
   s{i} =  v;
end

sizeM = size(M);
A = zeros(sizeM(1),2^n);

for i=1:ss(1)
    pos = M(i,:);
    m=1;
    for ii=1:n
        m=kron(m,s{pos(ii)});
    end
    v = m;
    A(i,:) = v;
end

hon  = 1;
col = [0.8 0.1 0.9];
fig = 24;
yo = 3;

vec = A(1,:);
blochPlot3(vec,col,0,fig,yo);

for i=2:4
    v=A(i,:);
    blochPlot3(v,col,hon,fig,yo) 
end