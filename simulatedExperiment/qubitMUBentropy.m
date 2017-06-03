%produces the 'modified' tensor product SIC-POVM projectors and stacks them in the
%measurement matrix. 
%inputs
%n: number of qubits
%theta: angle between the project |0> and any of the three others (see
%Figure 1 in the paper)

 function [entropyQubit] = qubitMUBentropy(n)

theta = 90*pi/180;

d=2^n;
[M] = permn([1 2 3 4 5 6], n);
ss = size(M);

 W = [1	  0;
     0 1;
     cos(theta/2)  sin(theta/2);
     sin(theta/2) -cos(theta/2) ;
     1*cos(theta/2)  1i*sin(theta/2)
      1*sin(theta/2)  -1i*cos(theta/2)];

sW = size(W);
for i=1:sW(1)
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

sA = size(A);


% k=0;
% for i=1:sA(1)
%     for ii=1:sA(1)
%         
%         va = A(i,:);
%         vb = A(ii,:);
%         ov = abs(va*vb')^2;
%         overlap(i,ii) = ov;
% 
%         k=k+1;
%         ovs(k) = ov;
%     end
% end
% 
% 
%     
% entovs = -sum(ovs.*log(ovs+eps));

entropyQubit = entropyMeasurements(A);

% overlapMatrix = abs(A*A').^2;
% ovs = reshape(overlapMatrix,[numel(overlapMatrix) 1]);
% 
% entovs = -sum(ovs.*log(ovs+eps));




