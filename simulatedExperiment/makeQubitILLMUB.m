%produces the n-qubit MUB bases

function [A entropy] = makeQubitILLMUB(n,theta)


% n = 2;
% theta = 90*pi/180;
% qq=0;
% THETA = linspace(0,90,20)*pi/180;
% THETA = ;
% for theta = THETA
% qq=qq+1;
    

d=2^n;
[M] = permn([1 2 3 4 5 6], n);
ss = size(M);

lama = exp(1i*2*pi/3);
theta2 = theta+pi;

 W = ( [1	  0;
     0 1;
     cos(theta/2)  sin(theta/2);
     sin(theta/2) -cos(theta/2) ;
     1*cos(theta/2)  1i*sin(theta/2)
      1*sin(theta/2)  -1i*cos(theta/2)]);
 
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
A = (zeros(sizeM(1),2^n));

for i=1:ss(1)
    pos = M(i,:);
    m=1;
    for ii=1:n
        m=(kron(m,s{pos(ii)}));
    end
    v = m;
    A(i,:) = v;
end

sA = size(A);

% % % [s e u] = svd(A)
% % 
% hon  = 1;
% col = [1 0 0];
% fig = 24;
% yo = 4;
% 
% vec = A(1,:);
% blochPlot3(vec,col,0,fig,yo);
% for i=2:sA(1)
% vec = A(i,:);
% blochPlot3(vec,col,hon,fig,yo);
% end

% k=0;
% for i=1:sA(1)
%     for ii=1:sA(1)
%         
%         va = A(i,:);
%         vb = A(ii,:);
%         ov = abs(va*vb')^2;
%         overlap(i,ii) = ov;
% 
%             k=k+1;
%             ovs(k) = ov;
%     end
% end
 


if nargout==2
    entropy = entropyMeasurements(A);
    entovsMUB = qubitMUBentropy(n);
    entropy = entropy/entovsMUB;
else
    entropy = 0;
end

% figure(5)
% hold on
% plot(THETA*180/pi ,entovs,'.r')
% 
% axis([0 90 0 1 ])
% 
% xlabel('Angle between bases')
% ylabel('Overlap Entropy')

