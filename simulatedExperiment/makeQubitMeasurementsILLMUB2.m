%produces the 'modified' tensor product SIC-POVM projectors and stacks them in the
%measurement matrix. 
%inputs
%n: number of qubits
%theta: angle between the project |0> and any of the three others (see
%Figure 1 in the paper)

% function [A] = makeQubitMeasurementsILLMUB(n,theta)



clear;clc

n = 1;
theta = 90*pi/180;
qq=0;
THETA = linspace(0,90,20)*pi/180;
THETA = 90*pi/180;


 for theta = THETA
qq=qq+1;
    
d=2^n;
[M] = permn([1 2 3 4 5 6], n);
ss = size(M);

lama = exp(1i*2*pi/3);
theta2 = theta+pi;

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


% % % % [s e u] = svd(A)
% % 
hon  = 1;
col = [0.8 0.1 0.9];
fig = 24;
yo = 3;

vec = A(1,:);
blochPlot3(vec,col,0,fig,yo);
for i=2:sA(1)
vec = A(i,:);
blochPlot3(vec,col,hon,fig,yo);
end

overlapMat = abs(A*A').^2;
ovs = reshape(overlapMat,[numel(overlapMat) 1]);
   
entovs(qq) = -sum(ovs.*log(ovs+eps));


end




% 
% figure(4)
% hold off
% plot(THETA*180/pi ,stdovs.^2)
% 
% axis([0 90 0 max(stdovs.^2)*1.2 ])
% 
% xlabel('Angle between bases')
% ylabel('Overlap variance')


% figure(5)
% hold off
% plot(THETA*180/pi ,entovs)
% 
% axis([0 90 0 max(entovs.^2)*1.2 ])
% 
% xlabel('Angle between bases')
% ylabel('Overlap Entropy')
% % 
