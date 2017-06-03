function [A1] = makeRandom2dRank1Projections(d,N)

twoD = rand(N,2).*exp(1i*2*pi*rand(N,2));
normD = 1./sqrt(sum(abs(twoD).^2,2));

twoDnormed = bsxfun(@times,twoD,normD);
a1 = twoDnormed;

[~, H] = sort(rand(N,d),2);
HH = H(:,1:2);

A1 = zeros(N,d);

idx1 = sub2ind(size(A1), 1:N , transpose(HH(:,1)));
idx2 = sub2ind(size(A1), 1:N , transpose(HH(:,2)));

A1(idx1) = a1(:,1);
A1(idx2) = a1(:,2);

% for ii=1:N
%     ai1 = rand(1,2).*exp(1i*2*pi*rand(1,2));
%     ai1 = ai1/sqrt(ai1*ai1');
%     aq1 = zeros(1,(d));
%     aq1pos = ceil(rand(1,2)*(d));
%     
%     while aq1pos(1) == aq1pos(2)
%         aq1pos = ceil(rand(1,2)*(d));
%     end
%     aq1(1,aq1pos(1)) = ai1(1);
%     aq1(1,aq1pos(2)) = ai1(2);
%     
%     aq1 = aq1+1E-10;
%     aq1 = aq1/sqrt(aq1*aq1');
%     ai = aq1;
%     A1(ii,1:d) = ai;
% end