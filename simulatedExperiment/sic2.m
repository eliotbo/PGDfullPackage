function sic = sic2()

d=2;
w= exp(1i*2*pi/d);

dd = d-1;
for j=0:dd
    for k = 0:dd
        D{j+1,k+1}=zeros(d) ;
        for m=0:dd
        v1 = zeros(d,1);
        v2 = zeros(1,d);
        pos1  = mod(k+m,d);
        pos2 = m;
        v1(pos1+1) = 1;
        v2(pos2+1) = 1;
        D{j+1,k+1} = D{j+1,k+1} + w^(j*k/2+j*m)*v1*v2;
        end
    end
end

phi = zeros(d,1);
phi(1) = 1;
phi = rand(d,1);
phhi = phi/sqrt(phi'*phi);

phi = [+.88807383397711526216076459641812180401e+0+.00000000000000000000000000000000000000e+0i
+.32505758367186814316112416777511970282e+0-.32505758367186814316112416777511970282e+0i];


kk=0;
sic = zeros(d^2,d);
for j=1:d
    for k = 1:d
        kk=kk+1;
        sic(kk,:) = D{j,k}*phi;
    end
end

% A = zeros(d^4,d^2);
% kk=0;
% for j=1:d^2
%     for k = 1:d^2
%         kk=kk+1;
%         A(kk,:) = kron(sic(j,:),sic(k,:));
%     end
% end


% overlap= sic*sic';

% figure(51)
% imagesc(abs(overlap).^2)