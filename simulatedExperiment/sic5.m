function sic = sic5()

d=5;
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

phi = [+.39104489402214774638257588694092612854e+0+.00000000000000000000000000000000000000e+0i
-.28486558319586666154004262263615905414e+0-.64712933282796239400249892482493465752e+0i
-.23188384736899577826443055554623498413e+0-.19820390755555243362222174127453043647e+0i
+.13193857997561206409072011157833671545e+0-.10939599651964327439560846264470517498e+0i
+.43096743921096438598841830212227390108e+0+.19747405581954685663309013848038222203e+0i];


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