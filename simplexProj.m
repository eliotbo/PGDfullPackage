function [rhoOut]= simplexProj(rho)

d=length(rho);

[V,Dunsort] = eig(rho);

Dun = real(diag(Dunsort));
[D,I] = sort(Dun,1,'descend');
V = V(:, I);

D2 = D;

% for k= 1:100
%     D = D + (1-sum(D))/d;
%     D(D<0) = 0;
% end

for k=1:d
    D2(1:(d+1-k)) = D2(1:(d+1-k)) + (1-sum(D2(1:(d+1-k))))/(d+1-k);
    if abs(sum(D2(D2<0)))== 0
        break; 
    else D2(d+1-k) = 0; 
    end
end

Dl2 = diag(D2);
rhoOut = V*Dl2*V';



