function [x] = simplex_proj(a)
%
% function [x] = simplex_proj(a)
%
% Projects a vector 'a' onto the unit simplex
%
% Reference: C. Michelot. JOTA, 1986 
%

% size of a
n=length(a);

x=a;

I = 1:n;

for t=1:n
    d=length(I);
    
    % projection on V_I
    x(I) = x(I) + (1/d)*(1 - sum(x(I)));
    
    N=find(x<0);
    
    % stopping criterion
    if isempty(N)
        break;
    end
    
    % projection on X_I
    x(N) = 0.0;
    
    % update I
    I=find(x>0);
    
end

