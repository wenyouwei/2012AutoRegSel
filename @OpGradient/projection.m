function q = projection(p,flag)
%function u = projection(p)
% Copyright Sep. 11, 2012, Dr.WEN You-Wei
% email: wenyouwei@gmail.google.com

if nargin == 1
    flag = 'infty';
end

px2 = p.x.^2;    py2 = p.y.^2;
dpxpy = sqrt( sum(px2,3) + sum(py2,3));

switch lower(flag)
    case 'infty'
        dpxpy(dpxpy<1) = 1;       
        q     = p;
        for i = 1 : size(p.x, 3)
            q.x(:,:,i) = p.x(:,:,i) ./ dpxpy;                 
            q.y(:,:,i) = p.y(:,:,i) ./ dpxpy;  
        end
    case 'l1'
        if sum(dpxpy)<=1, 
            q = p; 
        else
            [tmp,mu] = ProjLoneUpBound(dpxpy,1);
            q   = p;
            tmp = (1 - mu./ dpxpy);
            tmp(dpxpy<mu) = 0;
            for i = 1 : size(p.x, 3)
                q.x(:,:,i) = p.x(:,:,i) .* tmp;                 
                q.y(:,:,i) = p.y(:,:,i) .* tmp; 
            end
        end
        
end


function [f,mu] = ProjLoneUpBound(g,c)
%% min||f-g||_2^2 s.t. ||f||_1<=c

a = sum(abs(g(:)));mu = 0;
if a-c<1e-10, f = g*c/a; return; end
u = abs(g(:));
u = sort(u);
E = u;
L = length(u);
usum_k2end = sum(u(:));
for k = 1:L
    if k == 1
        E(1) = usum_k2end -  u(k)*(L-k+1);
    else
        usum_k2end = usum_k2end - u(k-1);
        E(k) = usum_k2end - u(k)*(L-k+1);
    end
    if E(k)<c,
        if k == 1
            if sum(abs(g(:)))<=c, 
                f = g; 
            else
                f = g/sum(abs(g(:)));
            end
            
            return;
        end
        a1 = u(k-1);
        b1 = E(k-1);
        a2 = u(k);
        b2 = E(k);
        mu = (a2-a1)*c + b2 * a1 - b1 * a2; 
        mu = mu/(b2-b1);
        f = abs(g);
        f = sign(g).*(f >= mu).*(f - mu); 
        return;
    end
end
