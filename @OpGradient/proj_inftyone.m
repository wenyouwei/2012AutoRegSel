function q = proj_inftyone(p)
%function u = projection(p)
% Copyright Sep. 11, 2012, Dr.WEN You-Wei
% email: wenyouwei@gmail.google.com


px2 = p.x.^2;    py2 = p.y.^2;
dpxpy = sqrt( sum(px2,3) + sum(py2,3));
dpxpy(dpxpy<1) = 1;

q     = p;
for i = 1 : size(p.x, 3)
    q.x(:,:,i) = p.x(:,:,i) ./ dpxpy;                 
    q.y(:,:,i) = p.y(:,:,i) ./ dpxpy;  
end