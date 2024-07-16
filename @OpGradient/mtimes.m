function q = mtimes(a,p)
%function y = mtimes(a, x)
% y = G * x	 
% Copyright Sep. 1, 2012, Dr.WEN You-Wei
% email: wenyouwei@gmail.google.com

if ~isa(a, ' OpGradient') 
    q = p;
    if numel(a) ~= 1
        error('Wrong, the paramter must be a scalar and a class');
    else
        q.x = a * p.x;
        q.y = a * p.y;
    end
else
    if numel(p) ~= 1
        error('Wrong, the paramter must be a scalar and a class');
    else
        q = a;
        q.x = a.x * p;
        q.y = a.y * p;
    end
end