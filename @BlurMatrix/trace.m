function c = trace(a)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

y = a.eigblurmatrix;
c = sum(y(:));
