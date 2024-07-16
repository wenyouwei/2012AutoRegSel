function y = inv(a)
% Copyright Jan. 25, 2008, Dr.WEN You-Wei
% email: wenyouwei@graduate.hku.hk

y = a;

switch a.boundarycond
    case {'cir','refl'}
        y.eigblurmatrix = 1./(a.eigblurmatrix);
    otherwise
        disp('Unknow method');
end
