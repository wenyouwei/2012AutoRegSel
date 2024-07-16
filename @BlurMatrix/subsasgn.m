function A = subsasgn(A,S,B)

switch S.subs
    case 'eigblurmatrix'
        A.eigblurmatrix = B;
end
        
