function b_tilde = tilde(b)
%TILDE returns corssproduct matrix of 3D vector b
assert( isequal(size(b),[3,1]) || isequal(size(b),[1,3]),"input must be 3 dimensional vector!");
b_tilde = [0 -b(3) b(2); b(3) 0 -b(1); -b(2) b(1) 0];
end

