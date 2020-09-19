function q = euler2Quaternions(d,phi)

assert( isequal(size(d),[3,1]) || isequal(size(d),[1,3]),"input must be 3 dimensional vector!");

% normalize direction vector
d = d/norm(d);

q0 = cos(phi/2);
q = [q0; d*sin(phi/2)];

end

