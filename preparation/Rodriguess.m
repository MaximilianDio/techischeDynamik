function S = Rodriguess(d,phi)
%RODRIGUESS returns the rotation matrix based on rotation axis d and angle
%phi. the axis has to be of norm 1!

% normalize
d = d/norm(d);

S = eye(3)+tilde(d)*tilde(d)*(1-cos(phi))+tilde(d)*sin(phi);

end

