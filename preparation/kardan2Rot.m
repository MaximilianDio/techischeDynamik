function S = kardan2Rot(alpha,beta,gamma)

S_AB = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
S_BC = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
S_CD = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];

S = S_AB*S_BC*S_CD;

end

