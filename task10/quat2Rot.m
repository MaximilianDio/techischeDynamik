function S = quat2Rot(q0,q)
%% transforms quaterions scalar q0 and vector q to 3D rotational matrix
        
    % crossproduct-matrix
    tilde = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
    
    % euler-rodrigues formula
    S = eye(3) + 2*q0*tilde(q) + 2*tilde(q)*tilde(q);
    
    % TODO: formulate expilicitly 
    
end