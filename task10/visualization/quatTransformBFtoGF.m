function [Xg,Yg,Zg] = quatTransformBFtoGF(Xb,Yb,Zb,rOg,q)
    % tranforms bodyframe to groundframe with translation vector and
    % rotation d (direction) + phi (angle in deg). rOg position to origin in ground frame. 
    q = reshape(q,[],1);
    assert(abs(q'*q - 1) <= 1e-4,'norm of direction has to be one!');
    
    % get original size of Vector (might be matrix for complex shape e.g.
    % cylinder)
    [n,m] = size(Xb);
    
    % convert vector to pos in body frame
    rPb = [reshape(Xb,1,[]);reshape(Yb,1,[]);reshape(Zb,1,[])];
    
    % rotation matrix (Euler-Rodrigues formula)
    tilde = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];
    Sgb = eye(3)+2*q(1)*tilde(q(2:4))+2*tilde(q(2:4))^2;
    
    % transform rPb to ground frame
    rPg = rOg + Sgb*rPb;
    
    Xg = reshape(rPg(1,:),n,m);
    Yg = reshape(rPg(2,:),n,m);
    Zg = reshape(rPg(3,:),n,m);
   
end
