function q = euler2quat(eul)
% returns 4,1 vector with q = [q0;q1;q2;q3];

type = 'ZYX';

switch (type)
    case 'ZYX'
        alpha = eul(1);
        beta = eul(2);
        gamma = eul(3);
        
        T = [cos(beta)*cos(gamma),-cos(beta)*sin(gamma), sin(beta);
         cos(alpha)*sin(gamma) - cos(gamma)*sin(alpha)*sin(beta), cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma), cos(beta)*sin(alpha);
         sin(alpha)*sin(gamma) - cos(alpha)*cos(gamma)*sin(beta), -cos(gamma)*sin(alpha) + cos(alpha)*sin(beta)*sin(gamma),  cos(alpha)*cos(beta)];
end
        % 1)
        S = trace(T);
        
        q02 = 1/4*(1+S);
        q12 = 1/4*(1+2*T(1,1)-S);
        q22 = 1/4*(1+2*T(2,2)-S);
        q23 = 1/4*(1+2*T(3,3)-S);
        
        % get largest 
        q2 = [q02,q12,q22,q23];
        [~,id] = max(abs(q2));
        
        % get rest 
        switch id
            case 1 %q0
                q0 = sqrt(q02);
                q1 = (T(3,2)-T(2,3))/(4*q0);
                q2 = (T(1,3)-T(3,1))/(4*q0);
                q3 = (T(2,1)-T(1,2))/(4*q0);
            case 2 %q1
                q1 = sqrt(q12);
                q0 = (T(3,2)-T(2,3))/(4*q1);
                q3 = (T(1,3)+T(3,1))/(4*q1);
                q2 = (T(2,1)+T(1,2))/(4*q1);
            case 3
                q2 = sqrt(q22);
                q3 = (T(3,2)+T(2,3))/(4*q2);
                q0 = (T(1,3)-T(3,1))/(4*q2);
                q1 = (T(2,1)+T(1,2))/(4*q2);
            case 4
                q3 = sqrt(q23);
                q2 = (T(3,2)+T(2,3))/(4*q3);
                q1 = (T(1,3)+T(3,1))/(4*q3);
                q0 = (T(2,1)-T(1,2))/(4*q3);
        end
        
        % q0 < 0 -> * -1 !!
        q = [q0;q1;q2;q3];
        
        if (q0 < 0) 
            q = -1*q;
        end
        
                
        
        
end

