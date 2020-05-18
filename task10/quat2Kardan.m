function [alpha,beta,gamma] = quat2Kardan(q0,q)
    %% transforms quaternions (q0 scalar and q vector) to corresponding 
    % Tait-Bryan angles - error is thrown if singularity configuration is
    % reached (No unique solition possible!)
       
    beta = asin(2*(q0*q(2)+q(1)*q(3)));
    
    % if beta = pi(1+2k) singularity occures!
    if (cos(beta) < 1e-2)
        error("singularity! angles can't be determined.");
    end
    
    % cosine is not bijective!
    alpha = acos((1-2*(q(1)^2+q(2)^2))/cos(beta));
    gamma = acos((1-2*(q(2)^2+q(3)^2))/cos(beta));
   
    
    % get right sign:
    if (alpha ~= 0)
        alpha = alpha*sign((q0*q(1)-q(2)*q(3))/(sin(alpha)));
    end
    if (gamma ~= 0)
        gamma = gamma*sign((q0*q(3)-q(2)*q(3))/(sin(gamma)));
    end
    
end