function [alpha,beta,gamma] = rot2Kardan(S)
    %% transforms rotation matrix S to corresponding 
    % Tait-Bryan angles - error is thrown if singularity configuration is
    % reached (No unique solition possible!)
       
    beta = asin(S(1,3));
    
    % if beta = pi(1+2k) singularity occures!
    if (cos(beta) < 1e-2)
        error("singularity! angles can't be determined.");
    end
    
    alpha = acos(S(3,3)/cos(beta));
    gamma = acos(S(1,1)/cos(beta));
    
end