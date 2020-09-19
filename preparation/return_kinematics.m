
function [v,a,JT,a_bar,Dy,DDy,De,DDe] = return_kinematics(y,e,r)
    syms t real;

    Dy = sym('Dy',size(y),'real');
    DDy = sym('DDy',size(y),'real');

    De = sym('De',size(e),'real');
    DDe = sym('DDe',size(e),'real');

    % create taylor expansion
    y_ = y + t*Dy + t^2/2*DDy;
    e_ = e + t*De + t^2/2*DDe;
    
    
    %% 
    r = subs(r,[y;e],[y_;e_]);

    % calculate time dreivatives
    v = diff(r,t);
    a = diff(v,t);

    % replace taylor expansion (t = 0)
    r = subs(r,t,0);
    v = subs(v,t,0);
    a = subs(a,t,0);
    
    % global jacobian
    JT = jacobian(r,y);

    a_bar = a - JT*DDy;

end