
function [Dc,DDc,C,ct,ctt,Dx,DDx,De_t,DDe_t] = boundary_condtions(c,x,e)
    syms t real;

    Dx = sym('x',size(x),'real');
    DDx = sym('Dx',size(x),'real');

    De_t = sym('De',size(e),'real');
    DDe_t = sym('DDe',size(e),'real');

    x_ = x + t*Dx + t^2/2*DDx;
    e_ = e + t*De_t + t^2/2*DDe_t;

    c = subs(c,[x;e],[x_;e_]);
    % calculate time dreivatives
    Dc = diff(c,t);
    DDc = diff(Dc,t);

    % replace taylor expansion (t = 0)
    c = subs(c,t,0);
    Dc = subs(Dc,t,0);
    DDc = subs(DDc,t,0);

    C = simplify(jacobian(c,x));

    ct = simplify(Dc - C*Dx);
    ctt = simplify(DDc - C*DDx);
end