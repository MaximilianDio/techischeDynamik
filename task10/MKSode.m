% Maximilian Dio 23.04.2020

function dxdt = MKSode(t,f)
% TODO: currently no time dependence 
% transforms system to statespace (system of first order ode)

    % get dimension of system
    n = size(t);
    n = n(1);
    
    % ode (first order)
    dxdt = zeros(n,1);
    dxdt(1:n/2) = t(n/2+1:end);
    dxdt(n/2+1:end) = f(t);
end