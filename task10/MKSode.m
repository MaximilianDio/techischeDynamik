% Maximilian Dio 23.04.2020

function dxdt = MKSode(x,ht)
% TODO: currently no time dependence 
% transforms system to statespace (system of first order ode)

    % get dimension of system
    n = size(x);
    n = n(1);
    
    % ode (first order)
    dxdt = zeros(n,1);
    dxdt(1:n/2) = x(n/2+1:end);
    dxdt(n/2+1:end) = ht(x);
end