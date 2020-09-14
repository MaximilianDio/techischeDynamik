%% Runge kutta method explicit

function [t,y] = explicitRKM(f,tspan,y0,varargin)
    
    tol = 1e-3; %standard relative tolerance if embedded method is chosen

    % check input arguments
    i = 1;
    while i <= length(varargin)
        switch cell2mat(varargin{i})
            case 'Method'
                i = i+1;
                switch cell2mat(varargin{i})
                    case 'expEuler' %explicit euler
                        p = 1;
                        embedded = false;
                        a = [0 ; 0];
                        b = [0 0; 0 0];
                        c1 = [1 0];
                    case 'Heun'
                        embedded = false;
                        p = 2;
                        a = [0 ; 1];
                        b = [0 0; 1 0];
                        c1 = [1/2 1/2];
                    case 'embedded'
                        embedded = true;
                        p = 2;
                        a = [0 ; 1];
                        b = [0 0; 1 0];
                        c1 = [1/2 1/2];
                        c2 = [1 0];
                    otherwise %improved euler
                        embedded = false;
                        a = [0 ; 1/2];
                        b = [0 0; 1/2 0];
                        c1 = [1/2 1/2];
                end
            case "Tol"
                i = i+1;
                tol = varargin{i};
                assert(isnumeric(tol),"tolerance must be numeric!");
            otherwise
        end
        i = i+1;
    end

    
    % function handle f of type f = @(t,y) dy(t,y)   
    % transform f into row vector output
    sizeF = size(f(0,y0));
    if sizeF(1) == length(f(0,y0))
        f = @(t,y) f(t,y)';
        N = sizeF(1);
    else
        N = sizeF(2);
    end
    
    
    % preallocate vectors
    tend = tspan(end);
    sizeT = size(tspan);
    tlength = length(tspan);
    if sizeT(1) == tlength
        t = tspan';
    else
        t = tspan;
    end
    y = zeros(tlength,N);
    
    % enter initial values:
    sizeY0 = size(y0);
    if size(sizeY0(1)) == N
        y0 = y0';
    end
    y_(1,:) = y0;
    
    %% run RKM 
    n = 1;
    t_(1) = t(1);
    h_ = t(2)-t(1);
    h = h_; % adapted step size - initially same;
    repeat_step = false;
    while t_(n) < tend
        if repeat_step == false
            H = h_; % change to adapted step size
        else 
            H = h; % change to new step size
        end
        t_(n+1) = t_(n) + H;
        
        % calc stage ki
        k1 = f(t_(n),y_(n,:));
        k2 = f(t_(n)+a(2)*H,y_(n,:)+b(2,1)*H*k1);
        
        % combine stages to new solution
        y_(n+1,:) = y_(n,:) + H*(c1(1)*k1 + c1(2)*k2);
        if embedded == true
           % second solution of lower order
           y_hat = y_(n,:) + H*(c2(1)*k1 + c2(2)*k2);
           % estimated error
           h = H*(tol/norm(y_(n+1,:)-y_hat))^(1/p);
           if H-h > 0.001
                % current step size not enough repeat step!
                repeat_step = true;
                continue;
           end
        end
        
        % step ok -> advance
        repeat_step = false;
        n = n+1;
    end
    
%     if embedded == true
%         disp("number of iterations: " + string(length(t_))); 
%     end
    % interpolate s.t. solution exists on input grid
    y = interp1(t_,y_,t);


end