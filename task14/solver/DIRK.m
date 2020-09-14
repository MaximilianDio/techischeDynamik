%% Diagonal implicit RKM

function [t,y] = DIRK(f,tspan,y0,varargin)
    
    % butcher table -> diagonally implicit
    % check input arguments
    i = 1;
    while i <= length(varargin)
        switch cell2mat(varargin{i})
            case 'Method'
                i = i+1;
                switch cell2mat(varargin{i})
                    case '2Stage' 
                        s = 2; % number of stages
                        a = [1/2+sqrt(3)/6 ; 1/2-sqrt(3)/18];
                        b = [1/2+sqrt(3)/6 0; 
                            -2/(3*sqrt(3)) 1/2+sqrt(3)/6];
                        c = [1/2 1/2];
                    case 'trapezoidal'
                        s = 2; % number of stages
                        a = [0 ; 1];
                        b = [0 0; 
                            1/2 1/2];
                        c = [1/2 1/2];   
                    case 'implicitEuler'
                        s = 1; % number of stages
                        a = [1 ; 0];
                        b = [1 0; 
                            0 0];
                        c = [1 0]; 
                    case "exam"
                        s = 2;
                        phi = 1/2;
                        a = [0; 1];
                        b = [0 0; (1-phi) phi];
                        c = [(1-phi) phi];
                end
            otherwise
        end
        i = i+1;
    end
    
   

    %% handle input
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
    sizeT = size(tspan);
    tlength = length(tspan);
    if sizeT(1) == tlength
        t = tspan;
    else
        t = tspan';
    end
    y = zeros(tlength,N);
    
    % enter initial values:
    sizeY0 = size(y0);
    if size(sizeY0(1)) == N
        y0 = y0';
    end
    y(1,:) = y0;
    
    %%
    
    h = t(2)-t(1); %step size
    options = optimoptions('fsolve','Display','none');
    for n = 1:tlength-1
       
        % calc stage ki and solve implicit function
        k1 = @(x) f(t(n)+a(1)*h,y(n,:) + b(1,1)*h*x)-x;
        k1 = fsolve(k1,y(n,:),options);
        
        if s >= 2
            k2 = @(x) f(t(n)+a(2)*h,y(n,:)+b(2,1)*h*k1+b(2,2)*h*x)-x;
            k2 = fsolve(k2,k1,options);
        else
            k2 = 0;
        end
        y(n+1,:) = y(n,:) + h*(c(1)*k1+c(2)*k2);
        
    end
end