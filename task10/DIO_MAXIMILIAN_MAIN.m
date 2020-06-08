%% Task 10 - bicycle wheel gyroscope
% Maximilian Dio - 21595892 -17.05.2020

close all; clear; clc;

%flags
ANIMATE = false;
SIMULATION = 2; % 1 or 2
SAVE_PLOTS = false;

%% parameters / constants
g = 9.81;   % [m/s^2] gravitational constant

R = 0.4;    % [m] radius of wheel
th = 0.02;  % [m] thickness of wheel
h = 0;      % [m] witdh of wheel
l = 0.25;   % [m] eccentricity
m = 10;     % [kg] mass 

%% initial conditions in Tait-Bryan angles
alpha0 = deg2rad(11);
beta0 = 0;
gamma0 = 0;

omega0_XB = 0;
omega0_YB = 0;
if (SIMULATION == 1)
    omgea0_ZB = 0; 
    t_end = 5; % [sec]
else
    omgea0_ZB = 10*pi;
    t_end = 30; % [sec]
end

%% helper functions 
% crossproduct-matrix
tilde = @(x) [0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

%% symbolic variables for analyitcal mechanics
% quaternions
syms q0 q1 q2 q3 Dq0 Dq1 Dq2 Dq3 D2q0 D2q1 D2q2 D2q3 t_ real
q = [q1;q2;q3];
Dq = [Dq1;Dq2;Dq3];
D2q = [D2q1;D2q2;D2q3];

y = [q0;q];
Dy = [Dq0;Dq];
D2y = [D2q0; D2q];
% taylor expansion for symbolic differentiation
y_ = y + t_*Dy + 1/2*t_^2*D2y;
Dy_ = Dy + t_*D2y;

%% rotation matices 
% euler-rodrigues formula
S_BG = eye(3) - 2*q0*tilde(q) + 2*tilde(q)*tilde(q);
S_GB = eye(3) + 2*q0*tilde(q) + 2*tilde(q)*tilde(q); %S_GB = S_BG'

%% moment of inertia
thet = m*(R^2+(R-th)^2);
Thet_B = diag([1/4*thet,1/4*thet,1/2*thet])+m*l^2*diag([1,1,0]);
% augment tensor
Thet_hat_B = [1 zeros(1,3); zeros(3,1) Thet_B];

%% location vector
r_B = [0;0;l];

%% external forces/ moments 
tau_B = cross(r_B,S_BG*m*[0;g;0]);
% augment vector
tau_hat_B = [0;tau_B];

%% kinematic ode (sytem matrix)
L_BG = [-q, q0*eye(3)-tilde(q)];
H_BG = 2*[q0, q';L_BG];

L_GB = [-q, q0*eye(3)+tilde(q)];
H_GB = 2*[q0, q';L_GB];

%% derive system of equation
omega_hat_B = H_BG*Dy;
omega_hat_B_ = subs(omega_hat_B,[y,Dy],[y_,Dy_]);

Domega_hat_B = simplify(subs(diff(omega_hat_B_,'t_'),t_,0));
alpha_hat_B = Domega_hat_B - H_BG*D2y;

D2y = (Thet_hat_B*H_BG)\(tau_hat_B - Thet_hat_B*alpha_hat_B - [1 zeros(1,3);zeros(3,1), tilde(L_BG*Dy)]*Thet_hat_B*omega_hat_B);

%% simulation
f = matlabFunction(D2y,'vars',{[q0;q1;q2;q3;Dq0;Dq1;Dq2;Dq3]});

% transform initial consitions
x0 = euler2quat([alpha0,beta0,gamma0]);
H0_BG = double(subs(H_BG,y,x0));
Dx0 = H0_BG\[0; omega0_XB; omega0_YB; omgea0_ZB]; 

% solving ODE numericaly
N = 1500; % number of samples
t_span = linspace(0,t_end,N);

opts = odeset('RelTol',1e-9,'AbsTol',1e-9);
[t_span,y_nonlin] = ode45(@(t,x) MKSode(x,f),t_span,[x0;Dx0],opts);

%% post processing 

% quaternions will drift
y_sol = y_nonlin(:,1:4);
unit = zeros(N,1);
for ii = 1:N
    unit(ii) = y_sol(ii,:)*y_sol(ii,:)'-1;
end

% location of COM
r_G = double(subs(S_GB*r_B,{q0,q1,q2,q3},{y_nonlin(:,1)',y_nonlin(:,2)',y_nonlin(:,3)',y_nonlin(:,4)'}));

% get Tait-Bryan angles
angles = zeros(N,3);
for ii = 1:N
    angles(ii,:) = quat2Kardan(y_nonlin(ii,1),y_nonlin(ii,2:4));
end

% get angular velocities in ground frame
omega_hat_G = double(subs(H_GB*Dy,{q0,q1,q2,q3,Dq0,Dq1,Dq2,Dq3},{y_nonlin(:,1)',y_nonlin(:,2)',y_nonlin(:,3)',y_nonlin(:,4)',y_nonlin(:,5)',y_nonlin(:,6)',y_nonlin(:,7)',y_nonlin(:,8)'}));
omega_G = omega_hat_G(2:4,:);

omega_hat_B = double(subs(omega_hat_B,{q0,q1,q2,q3,Dq0,Dq1,Dq2,Dq3},{y_nonlin(:,1)',y_nonlin(:,2)',y_nonlin(:,3)',y_nonlin(:,4)',y_nonlin(:,5)',y_nonlin(:,6)',y_nonlin(:,7)',y_nonlin(:,8)'}));
omega_B = omega_hat_B(2:4,:);

%% visualization

% plots
if SIMULATION == 1
    N_fig = 2;
    fig = gobjects(N_fig,1); 
    
    % plot yz location of com
    fig(1) = figure('Name','YZ_Location_of_COM');
    fig(1).Position = [680 676 560 302];
    axes; hold on; grid on; 
    xlabel('z [m]','FontSize', 12); 
    ylabel('y [m]','FontSize', 12); 
    title('yz-location of COM','FontSize', 14);
    set(gca,'YAxisLocation','origin','ydir','reverse');
    xlim([-l*1.1 l*1.1]); ylim([0 l*1.1]); pbaspect([2 1 1]);
    
    plot(r_G(3,:),r_G(2,:),'LineWidth',1.5);
    
     % COM in ground frame
    fig(2) = figure('Name','Location_COM_Groud_Frame');
    fig(2).Position = [ 680 476 1159 502];
    
    subplot 311; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('x_G','FontSize', 12); 
    title('x-position COM in ground frame','FontSize', 14);
    plot(t_span,r_G(1,:)','LineWidth',1.5,'Color',[0, 0.4470, 0.7410]); 
    
    subplot 312; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('y_G','FontSize', 12); 
    title('y-position COM in ground frame','FontSize', 14);
    plot(t_span,r_G(2,:)','LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980]);
    set(gca,'ydir','reverse');
    
    subplot 313; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('z_G','FontSize', 12); 
    title('z-position COM in ground frame','FontSize', 14);
    plot(t_span,r_G(3,:)','LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250]); 
        
elseif SIMULATION == 2
    N_fig = 4;
    fig = gobjects(N_fig,1); 
    
    % divergence of quaternions
    fig(1) = figure('Name','Quaternions_over_time'); 
    fig(1).Position = [ 680 476 1159 502];
    
    subplot 211; hold on; grid on; 
    xlabel('t in sec','FontSize', 12);
    ylabel('quaternions','FontSize', 12); 
    title('quaternions over time','FontSize', 14);
    plot(t_span,y_nonlin(:,1:4),'LineWidth',1.5); 
    legend('q0','q1','q2','q3');
    
    subplot 212; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('q_0^2+q^Tq-1','FontSize', 12); 
    title('divergence of q_0^2+q^Tq-1','FontSize', 14);
    plot(t_span,unit,'LineWidth',1.5); 
    
    % omega 
    fig(2) = figure('Name','angular_velocity_over_time'); 
    fig(2).Position = [ 680 476 1159 502];
    
    subplot 311; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('\omega_B','FontSize', 12); 
    title('angular velocity in body frame','FontSize', 14);
    plot(t_span,omega_B(1:2,:)','LineWidth',1.5); 
    legend('\omega_{xB}','\omega_{yB}');
    
    subplot 312; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('\omega_{zB}','FontSize', 12); 
    title('angular velocity \omega_{zB} in body frame','FontSize', 14);
    ylim([floor(min(omega_B(3,:))), ceil(max(omega_B(3,:)))]); 
    plot(t_span,omega_B(3,:)','LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250]);
    
    
    subplot 313; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('\omega_G','FontSize', 12); 
    title('angular velocity in ground frame','FontSize', 14);
    plot(t_span,omega_G','LineWidth',1.5); 
    legend('\omega_{xG}','\omega_{yG}','\omega_{zG}');
    
    % COM in ground frame
    fig(3) = figure('Name','Location_COM_Groud_Frame');
    fig(3).Position = [ 680 476 1159 502];
    
    subplot 311; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('x_G','FontSize', 12); 
    title('x-position COM in ground frame','FontSize', 14);
    plot(t_span,r_G(1,:)','LineWidth',1.5,'Color',[0, 0.4470, 0.7410]); 
    
    subplot 312; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('y_G','FontSize', 12); 
    title('y-position COM in ground frame','FontSize', 14);
    plot(t_span,r_G(2,:)','LineWidth',1.5,'Color',[0.8500, 0.3250, 0.0980]); 
    set(gca,'ydir','reverse');
    
    subplot 313; hold on; grid on; 
    xlabel('t in sec','FontSize', 12); 
    ylabel('z_G','FontSize', 12); 
    title('z-position COM in ground frame','FontSize', 14);
    plot(t_span,r_G(3,:)','LineWidth',1.5,'Color',[0.9290, 0.6940, 0.1250]); 
    
    
    % plot xz location of COM
    fig(4) = figure('Name','XZ_Location_of_COM'); 
    axes; hold on; grid on; 
    xlabel('x [m]','FontSize', 12); 
    ylabel('z [m]','FontSize', 12); 
    title('xz-location of COM','FontSize', 14);
    set(gca,'YAxisLocation','origin','ydir','reverse');
    xlim([-l*1.1 l*1.1]); ylim([-l*1.1 l*1.1]); pbaspect([1 1 1]);
    set(gca,'XAxisLocation','origin','YAxisLocation','origin');
    
    plot(r_G(1,:),r_G(3,:),'LineWidth',1.5);

end

%% save figures
fig_folder = 'figures/';
try
    if SAVE_PLOTS == true
        for ii = 1:N_fig
            print(fig(ii),[fig_folder 'SIM' num2str(SIMULATION) '_' fig(ii).Name ],'-dpng');
        end
    end
catch
    warning("an error occured while saving plots!");    
end

%% animation
if (ANIMATE == true)
    try
        addpath('visualization\');

        animation(R,th,r_G,y_nonlin(:,1:4)',t_span);
    catch
        warning("visualization path could not be found!");
    end
    
end
