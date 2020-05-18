%%
close all; clear; clc;

%%
syms alpha beta gamma l m g a b theta omega_1B omega_2B omega_3B real

SGK_ = [1 0 0; 0 cos(alpha) sin(alpha); 0 -sin(alpha) cos(alpha)];

SK_K__ = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];

SK__B = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];

SGB = SGK_*SK_K__*SK__B;
SBG = SGB';

F_G = [0; m*g; 0];
r_B = [0;0;l];

% torques
tau_B = cross( r_B,SBG*F_G)
tau_G = simplify(SGB*tau_B)

% moment of inertia
Theta_B = m*[a 0 0; 0 a 0; 0 0 b]

% 
omega_B = [omega_1B; omega_2B; omega_3B];

simplify(SGB*(Theta_B\tau_B-subs(Theta_B\(cross(omega_B,Theta_B*omega_B)),(a-b)/a,theta)))

