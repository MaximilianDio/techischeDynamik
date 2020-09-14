%% pendulum as DAE using redundant coordinates
% close all;
clear; clc;

l1 = 1;
m1 = .1;
l2 = 2;
m2 = 1;
g = 9.81;

alpha0 = pi/2;
beta0 = pi/2;

Dalpha0 = 0;
Dbeta0 = 0;
y0 = [0; l1*sin(alpha0); -l1*cos(alpha0);
      0; l1*sin(alpha0)+l2*sin(beta0); -l1*cos(alpha0)-l2*cos(beta0);
      0; 0; 0;
      0; 0; 0;
      0 ; 0; 0; 0; 0];

p = 2;
ei = 3;
qi = 2;

e = ei*p;
q = qi*p+1;

%%
syms x1 y1 z1 Dx1 Dy1 Dy1 Dz1 x2 y2 z2 Dx2 Dy2 Dy2 Dz2 F1 F2 F3 F4 F5 real

xI = [x1; y1; z1; x2; y2; z2];
xII = [Dx1; Dy1; Dz1; Dx2; Dy2; Dz2];
lambda = [F1; F2; F3; F4; F5];

M = blkdiag(eye(ei)*m1,eye(ei)*m2);
qc = [0;0;0;0;0;0];
qe = [0; 0; -m1*g;0; 0; -m2*g];

Z = eye(e);

%% boundary conditions
c = [x1;
    z1^2+y1^2-l1^2;
    x2;
    (y2-y1)^2+(z2-z1)^2-l2^2;
    z2];
C = jacobian(c,xI);

ct = zeros(e,1);
ctt = [0, 0, 0, 0, 0, 0;
       0, 2*Dy1, 2*Dz1,0,0,0;
       0, 0, 0, 0, 0, 0;
       0, 2*Dy1-2*Dy2, 2*Dz1-2*Dz2, 0, 2*Dy2-2*Dy1, 2*Dz2-2*Dz1;
       0, 0, 0, 0, 0, 0]*xII;

%% DAE formulation
DxI = matlabFunction(Z*xII, 'vars',{[xI;xII;lambda]});
DxII = matlabFunction(M\(qe-qc+C'*lambda),'vars',{[xI;xII;lambda]});
ctt = matlabFunction(ctt, 'vars',{[xI;xII;lambda]});
C = matlabFunction(C,'vars',{[xI;xII;lambda]});
M_ = @(t,x) [eye(e) zeros(e) zeros(e,q); 
           zeros(e) eye(e) zeros(e,q); 
           zeros(q,e) C(x)    zeros(q)];

f = @(t,x) [DxI(x);DxII(x);-ctt(x)];

options = odeset("Mass",M_,"RelTol",1e-4);
[t,y] = ode15s(f,0:0.01:5,y0,options);

figure
plot(y(:,2),y(:,3));
hold on
plot(y(:,5),y(:,6));
%%
figure
subplot 211; plot(t,y(:,1:3));
subplot 212; plot(t,y(:,4:6));
%%
% for ii = 1:2:length(t)
%     cla
%     plot([0,y(ii,2)],[0,y(ii,3)]);
%     hold on;
%     plot([y(ii,2),y(ii,5)],[y(ii,3),y(ii,6)]);
%     
%     t_swift = min(ii-1,40);
%     plot(y(ii-t_swift:ii,2),y(ii-t_swift:ii,3),'k:')
%     plot(y(ii-t_swift:ii,5),y(ii-t_swift:ii,6),'k:')
%     
%     
%     pause(0.01)
%     xlim([-2,4]); ylim([-2,2]);
%     
% end

