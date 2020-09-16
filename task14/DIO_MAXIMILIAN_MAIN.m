%% MAIN: start here - Task 14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Maximilian Dio - 21595892 - 
% Date:     15.09.2020
% Notes:    To animate the results set the ANIMATE flag true parameters can
% be changed in the Crankshaft class
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% clean workspace and cmd window
clear; clc;

%% FLAGS
ANIMATE = 0; % (0/1) 1: animate, 0: not
SAVE = 0; % (0/1) 1: save plots in figures folder, 0: not

%% PATHS
% add paths
addpath("subtasks");
addpath("visualization");
addpath("solver");

%% RUN  
% Simulation - run subtasks

% -- subtask d: simulation as DAE system
% -- subtask e1: simulation via manual coordinate partitioning by resolving
% algebraic connection
% -- subtask e2: simulation via manual coordinate partitioning by integration
% -- subtask f: simulation via coordinate partitioning by QR decomposition
% -- subtask i: simulation as DAE with redundant coordinates
% -- subtask j: simulation via coordinate partitioning by QR decomposition
% using redundant coordinates

alpha0 = 0.1;
Dalpha0 = 0.1;

Tend = 2;

crankshafts = { 
            CrankshaftTreeDAE(alpha0,Dalpha0);
            CrankshaftTreeMPAnalytical(alpha0,Dalpha0);
            CrankshaftTreeMPIntegration(alpha0,Dalpha0);
            CrankshaftTreeQR(alpha0,Dalpha0);
            CrankshaftRedundantDAE(alpha0,Dalpha0);
            CrankshaftRedundantQR(alpha0,Dalpha0);
            };

results = cell(size(crankshafts));

% run simulation
for ii = 1:length(crankshafts)
    results{ii} = crankshafts{ii}.solve(Tend);
end

relTol = Crankshaft.relTol;

%% VISUALIZE
close all;
vis_results(crankshafts,results,'animate',ANIMATE,'save',SAVE,'text',string(Tend) + "-" + relTol);
