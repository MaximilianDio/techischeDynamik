%% MAIN - Task 14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Maximilian Dio - 21595892 - 
% Date:     13.08.2020
% Notes:    change parameters of crankshaft in Crankshaft class
%           To animate the results set the ANIMATE flag true
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% clean workspace and cmd window
clear; clc;

%% FLAGS
ANIMATE = 0; % (0/1) 1: animate, 0: not 

%% PATHS
% add paths
addpath("subtasks");
addpath("visualization");
addpath("solver");

%% Simulation - run subtasks

% -- subtask d: simulation as DAE system
% -- subtask e1: simulation via manual coordinate partitioning by resolving
% algebraic connection
% -- subtask e2: simulation via manual coordinate partitioning by integration
% -- subtask f: simulation via coordinate partitioning by QR decomposition

alpha0 = 0.1;
Dalpha0 = 0.1;

Tend = 10;

crankshafts = { 
            CrankshaftTreeDAE(alpha0,Dalpha0);
            CrankshaftTreeMPAnalytical(alpha0,Dalpha0);
            CrankshaftTreeMPIntegration(alpha0,Dalpha0);
            CrankshaftTreeQR(alpha0,Dalpha0);
            CrankshaftRedundantDAE(alpha0,Dalpha0);
            CrankshaftRedundantQR(alpha0,Dalpha0);
            };

results = cell(size(crankshafts));
for ii = 1:length(crankshafts)
    results{ii} = crankshafts{ii}.solve(Tend);
end

%% visualize results
close all;
vis_results(crankshafts,results,'animate',ANIMATE);
