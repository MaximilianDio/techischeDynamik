%% MAIN: start here - Task 14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Maximilian Dio - 21595892 - 
% Date:     15.09.2020
% Notes:    To animate the results set the ANIMATE flag true parameters can
%           be changed in the Crankshaft class
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean workspace and cmd window
clear; clc;

%% FLAGS
ANIMATE = 0; % (0/1) 1: animate, 0: not
SAVE = 0; % (0/1) 1: save plots in figures folder, 0: not

%% select simulation parameters
% initial values (beta and Dbeta will be calculated to meet constraints)
alpha0 = 0.1;
Dalpha0 = 0.1;

% select integration time
% choose longer integration time to see divergence of boundary conditions
% -> use ode23t for DAEsolver when integration interval is long 
% - solver capitulates after about 20 seconds
Tend = 2; 

% timestep and integration tolerances
% tolerances below 10^-4 return unusable results over time!
relTol = 1e-4;
absTol = 1e-4;
tStep  = 0.01; 

% choose solver 
ode_solver = @ode45; % choose ode45, ode23 or ode113
dae_solver = @ode15s; % choose ode15s or ode23t

%% PATHS
% add paths
addpath("subtasks");
addpath("visualization");
addpath("solver");

%% RUN  
% Simulation - run subtasks

% -- subtask d:  simulation as DAE system
% -- subtask e1: simulation via manual coordinate partitioning by resolving
%                algebraic connection
% -- subtask e2: simulation via manual coordinate partitioning by integration
% -- subtask f:  simulation via coordinate partitioning by QR decomposition
%--------------------------------------------------------------------------
% -- subtask i:  simulation as DAE with redundant coordinates
% -- subtask j:  simulation via coordinate partitioning by QR decomposition
%                using redundant coordinates

% generate instances of crankshafts with different solution methods
% (subtasks)
crankshafts = { 
            CrankshaftTreeDAE(alpha0,Dalpha0,dae_solver);
            CrankshaftTreeMPAnalytical(alpha0,Dalpha0,ode_solver);
            CrankshaftTreeMPIntegration(alpha0,Dalpha0,ode_solver);
            CrankshaftTreeQR(alpha0,Dalpha0,ode_solver);
            CrankshaftRedundantDAE(alpha0,Dalpha0,dae_solver);
            CrankshaftRedundantQR(alpha0,Dalpha0,ode_solver);
};

% preallocate result cell
results = cell(size(crankshafts));

% run simulation
for ii = 1:length(crankshafts)
    crankshafts{ii}.relTol = relTol;
    crankshafts{ii}.absTol = absTol;
    crankshafts{ii}.tStep = tStep;
    results{ii} = crankshafts{ii}.solve(Tend);
end

%% VISUALIZE
close all;
vis_results(crankshafts,results,'animate',ANIMATE,'save',SAVE,'text',string(Tend) + "-" + string(relTol) + "-" + string(absTol) + "-" + string(tStep)  );
