%% MAIN - Task 14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author:   Maximilian Dio - 21595892 - 
% Date:     13.08.2020
% Notes:    change parameters of crankshaft in crankshaft.m
%           To animate the results set the ANIMATE flag true
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% clean figures, workspace and cmd window
close all; clear; clc;

%% FLAGS
ANIMATE = false;

%% PATHS
% add paths
addpath("subtasks");
addpath("visualization");
%%
% load system
syms alpha beta Dalpha Dbeta real
cs = crankshaft(alpha,beta,Dalpha,Dbeta);

%% Simulation - run subtasks

% -- subtask d: simulation as DAE system
% -- subtask e1: simulation via manual coordinate partitioning by resolving
% algebraic connection
% -- subtask e2: simulation via manual coordinate partitioning by integration
% -- subtask f: simulation via coordinate partitioning by QR decomposition

subtasks = { 
%             @(cs) subtask_d(cs);
            @(cs) subtask_e1(cs);
            @(cs) subtask_e2(cs);
%             @(cs) subtask_f(cs)
            };

results = cell(size(subtasks));
for ii = 1:length(subtasks)
%     try
    [results{ii}.y, results{ii}.Dy, results{ii}.c, results{ii}.Dc ,results{ii}.task_info] = subtasks{ii}(cs);
%     catch
%     warning
%     end
end

%% visualize results
vis_results(cs,results,'animate',ANIMATE);
