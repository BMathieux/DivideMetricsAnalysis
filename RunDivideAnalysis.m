%% RunDivideAnalysis.m
% -------------------------------------------------------------------------
% Drainage Divide Metrics Analysis
%
% This script performs a full analysis of drainage divide dynamics using 
% Topotoolbox. It extracts hillslope metrics, delineates paired basins, 
% computes hilltop curvature-derived denudation rates, and evaluates divide 
% stability and migration using morphometric and fluvial parameters.
%
% Workflow:
% 1. Prepare DEMs
% 2. Define parameters
% 3. Select divide and extract neighboring basins
% 4. Hilltop nodes filtering and measure hilltop curvature (optional,
% RunCHt)
% 5. Select window size for hilltop curvature (optional, iWs) 
% 6. Compute cross-divide metrics and create tables
% 7. Save tables as .csv
%
% Date: May 11, 2025
% Author: Bastien Mathieux (bastien.mathieux[at]gmail.com)


% Interactive options (true for interactive selection, false for default values set in DefineParams.m)
iA = true;   % Set Amin interactively or use default value
iWs = true;  % Set Ws interactively or use default value

% Toggle hilltop curvature analysis (CHt measurements, associated
% denudation rates and divide migration rates)

RunCHt = true;

%% Prepare DEMs

DEM1 = 'Thur_1m.tif';   % High-resolution DEM for hilltop curvature, basin extraction and basin metrics
DEMg = 'SVosges_10m.tif'; % Global DEM for ksn and chi calculation (can be the same than DEM1)

fprintf('Loading and preparing DEMs...\n'); tic
info = georasterinfo(DEM1); % Save spatial reference
DEM1 = GRIDobj(DEM1); 
DEMg = GRIDobj(DEMg); 

toc

%% Step 1: Define parameters 
params = DefineParams();

%% Step 2: Select divide and extract neighboring basins
fprintf('Selecting divide and extracting neighboring basins...\n');
BData = extractBasins(DEM1, params, 'Verbose', true, 'InteractiveAmin', iA);

%% Step 3: Select hilltop nodes and extract hilltop curvature
if RunCHt
    fprintf('Extracting hilltop metrics...\n');
    CHt_data = CHt_metrics(BData, DEM1, params, 'UseParallel', false, 'Verbose', true);
else
    CHt_data = nan;
end

%% Step 4: Select window size (ws) for hilltop curvature
fprintf('Selecting window size (ws) for curvature calculation...\n');
if RunCHt, params.wl = SelectWs(CHt_data, iWs, params.wl); end

%% Step 5: Compute cross-divide metrics and create tables

fprintf('Computing cross-divide metrics...\n');
[DTable, MTable] = CalcDivideTable(CHt_data, BData, DEMg, info, params, ...
    'Parallel', false, 'Verbose', true, 'RunCHt',RunCHt);

%% Step 5: Display and save the results
disp('Analysis complete.');
writetable(DTable, 'DTable.csv');

writetable(MTable, 'MTable.csv');
