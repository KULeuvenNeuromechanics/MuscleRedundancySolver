%% Example solve muscle redundancy with an MRI model

% clear variables and commond window
clear all; clc;

%% Input information
% Add here the paths of IK, ID , US and EMG data trials you want to work with
% As example we use trial 1 2 & 4 
Misc.IKfile  = {fullfile(pwd,'IK_gait.mot')};
Misc.IDfile  = {fullfile(pwd,'ID_gait.sto')};
Misc.EMGfile = {fullfile(pwd,'EMG_gait.mot')};
model_path   = {fullfile(pwd,'ScaledModel.osim')};
Out_path     = fullfile(pwd,'Results_LowerLimb');                    % folder to store results
time = [1.2 2.3]; 

%% Settings

% name of the resuls file
Misc.OutName ='gait_';

% select degrees of freedom
Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l'};    % select the DOFs you want to include in the optimization

% Set the tendon stifness of all muscles
Misc.kT = [];      % default way to set tendon stiffenss (default values is 35)

% Provide the correct headers int case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
Misc.EMGheaders = {'Time','bifemlh_r','tib_ant_r','per_long_r','lat_gas_r','bifemsh_r','soleus_r','vas_lat_r','vas_med_r','per_brev_l','tib_ant_l','per_long_l',...
    'lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l','rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l','glut_med2_r','rect_fem_r'};

% channels you want to use for EMG constraints
Misc.EMGSelection = {'per_brev_l','tib_ant_l','per_long_l','lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l','rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l'};

Misc.EMG_MuscleCopies = {};       %  use gastrocnemius medialis EMG to constrain activity of the lateral gastrocn

% select muscles
Misc.MuscleNames_Input = Misc.EMGSelection; % select muscles

% Settings for estimating tendon stiffness
Misc.Estimate_TendonStiffness = Misc.EMGSelection; % Names of muscles of which tendon stifness is estimated
Misc.lb_kT_scaling = 0.1; % Lower bound for scaling generic tendon stiffness
Misc.ub_kT_scaling = 2.2; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_TendonStiffness = Misc.EMGSelection; % Couple muscles that should have equal tendon stifness

% Settings for estimating optimal fiber length
Misc.Estimate_OptimalFiberLength = Misc.EMGSelection; % Names of muscles of which optimal fiber length is estimated - slack length is estimated for these muscles as well
Misc.lb_lMo_scaling = 0.1; % Lower bound for scaling optimal fiber length
Misc.ub_lMo_scaling = 2.2; % Upper bound for scaling optimal fiber length
Misc.lb_lTs_scaling = 0.9; % Lower bound for scaling tendon slack length
Misc.ub_lTs_scaling = 1.1; % Upper bound for scaling tendon slack length
Misc.Coupled_fiber_length = {'med_gas_l';'lat_gas_l'}; % Couple muscles that should have equal optimal fiber length
Misc.Coupled_slack_length = {'med_gas_l';'lat_gas_l'}; % Couple muscles that should have equal tendon slack length

% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 0;            % Boolean to select US tracking option

% information for the EMG constraint
Misc.EMGconstr  = 1;     		% Boolean to select EMG constrained option
Misc.EMGbounds  = [-0.01 0.01];  	% upper and lower bound for difference between simulated and measured muscle activity
Misc.BoundsScaleEMG = [0.01 3];  % maximal value to scale EMG

% Set weights
Misc.wEMG   = 0.1;   % weight on tracking EMG
Misc.wAct   = 0.1;
Misc.wTres  = 1000;
Misc.wVm    = 0.001;

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 1; 	% TO DO: we should report results of EMG driven simulation as well

%% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(model_path,time,Out_path,Misc);
