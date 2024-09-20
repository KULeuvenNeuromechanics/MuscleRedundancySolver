%% Example EMG driven simulation for ankle, knee and hip in the sagital plane

% In this example we estimate parameters of multiple lower-limb muscles
% using an EMG driven simulation of the ankle-knee and hip.
% Additionally, this example handles different trials with different number
% of degrees of freedom together.
% The same gait trial is used 4 times. First, using all 3 degrees of freedom
% (hip, knee, ankle). Second, using only ankle. Third, using only knee.
% Fourth, using ony hip. The idea is to provide an example to the user on
% how to implement the code with different trials that have different
% number of degrees of freedom.

% clear variables and command window
clear all; clc; close all;

%% Input information
% Add here the paths of IK, ID , US and EMG data trials you want to work with
for t=1:4
    copyfile('./IK_gait.mot',['./IK_gait_copy' num2str(t) '.mot'])
    copyfile('./ID_gait.sto',['./ID_gait_copy' num2str(t) '.sto'])
    copyfile('./EMG_gait.mot',['./EMG_gait_copy' num2str(t) '.mot'])
    Misc.IKfile{t}  = fullfile(pwd,['IK_gait_copy' num2str(t) '.mot']);
    Misc.IDfile{t}  = fullfile(pwd,['ID_gait_copy' num2str(t) '.sto']);
    Misc.EMGfile{t} = fullfile(pwd,['EMG_gait_copy' num2str(t) '.mot']);
    Misc.side{t} = 'l';
    time(t,:) = [1.2 2.3]; 
end
Misc.model_path   = {fullfile(pwd,'ScaledModel.osim')};
Misc.OutPath      = fullfile(pwd,'Results_LowerLimb_EMGconstrained_variableDOF');                    % folder to store results
Misc.AnalysisID = 'v1';

%% Settings

% name of the results file
Misc.OutName ='gait_';

% select degrees of freedom
Misc.DofNames_Input{1,1}={'ankle_angle_l','knee_angle_l','hip_flexion_l'};    % select the DOFs you want to include in the optimization
Misc.DofNames_Input{2,1}={'ankle_angle_l'};    % select the DOFs you want to include in the optimization
Misc.DofNames_Input{3,1}={'knee_angle_l'};    % select the DOFs you want to include in the optimization
Misc.DofNames_Input{4,1}={'hip_flexion_l'};    % select the DOFs you want to include in the optimization

% Set the tendon stifness of all muscles
Misc.kT = [];      % default way to set tendon stiffenss (default values is 35)

% Provide the correct headers int case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
% Misc.EMGheaders = {'Time','bifemlh_r','tib_ant_r','per_long_r','lat_gas_r',...
%     'bifemsh_r','soleus_r','vas_lat_r','vas_med_r','per_brev_l','tib_ant_l','per_long_l',...
%     'lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l',...
%     'rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l','glut_med2_r','rect_fem_r'};
% User can use either Misc.EMGheaders or Misc.EMGFileHeaderCorrespondence.
% Misc.EMGFileHeaderCorrespondence is a n x 2 cell array where first column
% is the header name in the EMG file and second column is the name of the
% muscle in the .osim model corresponding to that header name. The order of
% rows of this matrix does not have to match the order in the EMG file and
% not all rows in this matrix need to be in all the EMG files of all trials
Misc.EMGFileHeaderCorrespondence = {'Time', 'Time';...
    'HamL_r', 'bifemlh_r';...
    'TA_r', 'tib_ant_r';...
    'PerL_r', 'per_long_r';...
    'GL_r', 'lat_gas_r';...
    'HamM_r', 'bifemsh_r';...
    'Sol_r', 'soleus_r';...
    'VL_r', 'vas_lat_r';...
    'VM_r', 'vas_med_r';...
    'PerB_l', 'per_brev_l';...
    'TA_l', 'tib_ant_l';...
    'PerL_l', 'per_long_l';...
    'GL_l', 'lat_gas_l';...
    'GM_l', 'med_gas_l';...
    'Sol_l', 'soleus_l';...
    'VL_l', 'vas_lat_l';...
    'VM_l', 'vas_med_l';...
    'AddL_l','add_long_l';...
    'RF_l', 'rect_fem_l';...
    'TFL_l', 'tfl_l';...
    'GluMed_l', 'glut_med2_l';...
    'HamL_l', 'bifemsh_l';...
    'HamM_l', 'bifemlh_l';...
    'GluMed_r', 'glut_med2_r';...
    'RF_r', 'rect_fem_r'};

% channels you want to use for EMG constraints
Misc.EMGSelection = {'per_brev_l','tib_ant_l','per_long_l','lat_gas_l',...
    'med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l',...
    'rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l'};

%  use vastus lateralis EMG to constrain activity of the vastus intermedius
Misc.EMG_MuscleCopies = {'vas_lat_l','vas_int_l'};       

% parameter estimation
Misc.Estimate_TendonStiffness = {'med_gas_l';'lat_gas_l';'soleus_l';'per_brev_l';'per_long_l'}; % Names of muscles of which tendon stifness is estimated
Misc.lb_kT_scaling = 0.8; % Lower bound for scaling generic tendon stiffness
Misc.ub_kT_scaling = 1.2; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_TendonStiffness = {'med_gas_l','lat_gas_l';...
    'med_gas_l','soleus_l';...
    'per_brev_l','per_long_l'}; % Couple muscles that should have equal tendons stiffness

Misc.Estimate_OptimalFiberLength= {'med_gas_l';'lat_gas_l';'soleus_l';'per_brev_l';'per_long_l'}; % Names of muscles of which tendon stifness is estimated
Misc.lb_lMo_scaling= 0.8; % Lower bound for scaling generic tendon stiffness
Misc.ub_lMo_scaling = 1.2; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_fiber_length = {'med_gas_l','lat_gas_l';...
    'med_gas_l','soleus_l';...
    'per_brev_l','per_long_l'}; % Couple muscles that should have equal tendons stiffness


% select muscles
Misc.MuscleNames_Input = []; % select muscles

% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 0;            % Boolean to select US tracking option

% information for the EMG constraint
Misc.EMGconstr  = 1;     		% Boolean to select EMG constrained option
Misc.EMGbounds  = [-0.1 0.1];  	% upper and lower bound for difference between simulated and measured muscle activity
Misc.BoundsScaleEMG = [0.01 5];  % maximal value to scale EMG

% Set weights
Misc.wEMG   = 0.1;   % weight on tracking EMG
Misc.wAct   = 0.001;
Misc.wTres  = 10;
Misc.wVm    = 0.00001;

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 1; 	% TO DO: we should report results of EMG driven simulation as well

%% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(time,Misc);

for t=1:4
    delete(['./IK_gait_copy' num2str(t) '.mot'])
    delete(['./ID_gait_copy' num2str(t) '.sto'])
    delete(['./EMG_gait_copy' num2str(t) '.mot'])
end
