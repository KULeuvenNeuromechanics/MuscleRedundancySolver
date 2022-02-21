%% Example EMG driven simulation for ankle, knee and hip in the sagital plane

% In this example we estimate parameters of multiple lower-limb muscles
% using an EMG driven simulation of the ankle-knee and hip.

% Clear variables and command window
clear all; clc; close all;

%% Input information

% select datafolder
AnalysisPath = pwd;
MainPath = pwd;
DataPath = fullfile(MainPath,'Data');

Misc.subjectName = 'CP2';
model_path  = fullfile(DataPath,'Model','C2_CP2_T0_scaled_sf.osim');
OutPath    = fullfile(AnalysisPath,'Results_3');                    % folder to store results
Misc.model_path = model_path;
Misc.OutPath = OutPath;
Misc.MRS_validation_together = 0;

trialTypes = {'s2s','walking','squat','cmj'};
trialInfo = readtable(fullfile(DataPath,'Timepoints.xlsx'));
sides = {'L','R'};
trialInfo = trialInfo(ismember(trialInfo.Var3,trialTypes),:);
trialInfo = trialInfo(ismember(trialInfo.Var4,sides),:);
% Get start and end time of the different files
time = [trialInfo.Var5 trialInfo.Var6];
% time=[0.516 1.95]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)
Misc.trialName = trialInfo.Var2';
Misc.IKfile = strcat({fullfile(DataPath,'IK\')},trialInfo.Var2,{'_IK.mot'})';
Misc.IDfile = strcat({fullfile(DataPath,'ID\')},trialInfo.Var2,{'_IK_ID.sto'})';
Misc.EMGfile = cell(1,length(Misc.IKfile));
Misc.side = lower(trialInfo.Var4)';
Misc.trialType = trialInfo.Var3';
Misc.trialFP = trialInfo.Var7;

% select the DOFs you want to include in the optimization
% Add IPSA conditions
trialTypes_nonIPSA = {'s2s','walking','squat','cmj'};
DOF_nonIPSA = {'hip_flexion_','hip_adduction_','hip_rotation_','knee_angle_','ankle_angle_'};
for t=1:length(Misc.trialType)
    if ismember(Misc.trialType{t},trialTypes_nonIPSA)
        Misc.DofNames_Input{t}=strcat(DOF_nonIPSA,Misc.side{t});
    end
end

% Name of the results file
Misc.OutName = strcat(trialInfo.Var2,{'_'},trialInfo.Var3,{'_'},trialInfo.Var4,{'_'})';

%% Settings

% Set the tendon stifness of all muscles
Misc.kT = [];      % default way to set tendon stiffenss (default values is 35)

% In case the headers in your EMG file differ from the muscle names in your OpenSim model, 
% assign here the correct names to your EMG file (leave empty when you don't want to use this)
% Misc.EMGheaders = {'Time','bifemlh_r','tib_ant_r','per_long_r','lat_gas_r','bifemsh_r','soleus_r','vas_lat_r','vas_med_r','per_brev_l','tib_ant_l','per_long_l',...
%     'lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l','rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l','glut_med2_r','rect_fem_r'};
Misc.EMGFileHeaderCorrespondence = {'Time', 'Time';...
    'RREF', 'rect_fem_r';...
    'RVAL', 'vas_lat_r';...
    'RBIF', 'bifemlh_r';...
    'RMEH', 'semiten_r';...
    'RTIA', 'tib_ant_r';...
    'RGAS', 'lat_gas_r';...
    'RSOL', 'soleus_r';...
    'RGLU', 'glut_med1_r';...
    'LREF', 'rect_fem_l';...
    'LVAL', 'vas_lat_l';...
    'LBIF', 'bifemlh_l';...
    'LMEH', 'semiten_l';...
    'LTIA', 'tib_ant_l';...
    'LGAS', 'lat_gas_l';...
    'LSOL', 'soleus_l';...
    'LGLU', 'glut_med1_l'};

% Channels you want to use for EMG constraints
Misc.EMGSelection = {'rect_fem_r','vas_lat_r','bifemlh_r','semiten_r','tib_ant_r','lat_gas_r','soleus_r','glut_med1_r',...
    'rect_fem_l','vas_lat_l','bifemlh_l','semiten_l','tib_ant_l','lat_gas_l','soleus_l','glut_med1_l'};

Misc.EMG_MuscleCopies = {'lat_gas_r','med_gas_r';...
    'vas_lat_r',   'vas_med_r';...
    'vas_lat_r',   'vas_int_r';...
    'glut_med1_r', 'glut_med2_r';...
    'glut_med1_r', 'glut_med3_r';...
    'lat_gas_l',   'med_gas_l';...
    'vas_lat_l',   'vas_med_l';...
    'vas_lat_l',   'vas_int_l';...
    'glut_med1_l', 'glut_med2_l';...
    'glut_med1_l', 'glut_med3_l'};
    
% % Select muscles
% Misc.MuscleNames_Input = Misc.EMGSelection; % select muscles

% NOTE: The coupling here means that the parameters changes by the same
% factor, not that the value itself is coupled
Misc.Coupled_fiber_length = {'lat_gas_r', 'med_gas_r';...
    'lat_gas_l', 'med_gas_l'};
%     'vas_med_r', 'vas_lat_r';...
%     'vas_lat_r', 'vas_int_r'};%...
    % 'glut_med1_r','glut_med2_r';...
    % 'glut_med1_r','glut_med3_r';...
    % 'glut_max1_r','glut_max2_r';...
    % 'glut_max1_r','glut_max3_r';...
    % 'glut_min1_r','glut_min2_r';...
    % 'glut_min1_r','glut_min3_r'}; % Couple muscles that should have equal optimal fiber length
Misc.Coupled_slack_length = {'lat_gas_r', 'med_gas_r';...
    'lat_gas_l', 'med_gas_l'};
%     'vas_med_r', 'vas_lat_r';...
%     'vas_lat_r', 'vas_int_r'};%...
    % 'glut_med1_r','glut_med2_r';...
    % 'glut_med1_r','glut_med3_r';...
    % 'glut_max1_r','glut_max2_r';...
    % 'glut_max1_r','glut_max3_r';...
    % 'glut_min1_r','glut_min2_r';...
    % 'glut_min1_r','glut_min3_r'}; % Couple muscles that should have equal tendon slack length

% Settings for estimating optimal fiber length
Misc.lb_lMo_scaling = 0.5; % Lower bound for scaling optimal fiber length
Misc.ub_lMo_scaling = 1.5; % Upper bound for scaling optimal fiber length
Misc.lb_lTs_scaling = 0.9; % Lower bound for scaling tendon slack length
Misc.ub_lTs_scaling = 1.1; % Upper bound for scaling tendon slack length

% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 0;            % Boolean to select US tracking option

% Information for the EMG constraint
Misc.EMGconstr  = 1;     		% Boolean to select EMG constrained option
Misc.EMGfile = strcat({fullfile(DataPath,'EMG\EMG_norm_')},trialInfo.Var2,{'.mot'})';
Misc.EMGbounds  = [-0.1 0.1];  	% upper and lower bound for difference between simulated and measured muscle activity
Misc.BoundsScaleEMG = [0.5 1.1];  % maximal value to scale EMG

% Set weights
Misc.wEMG   = 0.5;   % weight on tracking EMG
Misc.wAct   = 1;
Misc.wTres  = 1000;
Misc.wVm    = 0.001;

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 1; 	% TO DO: we should report results of EMG driven simulation as well

Misc = getMuscleProperties(model_path,Misc); % NOTE: getShift was being used before adjusting kT, it has been corrected here

% Settings for estimating optimal fiber length
Misc.Estimate_OptimalFiberLength = Misc.allMuscleList(1:86);% Manually removed trunk muscles

% Settings for estimating tendon slack length
Misc.Estimate_TendonSlackLength = Misc.allMuscleList(1:86);% Manually removed trunk muscles

% Settings for estimating tendon stiffness
Misc.Estimate_TendonStiffness = {'rect_fem_r','vas_lat_r','bifemlh_r','semiten_r','tib_ant_r','lat_gas_r','soleus_r','glut_med1_r',...
    'rect_fem_l','vas_lat_l','bifemlh_l','semiten_l','tib_ant_l','lat_gas_l','soleus_l','glut_med1_l'}; % Names of muscles of which tendon stifness is estimated
Misc.lb_kT_scaling = 0.1; % Lower bound for scaling generic tendon stiffness
Misc.ub_kT_scaling = 2.2; % Upper bound for scaling generic tendon stiffness
% NOTE: The coupling here means that the parameters changes by the same
% factor, not that the value itself is coupled
Misc.Coupled_TendonStiffness = {'lat_gas_r', 'med_gas_r';...
    'lat_gas_l', 'med_gas_l'};
%     'vas_lat_r', 'vas_med_r';...
%     'vas_lat_l', 'vas_med_l'} % Couple muscles that should have equal tendon stifness

% For testing purposes
Misc.Mesh_Frequency=20;

% Trials to analyze
% Misc.trials_sel = 1:length(Misc.IKfile);
% Misc.trials_sel = find(ismember(Misc.side,'l'));
% Misc.trials_sel = find(ismember(Misc.side,'r'));
Misc.trials_sel = [12 14 15 16];
Misc.analysisName = 'trials_12_14_15_16';

%% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(time,Misc);
