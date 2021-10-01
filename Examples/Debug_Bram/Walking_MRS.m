%% Default example solve muscle redundancy 
% (as in DeGroote2016)

clear all;
%% Input information

% select datafolder
ExamplePath = pwd;
DataPath = fullfile(pwd,'WalkingData');

% Add here the paths of IK, ID , US and EMG data trials you want to work with
Misc.IKfile = {fullfile(DataPath,'SimCP-Prosp29-T0_05_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'SimCP-Prosp29-T0_05_ID.sto')};
model_path  = fullfile(DataPath,'simCP-SEMLS-1_pre_PROSP29_scaled03.osim');
Out_path    = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files
time=[5.065 6.125]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)

% Settings
Misc.DofNames_Input={'ankle_flex_r','knee_flex_r','hip_flex_r','hip_add_r','hip_rot_r'};    % select the DOFs you want to include in the optimization

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;

% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;

% name output
Misc.OutName = 'Walking3_';

% adapt the stiffness of the achilles tendon (optional input argument)
% Misc.Set_kT_ByName = {'soleus_r',20;
%     'med_gas_r',20;
%     'lat_gas_r',20};

%% Import CasADi
% import casadi.*

%% Run muscle tendon estimator:
[Results,DatStore] = solveMuscleRedundancy(model_path,time,Out_path,Misc);
