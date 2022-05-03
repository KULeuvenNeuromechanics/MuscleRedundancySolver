%% Default example solve muscle redundancy 
% (as in DeGroote2016)

clear all; clc; close all;
%% Input information

% select datafolder
ExamplePath = pwd;
DataPath = [pwd '\WalkingData'];

% Add here the paths of IK, ID , US and EMG data trials you want to work with
Misc.IKfile = {fullfile(DataPath,'Walking_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'Walking_ID.sto')};
Misc.model_path  = fullfile(DataPath,'subject1.osim');
Misc.OutPath     = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files
time=[0.516 1.95]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)

% Settings
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};    % select the DOFs you want to include in the optimization

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = true;

% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;

% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 1; 	% TO DO: we should report results of EMG driven simulation as well

% name output
Misc.OutName = 'Walking3_';

% adapt the stiffness of the achilles tendon (optional input argument)
Misc.Set_kT_ByName = {'soleus_r',20;
    'med_gas_r',20;
    'lat_gas_r',20};

%% Run muscle tendon estimator:
[Results,DatStore] = solveMuscleRedundancy(time,Misc);
