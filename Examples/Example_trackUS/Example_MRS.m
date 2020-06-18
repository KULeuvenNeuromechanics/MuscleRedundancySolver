%% Example solve muscle redundancy without parameter estimation

% default muscle redundancy solver as in DeGroote 2016

clear all;

%% Input information
ExamplePath = pwd;
DataPath = [pwd '\Data'];

% Add here the paths of IK, ID trials you want to work with
% As example we use trial 1 2 & 4 
Misc.IKfile = {fullfile(DataPath,'trial_1_IK.mot'), fullfile(DataPath,'trial_2_IK.mot'), fullfile(DataPath,'trial_4_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'trial_1_ID.mot'), fullfile(DataPath,'trial_2_ID.mot'), fullfile(DataPath,'trial_4_ID.mot')};
model_path  = fullfile(DataPath,'model.osim');
Out_path    = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files
time = [0 9999];  % (select all time frames between 0 and 9999) 

% degrees of freedom
Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};    % select the DOFs you want to include in the optimization

% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;

% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(model_path,time,Out_path,Misc);
