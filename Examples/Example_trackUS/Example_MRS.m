%% Example solve muscle redundancy with an MRI model

clear all;

%% Input information
% Install instructions:
%   Add the main path ...\solvemuscleredundancy_dev to your matlab folder,
%   using the following command (in my case)
ExamplePath = pwd;
DataPath = [pwd '\Data'];
% Path up 2 folders is main path
idcs   = strfind(ExamplePath,'\');
newdir = ExamplePath(1:idcs(end)-1);
idcs   = strfind(newdir,'\');
MainPath = newdir(1:idcs(end)-1);

% Add all folders under main path
addpath(genpath(MainPath));

% Add here the paths of IK, ID , US and EMG data trials you want to work with
% As example we use trial 1 2 & 4 
Misc.IKfile = {fullfile(DataPath,'trial_1_IK.mot'), fullfile(DataPath,'trial_2_IK.mot'), fullfile(DataPath,'trial_4_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'trial_1_ID.mot'), fullfile(DataPath,'trial_2_ID.mot'), fullfile(DataPath,'trial_4_ID.mot')};
model_path  = fullfile(DataPath,'model.osim');
Out_path    = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files
time = [0 9999];  % (selects the full file) 

% degrees of freedom
Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};    % select the DOFs you want to include in the optimization

% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;

% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(model_path,time,Out_path,Misc);
