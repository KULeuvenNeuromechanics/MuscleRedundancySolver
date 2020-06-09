%% Example solve muscle redundancy with an MRI model

clear all;
%% Input information
% Install instructions:
%   Add the main path ...\solvemuscleredundancy_dev to your matlab folder,
%   using the following command (in my case)
ExamplePath = pwd;
DataPath = [pwd '\WalkingData'];
% Path up 2 folders is main path
idcs   = strfind(ExamplePath,'\');
newdir = ExamplePath(1:idcs(end)-1);
idcs   = strfind(newdir,'\');
MainPath = newdir(1:idcs(end)-1);
% Add all folders under main path
addpath(genpath(MainPath));

% Add here the paths of IK, ID , US and EMG data trials you want to work with
% As example we use trial 1 2 & 4 
Misc.IKfile = {fullfile(DataPath,'Walking_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'Walking_ID.sto')};
Misc.USfile = {};
Misc.EMGfile = {};

model_path  = fullfile(DataPath,'subject1.osim');
Out_path    = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files (you can also specify this
% manually)
time=[0.516 1.95]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)

%% Settings
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_adduction_r','hip_rotation_r'};    % select the DOFs you want to include in the optimization

% Set the tendon stifness of all muscles
Misc.ATendon = [];      % default way to set tendon stiffenss (default values is 35)

% settings related to parameter estimation
Misc.EMGconstr   = 0;     		% Boolean to select EMG constrained option
Misc.UStracking  = 0;            % Boolean to select US tracking option

% Provide the correct headers int case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
Misc.EMGheaders = {'time','med_gas_l', 'soleus_l', 'vas_lat_l'};

% channels you want to use for EMG constraints
Misc.EMGSelection = {'med_gas_l', 'soleus_l'};

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 0;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 0;
% set the mesh frequency
Misc.Mesh_Frequency = 100;
% name output
Misc.OutName = 'Walking_';

%% Run muscle tendon estimator:
[Results,DatStore] = MuscleTendonEstimator(model_path,time,Out_path,Misc);

% Save the results structure where you want
save('Results.mat','Results');
