%% Example solve muscle redundancy with an MRI model

% clear variables and commond window
clear all; clc;

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
Misc.IKfile = {fullfile(DataPath,'Walking_IK.mot')};% fullfile(DataPath,'trial_2_IK.mot'); fullfile(DataPath,'trial_4_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'Walking_ID.sto')};% fullfile(DataPath,'trial_2_ID.mot'); fullfile(DataPath,'trial_4_ID.mot')};
Misc.USfile = {};% fullfile(DataPath,'trial_2_US.mot'); fullfile(DataPath,'trial_4_US.mot')}; %
Misc.EMGfile = {}; %fullfile(DataPath,'trial_2_emg.mot'); fullfile(DataPath,'trial_4_emg.mot')};

model_path  = fullfile(DataPath,'subject1.osim');
Out_path    = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files (you can also specify this
% manually)
time = zeros(size(Misc.IKfile,1),2);
for i = 1:size(Misc.IKfile,1)
    IK = importdata(Misc.IKfile{i});
    time(i,:) = [IK.data(1,1) IK.data(end,1)];
end

%% Settings
Misc.DofNames_Input={'ankle_angle_l'};%,'hip_flexion_l','hip_adduction_l','hip_rotation_l'};    % select the DOFs you want to include in the optimization

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
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 0; 	% TO DO: we should report results of EMG driven simulation as well
% set the mesh frequency
Misc.Mesh_Frequency = 200;
% bounds
Bounds = [];
% name output
Misc.OutName = 'Walking_';

%% Run muscle tendon estimator:
[Results,Parameters,DatStore] = MuscleTendonEstimator(model_path,time,Bounds,Out_path,Misc);

% Save the results structure where you want
save('Results.mat','Results');
