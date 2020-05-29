%% Example solve muscle redundancy with an MRI model

% clear variables and commond window
clear all; clc;

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
Misc.IKfile = {fullfile(DataPath,'trial_1_IK.mot')};% fullfile(DataPath,'trial_2_IK.mot'); fullfile(DataPath,'trial_4_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'trial_1_ID.mot')};% fullfile(DataPath,'trial_2_ID.mot'); fullfile(DataPath,'trial_4_ID.mot')};
% Misc.USfile = {fullfile(DataPath,'trial_1_US.mot')};% fullfile(DataPath,'trial_2_US.mot'); fullfile(DataPath,'trial_4_US.mot')}; %
% Misc.EMGfile = {fullfile(DataPath,'trial_1_emg.mot')}; %fullfile(DataPath,'trial_2_emg.mot'); fullfile(DataPath,'trial_4_emg.mot')};
Misc.USfile = {};% fullfile(DataPath,'trial_2_US.mot'); fullfile(DataPath,'trial_4_US.mot')}; %
Misc.EMGfile = {}; %fullfile(DataPath,'trial_2_emg.mot'); fullfile(DataPath,'trial_4_emg.mot')};
model_path  = fullfile(DataPath,'model.osim');
Out_path    = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files (you can also specify this
% manually)
% time = zeros(size(Misc.IKfile,1),2);
% for i = 1:size(Misc.IKfile,1)
%     IK = importdata(Misc.IKfile{i});
%     time(i,:) = [IK.data(1,1) IK.data(end,1)];
% end
time = [0 9999];  % (selects the full file) 

%% Settings
Misc.DofNames_Input={'ankle_angle_l'};    % select the DOFs you want to include in the optimization

Misc.MuscleNames_Input = {'med_gas_l','lat_gas_l','soleus_l','tib_ant_l'}; % select muscles
Misc.Atendon = [];
% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 0;            % Boolean to select US tracking option
Bounds = [];		% currently still empty
% information for the EMG constraint
Misc.EMGconstr  = 0;     		% Boolean to select EMG constrained option
% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 0; 	% TO DO: we should report results of EMG driven simulation as well
% change mesh frequency
Misc.Mesh_Frequency = 100;
% output name
Misc.OutName = 'Minimal_';
%% Run muscle tendon estimator:
[Results,Parameters,DatStore,Misc] = MuscleTendonEstimator(model_path,time,Bounds,Out_path,Misc);

% Save the results structure where you want
save('Results.mat','Results');

%% evaluate results
F = Results.TForce.genericMRS';
dM = DatStore.MAinterp;

Tm = F.*dM;
figure(); plot(Tm); legend(DatStore.MuscleNames);
