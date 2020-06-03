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
Misc.IKfile = {fullfile(DataPath,'trial_1_IK.mot'); fullfile(DataPath,'trial_2_IK.mot'); fullfile(DataPath,'trial_4_IK.mot')};
Misc.IDfile = {fullfile(DataPath,'trial_1_ID.mot'); fullfile(DataPath,'trial_2_ID.mot'); fullfile(DataPath,'trial_4_ID.mot')};
Misc.USfile = {fullfile(DataPath,'trial_1_US.mot'); fullfile(DataPath,'trial_2_US.mot'); fullfile(DataPath,'trial_4_US.mot')}; %
Misc.EMGfile = {fullfile(DataPath,'trial_1_emg.mot'); fullfile(DataPath,'trial_2_emg.mot'); fullfile(DataPath,'trial_4_emg.mot')};
model_path  = fullfile(DataPath,'model.osim');
Out_path    = fullfile(ExamplePath,'Results');                    % folder to store results

% Get start and end time of the different files (you can also specify this
% manually)
time = zeros(size(Misc.IKfile,1),2);
for i = 1:size(Misc.IKfile,1)
    IK = importdata(Misc.IKfile{i});
    time(i,:) = [IK.data(1,1) IK.data(end,1)];
end

%% Settings
Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};    % select the DOFs you want to include in the optimization

% Set the tendon stifness of all muscles
Misc.ATendon = [];      % default way to set tendon stiffenss (default values is 35)

% Settings for estimating tendon stiffness
Misc.Estimate_TendonStifness = {'med_gas_l';'lat_gas_l';'soleus_l'}; % Names of muscles of which tendon stifness is estimated
Misc.lb_kT_scaling = 0.2; % Lower bound for scaling generic tendon stiffness
Misc.ub_kT_scaling = 1.2; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_TendonStifness = {'med_gas_l';'lat_gas_l';'soleus_l'}; % Couple muscles that should have equal tendon stifness
Misc.Coupled_fiber_length = {'med_gas_l';'lat_gas_l'};
Misc.Coupled_slack_length = {}; %{'med_gas_l';'lat_gas_l'};

% Settings for estimating optimal fiber length
Misc.Estimate_OptFL = {'med_gas_l';'lat_gas_l'};%;'lat_gas_l';'soleus_l'}; % Names of muscles of which optimal fiber length is estimated - slack length is estimated for these muscles as well
Misc.lb_lMo_scaling = 0.7; % Lower bound for scaling optimal fiber length
Misc.ub_lMo_scaling = 1.5; % Upper bound for scaling optimal fiber length
Misc.lb_lTs_scaling = 0.7; % Lower bound for scaling tendon slack length
Misc.ub_lTs_scaling = 1.5; % Upper bound for scaling tendon slack length

% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 1;            % Boolean to select US tracking option
Misc.USSelection = {'med_gas_l'};

% Provide the correct headers int case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
Misc.EMGheaders = {'time','med_gas_l', 'soleus_l', 'vas_lat_l'};

% channels you want to use for EMG constraints
Misc.EMGSelection = {'med_gas_l', 'soleus_l'};

% Use this structure if you want to use one EMG channel for multiple
% muscles in the opensim model. The first name of each row is the reference
% and should always be in the header of the EMGfile or in the  EMGheaders.
Misc.EMG_MuscleCopies = {'med_gas_l','lat_gas_l'};       %  use gastrocnemius medialis EMG to constrain activity of the lateral gastrocn

% information for the EMG constraint
Misc.EMGconstr  = 0;     		% Boolean to select EMG constrained option
Misc.EMGbounds  = [-0.3 0.3];  	% upper and lower bound for deviation simulated and measured muscle activity
Misc.MaxScaleEMG = 10; 			% maximal value to scale EMG

% Set weights
Misc.wEMG 		= 0.001;			% weight on tracking EMG
Misc.wlM    = 10;               % weight on tracking fiber length

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 1; 	% TO DO: we should report results of EMG driven simulation as well

% change mesh frequency
Misc.Mesh_Frequency = 100;
%% Run muscle tendon estimator:
[Results,DatStore,Misc] = MuscleTendonEstimator(model_path,time,Out_path,Misc);

% Save the results structure where you want
save('Results.mat','Results');
