%% Example solve muscle redundancy with an MRI model

clear all;
%% Input information
% NEED TO ADD TO PATH HOW EVER YOU DO IT 

ID_DataPath         = [pwd '\ID'];
IK_DataPath         = [pwd '\IK'];
EMG_DataPath        = [pwd '\EMG'];
Models_DataPath     = [pwd '\Models'];
Results_DataPath    = [pwd '\Results2'];

% Add here the paths of IK, ID , US and EMG data trials you want to work with
% As example we use trial 1 2 & 4 
Misc.IKfile = {fullfile(IK_DataPath,'ng10203_IK.mot'); fullfile(IK_DataPath,'ng10220_IK.mot')};% fullfile(IK_DataPath,'ng10217_IK.mot'); fullfile(IK_DataPath,'ng10220_IK.mot'); fullfile(IK_DataPath,'ng10224_IK.mot')};
Misc.IDfile = {fullfile(ID_DataPath,'ng10203_ID.sto'); fullfile(ID_DataPath,'ng10220_ID.sto')};% fullfile(ID_DataPath,'ng10217_ID.sto'); fullfile(ID_DataPath,'ng10220_ID.sto'); fullfile(ID_DataPath,'ng10224_ID.sto')};
% Misc.USfile = {fullfile(DataPath,'trial_1_US.mot')};% fullfile(DataPath,'trial_2_US.mot'); fullfile(DataPath,'trial_4_US.mot')}; %
Misc.USfile = [];
Misc.EMGfile = {fullfile(EMG_DataPath,'ng10203_crop_emg.mot'); fullfile(EMG_DataPath,'ng10220_crop_emg.mot')};% fullfile(EMG_DataPath,'ng10217_emg.mot'); fullfile(EMG_DataPath,'ng10220_emg.mot'); fullfile(EMG_DataPath,'ng10224_emg.mot')};

% model_path  = fullfile(Models_DataPath,'gait2392_arms_ORLAU_scaled.osim');
model_path  = fullfile(Models_DataPath,'model.osim');

Out_path    = fullfile(Results_DataPath);                    % folder to store results

% Get start and end time of the different files 
% time = zeros(size(Misc.IKfile,1),2);
% for i = 1:size(Misc.IKfile,1)
%     IK = importdata(Misc.IKfile{i});
%     time(i,:) = [IK.data(1,1) IK.data(end,1)];
% end
% time = [1.23 2.3; 
%     4.5 6.2];
time = [1.2 2; 
    4.5 5];
%% Settings
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r'};%,'hip_flexion_l','hip_adduction_l','hip_rotation_l'};    % select the DOFs you want to include in the optimization

% Set the tendon stifness of all muscles
Misc.ATendon = [];      % default way to set tendon stiffenss (default values is 35)

% Settings for estimating tendon stiffness
Misc.Estimate_TendonStifness = {'med_gas_r';'lat_gas_r';'soleus_r'}; % Names of muscles of which tendon stifness is estimated
Misc.lb_kT_scaling = 0.1; % Lower bound for scaling generic tendon stiffness
Misc.ub_kT_scaling = 2.2; % Upper bound for scaling generic tendon stiffness
Misc.Coupled_TendonStifness = {'med_gas_r';'lat_gas_r';'soleus_r'}; % Couple muscles that should have equal tendon stifness

% Settings for estimating optimal fiber length
Misc.Estimate_OptFL = {'med_gas_r';'soleus_r';'lat_gas_r'}; % Names of muscles of which optimal fiber length is estimated - slack length is estimated for these muscles as well
Misc.lb_lMo_scaling = 0.1; % Lower bound for scaling optimal fiber length
Misc.ub_lMo_scaling = 2.2; % Upper bound for scaling optimal fiber length
Misc.lb_lTs_scaling = 0.9; % Lower bound for scaling tendon slack length
Misc.ub_lTs_scaling = 1.1; % Upper bound for scaling tendon slack length
Misc.Coupled_fiber_length = {'med_gas_r';'lat_gas_r'}; % Couple muscles that should have equal optimal fiber length
Misc.Coupled_slack_length = {'med_gas_r';'lat_gas_r'}; % Couple muscles that should have equal tendon slack length

% Select muscle for which you want the fiberlengths to track the US data
Misc.UStracking  = 0;            % Boolean to select US tracking option
% Misc.USSelection = {'med_gas_l'};

% Provide the correct headers int case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
% Misc.EMGheaders = {'time','soleus_r', 'rect_fem_r',  'lat_gas_r', 'med_gas_r', 'vas_lat_r', 'vas_med_r', 'tib_ant_r'};
Misc.EMGheaders = {};

% channels you want to use for EMG constraints
Misc.EMGSelection = {'soleus_r', 'rect_fem_r',  'lat_gas_r', 'med_gas_r', 'vas_lat_r', 'vas_med_r', 'bifemlh_r', 'semimem_r', 'tib_ant_r', 'per_brev_r', 'ext_dig_r', 'ext_hal_r'};

% Use this structure if you want to use one EMG channel for multiple
% muscles in the opensim model. The first name of each row is the reference
% and should always be in the header of the EMGfile or in the  EMGheaders.
% Need to make clear to only include muscles in Misc.EMGSelection
Misc.EMG_MuscleCopies = {'semimem_r','semiten_r'};       %  use gastrocnemius medialis EMG to constrain activity of the lateral gastrocn
% Misc.EMG_MuscleCopies = [];

% information for the EMG constraint
Misc.EMGconstr  = 1;     		% Boolean to select EMG constrained option
Misc.EMGbounds  = [-0.3 0.3];  	% upper and lower bound for deviation simulated and measured muscle activity
Misc.BoundsScaleEMG = [0.8 1.2]; 			% maximal value to scale EMG

% Set weights
Misc.wEMG   = 100;   % weight on tracking EMG
Misc.wlM    = 10;   % weight on tracking fiber length

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 1; 	% TO DO: we should report results of EMG driven simulation as well

%% Run muscle tendon estimator:
[Results,DatStore,Misc] = MuscleTendonEstimator(model_path,time,Out_path,Misc);

% Save the results structure where you want
save('Results.mat','Results');
