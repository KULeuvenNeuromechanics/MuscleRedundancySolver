%% Example solve muscle redundancy with an MRI model
clear all; clc;
% Install instructions:
%   Add the main path ...\solvemuscleredundancy_dev to your matlab folder,
%   using the following command (in my case)
MainPath = 'C:\Users\u0088756\Documents\PostDoc\Software\Original\solvemuscleredundancy_dev';
addpath(genpath(MainPath));

% Input arguments
Datapath =fullfile(MainPath,'Examples','EMG_constraintMRI');
IK_path = fullfile(Datapath,'gait1_kinematics.mot');      % point to the IK file
ID_path = fullfile(Datapath,'inverse_dynamics.sto');      % point to the ID file
model_path = fullfile(Datapath,'SEMLS_2_gen_final_2392_Fmax.osim');         % point to the model
Out_path=fullfile(MainPath,'Results_MRI');                    % folder to store results

% Example for Hans: Note if you want to analyse all the data in the IK file
IK = importdata(IK_path);
time = [IK.data(1,1) IK.data(end,1)];

% settings
Misc.DofNames_Input={'ankle_flex_r','knee_flex_r'};    % select the DOFs you want to include in the optimization
% Misc.DofNames_Input={'ankle_flex_r','knee_flex_r','hip_flex_r','hip_add_r','hip_rot_r'};    % select the DOFs you want to include in the optimization
Misc.RunAnalysis = 0;   % boolean to select if you want to run the muscle analysis

% settings related to tendon stiffness
Misc.ATendon = [];      % default way to set tendon stiffenss (default values is 35)
Misc.Set_ATendon_ByName = {'gas_med_r',15;
    'gas_lat_r',15;
    'soleus_r',15};

% information for the EMG constraint
Misc.EMGconstr  = 1;     % Boolean to select EMG constrained option
Misc.EMGfile    = fullfile(Datapath,'gait1_EMG.mot');
Misc.EMGbounds  = [-0.3 0.3];    % upper and lower bound for deviation simulated and measured muscle activity
Misc.MaxScale   = 10;  % maximal value to scale EMG 
Misc.ActDynEMG = 1;     % select if you want to process the EMG through activation dynamics

% Provide the correct headers int case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
Misc.EMGheaders = {'time','rectus_fem_l', 'vas_lat_l', 'bi_fem_lh_l', 'semiten_l', 'tib_ant_l', 'gas_med_l', 'soleus_l', 'glut_med2_l',...
    'rectus_fem_r', 'vas_lat_r', 'bi_fem_lh_r', 'semiten_r', 'tib_ant_r', 'gas_med_r', 'soleus_r', 'glut_med2_r'}; 

% channels you want to use for EMG constraints
Misc.EMGSelection = {'rectus_fem_r', 'vas_lat_r', 'bi_fem_lh_r', 'semiten_r', 'tib_ant_r', 'gas_med_r', 'soleus_r', 'glut_med2_r'}; 

% Use this structure if you want to use one EMG channel for multiple
% muscles in the opensim model. The first name of each row is the reference 
% and should always be in the header of the EMGfile or in the  EMGheaders.
Misc.EMG_MuscleCopies = {'gas_med_r','gas_lat_r'};       %  use gastrocnemius medialis EMG to constrain activity of the lateral gastrocn




% Use this weird implementation with adapting to bounds based on activity
% of a static optimization without constraints
Misc.BoundsScale_EMG = 1;
% or the implementation when the bounds depend on the current simulated
% activity: Note currently not working

% Idea: a(t) bmin < a(t) - s EMG(t) < a(t) bmax
%   => we don't want to divide by zero. Hence implementation:
%           (1) 0 < a(t) (1-bmin) - sEMG(t) 
%           (2) a(t) (1-bmax) - sEMG(t) <0
%   => implemented this in constraints and jacobian, but apparantly not
%   solution possible. Ask Friedl about this.
Misc.ActBound = 1;


% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;


% Run muscle tendon estimator:
% ...... (Run function here)



