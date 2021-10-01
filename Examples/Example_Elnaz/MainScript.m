%% Example solve muscle redundancy with an MRI model

% clear variables and commond window
clear all; clc;

%% Input information
% Add here the paths of IK, ID , US and EMG data trials you want to work with
% As example we use trial 1 2 & 4 
Misc.IKfile  = {fullfile(pwd,'TRC_gait_6_IK.mot')};
Misc.IDfile  = {fullfile(pwd,'ID_6.sto')};
Misc.EMGfile = {fullfile(pwd,'EMG_gait_6.mot')};
model_path   = {fullfile(pwd,'gait2992_scaled.osim')};
Out_path     = fullfile(pwd,'Results2');                    % folder to store results
time         = [1.2 2.3]; 
% time = [1.7 1.9];

%% Settings

% name of the resuls file
Misc.OutName ='gait_6_EMGc';

% select degrees of freedom
Misc.DofNames_Input={'ankle_angle_l','knee_angle_l','hip_flexion_l','hip_adduction_l','hip_rotation_l'};    % select the DOFs you want to include in the optimization

% Provide the correct headers int case you EMG file has not the same
% headers as the muscle names in OpenSim (leave empty when you don't want
% to use this)
Misc.EMGheaders = {'Time','bifemlh_r','tib_ant_r','per_long_r','lat_gas_r','bifemsh_r','soleus_r','vas_lat_r','vas_med_r','per_brev_l','tib_ant_l','per_long_l',...
    'lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l','rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l','glut_med2_r','rect_fem_r'};

% channels you want to use for EMG constraints
% Misc.EMGSelection = {'per_brev_l','tib_ant_l','per_long_l','lat_gas_l','med_gas_l','soleus_l','vas_lat_l','vas_med_l','add_long_l','rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l'};

% Use this structure if you want to use one EMG channel for multiple
% muscles in the opensim model. The first name of each row is the reference
% and should always be in the header of the EMGfile or in the Misc.EMGheaders.
% and the selected EMG channels (Misc.EMGSelection).
Misc.EMG_MuscleCopies = {'glut_med2_l','glut_med1_l';
    'glut_med2_l','glut_med3_l'};       %  use gastrocnemius medialis EMG to constrain activity of the lateral gastrocn

% information for the EMG constraint
Misc.EMGconstr  = 1;                % Boolean to select EMG constrained option
Misc.BoundsScaleEMG = [0.1 1.2];    % maximal value to scale EMG

% bounds on EMG
Misc.EMGbounds    = [-0.05 0.05];   % upper and lower bound for deviation simulated and measured muscle activity (30% of EMG value)

% EMG part in the objective function
Misc.wEMG   = 0.1;   % weight on tracking EMG
Misc.wAct   = 0.001;
Misc.wTres  = 10;
Misc.wVm    = 0.00001;

% set the tendon stiffness
Misc.Set_kT_ByName =...
    {'soleus_l',20;
    'lat_gas_l',20;
    'med_gas_l',20};

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% no validation 
Misc.ValidationBool = 0;

%% Solve the muscle redundancy problem with EMG constraints

% here without soleus
Misc.EMGSelection = {'per_brev_l','tib_ant_l','per_long_l','lat_gas_l','med_gas_l','vas_lat_l','vas_med_l','add_long_l','rect_fem_l','tfl_l','glut_med2_l','bifemsh_l','bifemlh_l'};
[Results,DatStore,Misc] = solveMuscleRedundancy(model_path,time,Out_path,Misc);

% Compute the correlation between measured and simulated muscle activity

% step 1: get index of soleus (collom)
iSol_sim = find(strcmp(Results.Misc.MuscleNames_Input,'soleus_l'));

% step 2: get simulated activity soleus
% from Results
eSol_SimEMG = Results.MExcitation.genericMRS(iSol_sim,:)';    % simulated soleus activity with EMG constraints
eSol_Sim = Results.MExcitation.MTE(iSol_sim,:)';       % simulated without EMG constraints

% step 3: get measured activity
EMG = importdata(Misc.EMGfile{1});
EMG_headers = strsplit(EMG.textdata{end});
iSol = find(strcmp(EMG_headers,'Sol_l'));
eSol_EMG = EMG.data(:,iSol);

% interpolate measured EMG data to sampling rate of simulation
t_EMG = EMG.data(:,1);
tSim = Results.Time.genericMRS;
eSol_EMG_int = interp1(t_EMG,eSol_EMG,tSim(1:end-1));

% step 4:  bereken corrlatie opgementen EMG en simulatie
Corr_NoEMG = corr(eSol_EMG_int(1:end-1),eSol_Sim);
Corr_EMG = corr(eSol_EMG_int,eSol_SimEMG);

% check if this makes sense:
%%
figure();
plot(tSim(1:end-1), eSol_SimEMG,'b'); hold on;
plot(tSim(1:end-1), eSol_Sim,'r');
plot(tSim(1:end-1), eSol_EMG_int,'--k');
legend('EMG constr','Default','Measured');





