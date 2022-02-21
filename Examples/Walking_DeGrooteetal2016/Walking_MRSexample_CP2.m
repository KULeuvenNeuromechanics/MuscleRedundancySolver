%% Default example solve muscle redundancy 
% (as in DeGroote2016)

clear all; clc; close all;
%% Input information

% select datafolder
MainPath = pwd;
DataPath = [pwd '/Data'];

% Add here the paths of IK, ID , US and EMG data trials you want to work with
% Misc.IKfile = {fullfile(DataPath,'Walking_IK.mot')};
% Misc.IDfile = {fullfile(DataPath,'Walking_ID.sto')};
% Misc.IKfile = strcat(extractfield(dir([DataPath '/IK/*.mot']),'folder'),...
%     {'\'},extractfield(dir([DataPath '/IK/*.mot']),'name'));
% Misc.IDfile = strcat(extractfield(dir([DataPath '/ID/*.sto']),'folder'),...
%     {'\'},extractfield(dir([DataPath '/ID/*.sto']),'name'));
Misc.subjectName = 'CP2';
model_path  = fullfile(DataPath,'Model','C2_CP2_T0_scaled_sf.osim');
OutPath    = fullfile(MainPath,'Results_3');                    % folder to store results
Misc.model_path = model_path;
Misc.OutPath = OutPath;

% Get start and end time of the different files
trialTypes = {'s2s','walking','squat','cmj'};
trialInfo = readtable(fullfile(DataPath,'Timepoints.xlsx'));
sides = {'L','R'};
trialInfo = trialInfo(ismember(trialInfo.Var3,trialTypes),:);
trialInfo = trialInfo(ismember(trialInfo.Var4,sides),:);
time = [trialInfo.Var5 trialInfo.Var6];
% time=[0.516 1.95]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)
Misc.IKfile = strcat({fullfile(DataPath,'IK\')},trialInfo.Var2,{'_IK.mot'})';
Misc.IDfile = strcat({fullfile(DataPath,'ID\')},trialInfo.Var2,{'_IK_ID.sto'})';
Misc.side = lower(trialInfo.Var4)';
Misc.trialType = trialInfo.Var3';
Misc.trialFP = trialInfo.Var7;

% select the DOFs you want to include in the optimization
% Add IPSA conditions
trialTypes_nonIPSA = {'s2s','walking','squat','cmj'};
DOF_nonIPSA = {'hip_flexion_','hip_adduction_','hip_rotation_','knee_angle_','ankle_angle_'};
for t=1:length(Misc.trialType)
    if ismember(Misc.trialType{t},trialTypes_nonIPSA)
        Misc.DofNames_Input{t}=strcat(DOF_nonIPSA,Misc.side{t});
    end
end

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;

% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;

% name output
Misc.OutName = strcat(trialInfo.Var3,{'_'})';

% adapt the stiffness of the achilles tendon (optional input argument)
% CHECK IF BOTH LEGS SHOULD BE CHANGED OR JUST THE SPASTIC LEG
Misc.Set_kT_ByName = {'soleus_r',20;'med_gas_r',20;'lat_gas_r',20;...
    'soleus_l',20;'med_gas_l',20;'lat_gas_l',20};

Misc = getMuscleProperties(model_path,Misc);
% NOTE: getShift was being used before adjusting kT, it has been corrected here

Misc.Mesh_Frequency=20;

%% Run muscle tendon estimator:
[Results,DatStore] = solveMuscleRedundancy_CP2(time,Misc);
