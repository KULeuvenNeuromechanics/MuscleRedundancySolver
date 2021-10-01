
clear all; close all
clc;
%% Input information
% Add here the paths of IK, ID , US and EMG data trials you want to work with
Misc.IKfile  = {fullfile(pwd,'Anterior_1_1.mot')};
Misc.IDfile  = {fullfile(pwd,'inverse_dynamics.sto')};
model_path   = {fullfile(pwd,'ModelNoArms_scaled.osim')};
Out_path     = fullfile(pwd,'Results_MRS');                    % folder to store results
time = [1.5 5]; 

% Not tracking EMG/US
Misc.UStracking = 0; 
Misc.EMGconstr = 0;

%% Settings

% name of the results file
Misc.OutName ='PerturbedBalance_Anterior_1_1';

% select degrees of freedom
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r',...
    'hip_adduction_r','hip_rotation_r','ankle_angle_l','knee_angle_l',...
    'hip_flexion_l','hip_adduction_l','hip_rotation_l'};    % select the DOFs you want to include in the optimization

% Set the tendon stifness of all muscles
Misc.kT = [];      % default way to set tendon stiffenss (default values is 35)

% select muscles
Misc.MuscleNames_Input = []; % select muscles

% Plotter Bool: Boolean to select if you want to plot lots of output information of intermediate steps in the script
Misc.PlotBool = 1;
% MRS Bool: Select if you want to run the generic muscle redundancy solver
Misc.MRSBool = 1;
% Validation Bool: Select if you want to run the muscle redundancy solver with the optimized parameters
Misc.ValidationBool = 0; 	% TO DO: we should report results of EMG driven simulation as well

%% Run muscle tendon estimator:
[Results,DatStore,Misc] = solveMuscleRedundancy(model_path,time,Out_path,Misc);
