function [Results,Parameters,DatStore] = MuscleTendonEstimator(model_path,IK_path,ID_path,time,Bounds,OutPath,Misc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% -----------------------------------------------------------------------%
% INPUTS:
%           model_path: path to the .osim model
%           IK_path: path to the inverse kinematics results
%           ID_path: path to the inverse dynamics results
%           time: time window
%           Bounds: structure with bounds on states, controls, static param
%           OutPath: path to folder where results will be saved
%           Misc: structure of input data (see manual for more details)
%
% OUTPUTS:
%           Results:    structure with outputs (states, controls, ...)
%           Parameters: structure with static parameters
%           DatStore:   structure with data used for solving the optimal
%           control problem
% -----------------------------------------------------------------------%

% Update Misc Input to latest format
Misc = UpdateMiscInput(Misc);

% ----------------------------------------------------------------------- %
% Muscle analysis ------------------------------------------------------- %
Misc.time=time;
MuscleAnalysisPath=fullfile(OutPath,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end
disp('MuscleAnalysis Running .....');
OpenSim_Muscle_Analysis(IK_path,model_path,MuscleAnalysisPath,[time(1) time(end)],Misc.DofNames_Input)
disp('MuscleAnalysis Finished');
Misc.MuscleAnalysisPath=MuscleAnalysisPath;

% ----------------------------------------------------------------------- %
% Extract muscle information -------------------------------------------- %
% Get number of degrees of freedom (dofs), muscle-tendon lengths and moment
% arms for the selected muscles.
[~,Misc.trialName,~]=fileparts(IK_path);
if ~isfield(Misc,'MuscleNames_Input') || isempty(Misc.MuscleNames_Input)
    Misc=getMuscles4DOFS(Misc);
end
% Shift tendon force-length curve as a function of the tendon stiffness
Misc.shift = getShift(Misc.Atendon);
[DatStore] = getMuscleInfo(IK_path,ID_path,Misc);

% update Tendon stiffness for specific muscles based on input arguments (we
% should replace this function with bounds (we can use it as inspiration
%if isfield(Misc,'Set_ATendon_ByName') && ~isempty(Misc.Set_ATendon_ByName)
%   [Misc,DatStore] = set_ATendon_ByName(Misc,DatStore);
%end

% get the EMG information
[DatStore] = GetEMGInfo(Misc,DatStore);   

% Extract the muscle-tendon properties
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=...
    ReadMuscleParameters(model_path,DatStore.MuscleNames);

% ....... To be continued......


end

