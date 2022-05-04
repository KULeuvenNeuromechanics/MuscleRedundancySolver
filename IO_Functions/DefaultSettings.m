function [Misc] = DefaultSettings(Misc)
% --------------------------------------------------------------------------
%DefaultSettings
%     If user does not specify values for specific Misc fields, these are
%     filled out with default values. This function also corrects the
%     format of certain inputs, as required by the code.
% 
% INPUT:
%     Misc
%     Miscellaneous info used through the code
% 
% OUTPUT:
%     Misc
%     Miscellaneous info used through the code
% 
% Original author: 
% Original date: 
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

%% Filters
if ~isfield(Misc,'f_cutoff_ID') || isempty(Misc.f_cutoff_ID)
    Misc.f_cutoff_ID=6;
end
if ~isfield(Misc,'f_order_ID') || isempty(Misc.f_order_ID)
    Misc.f_order_ID=6;
end
% Muscle-tendon lengths
if ~isfield(Misc,'f_cutoff_lMT') || isempty(Misc.f_cutoff_lMT)
    Misc.f_cutoff_lMT=6;
end
if ~isfield(Misc,'f_order_lMT') || isempty(Misc.f_order_lMT)
    Misc.f_order_lMT=6;
end
% Moment arms
if ~isfield(Misc,'f_cutoff_dM') || isempty(Misc.f_cutoff_dM)
    Misc.f_cutoff_dM=6;
end
if ~isfield(Misc,'f_order_dM') || isempty(Misc.f_order_dM)
    Misc.f_order_dM=6;
end
% Inverse Kinematics
if ~isfield(Misc,'f_cutoff_IK') || isempty(Misc.f_cutoff_IK)
    Misc.f_cutoff_IK=6;
end
if ~isfield(Misc,'f_order_IK') || isempty(Misc.f_order_IK)
    Misc.f_order_IK=6;
end

%% Mesh Frequency
if ~isfield(Misc,'Mesh_Frequency') || isempty(Misc.Mesh_Frequency)
    Misc.Mesh_Frequency=100;
end

%% Run muscle analysis (default true)
if ~isfield(Misc,'RunAnalysis') || isempty(Misc.RunAnalysis)
    Misc.RunAnalysis = 1;
end

%% EMG
if isfield(Misc,'EMGheaders')
    if size(Misc.EMGheaders,1) == 1
        tempHeaders = Misc.EMGheaders;
        Misc.EMGheaders = [];
        for i=1:length(Misc.IKfile)
            Misc.EMGheaders{i} = tempHeaders;
        end
    end
end

if isfield(Misc,'EMGheaders') && isfield(Misc,'EMGFileHeaderCorrespondence')
    error('Only one of Misc.EMGheaders or Misc.EMGFileHeaderCorrespondence is allowed. Please remove either Misc.EMGheaders or Misc.EMGFileHeaderCorrespondence')
end
    
if ~isfield(Misc,'EMGconstr') 
    Misc.EMGconstr  = 0;
    boolEMG = 0;
elseif  isempty(Misc.EMGconstr)
    Misc.EMGconstr  = 0;
    boolEMG = 0;
elseif Misc.EMGconstr == 0
    boolEMG = 0;
elseif Misc.EMGconstr == 1
    boolEMG = 1;
end
Misc.boolEMG = boolEMG;

% bounds on scaling EMG
if ~isfield(Misc,'BoundsScaleEMG') || isempty(Misc.BoundsScaleEMG)
    Misc.BoundsScaleEMG = [0.9 1.1];
end
% bounds on EMG signal
if ~isfield(Misc,'EMGbounds')
    Misc.EMGbounds = [];
end
% copies of EMG
if isfield(Misc,'EMG_MuscleCopies')
    if size(Misc.EMG_MuscleCopies,2) > 2
        tempEMG_MuscleCopies = Misc.EMG_MuscleCopies;
        Misc.EMG_MuscleCopies = [];
        ct = 0;
        for r = 1:size(tempEMG_MuscleCopies,1)
            for c = 1:size(tempEMG_MuscleCopies,2)-1
                ct = ct + 1;
                Misc.EMG_MuscleCopies{ct,1} = tempEMG_MuscleCopies{r,c};
                Misc.EMG_MuscleCopies{ct,2} = tempEMG_MuscleCopies{r,c+1};
            end
        end
    end
else
    Misc.EMG_MuscleCopies =[];
end
if ~isfield(Misc,'EMG_MuscleCopies')
    Misc.EMG_MuscleCopies = [];
end
% selected EMG information
if ~isfield(Misc,'EMGSelection')
    Misc.EMGSelection = [];
end

%% Info related to parameter optimizatiEMGbounds on
if isfield(Misc,'Estimate_OptimalFiberLength')
    [r c] = size(Misc.Estimate_OptimalFiberLength);
    if ((r == 1) && (c > 1))
        Misc.Estimate_OptimalFiberLength = Misc.Estimate_OptimalFiberLength';
    end
else
    Misc.Estimate_OptimalFiberLength =[];
end
if isfield(Misc,'Estimate_TendonSlackLength')
    [r c] = size(Misc.Estimate_TendonSlackLength);
    if ((r == 1) && (c > 1))
        Misc.Estimate_TendonSlackLength = Misc.Estimate_TendonSlackLength';
    end
else
    Misc.Estimate_TendonSlackLength =[];
end
if isfield(Misc,'Estimate_TendonStiffness')
    [r c] = size(Misc.Estimate_TendonStiffness);
    if ((r == 1) && (c > 1))
        Misc.Estimate_TendonStiffness = Misc.Estimate_TendonStiffness';
    end
else
    Misc.Estimate_TendonStiffness =[];
end

if isfield(Misc,'Coupled_fiber_length')
    if size(Misc.Coupled_fiber_length,2) > 2
        tempCoupled_fiber_length = Misc.Coupled_fiber_length;
        Misc.Coupled_fiber_length = [];
        ct = 0;
        for r = 1:size(tempCoupled_fiber_length,1)
            for c = 1:size(tempCoupled_fiber_length,2)-1
                ct = ct + 1;
                Misc.Coupled_fiber_length{ct,1} = tempCoupled_fiber_length{r,c};
                Misc.Coupled_fiber_length{ct,2} = tempCoupled_fiber_length{r,c+1};
            end
        end
    end
else
    Misc.Coupled_fiber_length =[];
end

if isfield(Misc,'Coupled_slack_length')
    if size(Misc.Coupled_slack_length,2) > 2
        tempCoupled_slack_length = Misc.Coupled_slack_length;
        Misc.Coupled_slack_length = [];
        ct = 0;
        for r = 1:size(tempCoupled_slack_length,1)
            for c = 1:size(tempCoupled_slack_length,2)-1
                ct = ct + 1;
                Misc.Coupled_slack_length{ct,1} = tempCoupled_slack_length{r,c};
                Misc.Coupled_slack_length{ct,2} = tempCoupled_slack_length{r,c+1};
            end
        end
    end
else
    Misc.Coupled_slack_length =[];
end

if isfield(Misc,'Coupled_TendonStiffness')
    if size(Misc.Coupled_TendonStiffness,2) > 2
        tempCoupled_TendonStiffness = Misc.Coupled_TendonStiffness;
        Misc.Coupled_TendonStiffness = [];
        ct = 0;
        for r = 1:size(tempCoupled_TendonStiffness,1)
            for c = 1:size(tempCoupled_TendonStiffness,2)-1
                ct = ct + 1;
                Misc.Coupled_TendonStiffness{ct,1} = tempCoupled_TendonStiffness{r,c};
                Misc.Coupled_TendonStiffness{ct,2} = tempCoupled_TendonStiffness{r,c+1};
            end
        end
    end
else
    Misc.Coupled_TendonStiffness =[];
end

if ~isfield(Misc,'lb_kT_scaling')
    Misc.lb_kT_scaling = 0.2;
end

if ~isfield(Misc,'ub_kT_scaling')
    Misc.ub_kT_scaling = 1.2;
end

if ~isfield(Misc,'lb_lMo_scaling')
    Misc.lb_lMo_scaling = 0.7; % Lower bound for scaling optimal fiber length
end
if ~isfield(Misc,'ub_lMo_scaling')
    Misc.ub_lMo_scaling = 1.5; % Upper bound for scaling optimal fiber length
end
if ~isfield(Misc,'lb_lTs_scaling')
    Misc.lb_lTs_scaling = 0.7; % Lower bound for scaling tendon slack length
end
if ~isfield(Misc,'ub_lTs_scaling')
    Misc.ub_lTs_scaling = 1.5; % Upper bound for scaling tendon slack length
end

%% ResultsName
if ~isfield(Misc,'OutName') || isempty(Misc.OutName)
    Misc.OutName = '';
end
% filename of the osim model with updated parameters
if ~isfield(Misc,'newModelFile')
    file_path = char(Misc.model_path);
    [~,oldModelFile,~] = fileparts(file_path);
    Misc.newModelFile = [oldModelFile '_newParams.osim']; 
end

%% weights
if ~isfield(Misc,'wlM') || isempty(Misc.wlM)
    Misc.wlM = 0;
end
if ~isfield(Misc,'wEMG')|| isempty(Misc.wEMG)
    Misc.wEMG = 0;
end
if ~isfield(Misc,'wAct')|| isempty(Misc.wAct)
    Misc.wAct = 1;
end
if ~isfield(Misc,'wTres')|| isempty(Misc.wTres)
    Misc.wTres = 1000;
end
if ~isfield(Misc,'wVm')|| isempty(Misc.wVm)
    Misc.wVm = 0.01;
end

%% reserve actuator
if ~isfield(Misc,'Topt')|| isempty(Misc.Topt)
    Misc.Topt = 150;
end
if ~isfield(Misc,'USSelection')
    Misc.USSelection = {};
end

%% ultrasound
if ~isfield(Misc,'UStracking') || isempty(Misc.UStracking)
    Misc.UStracking = 0;
end
if ~isfield(Misc,'USfile')
    Misc.USfile = [];
end

%% Plotter
if ~isfield(Misc,'PlotBool')
    Misc.PlotBool = 0;
end

%% Validation
if ~isfield(Misc,'ValidationBool')
    Misc.ValidationBool = 1;
end

%% Solve muscle redundancy
if ~isfield(Misc,'MRSbool')
    Misc.MRSbool = 1;
end

%% Muscle parameters
if ~isfield(Misc,'MuscleNames_Input')
    Misc.MuscleNames_Input = cell(length(Misc.IKfile),1);
else
    [ntr_muscles, nmuscles] = size(Misc.MuscleNames_Input);
    if ntr_muscles == 1 || ntr_muscles~=length(Misc.IKfile)
        % select the same input dofs for all trials
        Misc.MuscleNames_Input_copy = Misc.MuscleNames_Input;
        Misc.MuscleNames_Input =cell(0);
        for t = 1:length(Misc.IKfile)
            Misc.MuscleNames_Input{t} = Misc.MuscleNames_Input_copy;
        end
    end
end

%% Activation and contraction dynamics
% ----------------------------------------------------------------------- %
if ~isfield(Misc,'tau_act')
    Misc.tau_act = 0.015;
end
if ~isfield(Misc,'tau_deact')
    Misc.tau_deact = 0.06;
end
if ~isfield(Misc,'b') % tanh coefficient for smooth activation dynamics
    Misc.b = 0.1;
end

%% Input Dofs

[ntr_dof, ndofs] = size(Misc.DofNames_Input);
if ntr_dof == 1 || ntr_dof~=length(Misc.IKfile)
    % select the same input dofs for all trials
    DofNames_Input_copy = Misc.DofNames_Input;
    Misc.DofNames_Input =cell(0);
    for t = 1:length(Misc.IKfile)
        Misc.DofNames_Input{t,1} = DofNames_Input_copy;
    end
end

%% update size EMG headers if needed
if ~isfield(Misc,'normalizeToMRS')
    Misc.normalizeToMRS = false;
end

