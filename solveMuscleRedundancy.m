function [Results,DatStore,Misc] = solveMuscleRedundancy(time,Misc)
% -----------------------------------------------------------------------%
% INPUTS:
%           time: time window
%           Misc: structure with general input data (see manual for more details)
%
% OUTPUTS:
%           Results:    structure with outputs (states, controls, ...)
%           DatStore:   structure with data used for solving the optimal control problem
%           Misc:       structure with general input data (see manual for more details)
% -----------------------------------------------------------------------%

% update default settings
if isfield(Misc,'EMGconstr') && Misc.EMGconstr == 1
    % initialize certain settings
    Misc.EMGheaders = cell(1,length(Misc.IKfile));
end
Misc = DefaultSettings(Misc);

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%check if we have to adapt the start and end time so that it corresponds to
%time frames in the IK solution
[time] = Check_TimeIndices(Misc,time);
Misc.time=time;

% ----------------------------------------------------------------------- %
%% Perform muscle analysis for all trials
DatStore = struct;
MuscleAnalysisPath=fullfile(Misc.OutPath,'MuscleAnalysis');
Misc.MuscleAnalysisPath=MuscleAnalysisPath;
if ~exist(MuscleAnalysisPath,'dir')
    mkdir(MuscleAnalysisPath);
    for i = 1:Misc.nTrials
        % select the IK and ID file
        IK_path_trial = Misc.IKfile{i};
        % Run muscle analysis    
        if Misc.RunAnalysis
            disp('MuscleAnalysis Running .....');
            OpenSim_Muscle_Analysis(IK_path_trial,Misc.model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input{i})
            disp('MuscleAnalysis Finished');
        end
    end
end

%% From here only use the trails selected
% ----------------------------------------------------------------------- %
% Extract muscle information -------------------------------------------- %
% Get number of degrees of freedom (dofs), muscle-tendon lengths, moment
% arms, stiffness and shift for the selected muscles.
for trial = Misc.trials_sel
    [~,Misc.MAtrialName{trial},~]=fileparts(Misc.IKfile{trial});
end

Misc = getMuscles4DOFS(Misc);

[Misc,DatStore] = getMuscleInfo(Misc,DatStore);
    
% display warnings in muscle selection
[Misc] = Warnings_MuscleNames(DatStore,Misc);

% get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
[DatStore,Misc] = GetIndices(DatStore,Misc);

% get the EMG information
[Misc,DatStore] = GetEMGInfo(Misc,DatStore);
% [DatStore] = GetEMGInfo_DG(Misc,DatStore);
[Misc,DatStore] = GetUSInfo(Misc,DatStore);

% get the number of muscles
for trial = Misc.trials_sel
    NMuscles(trial) = DatStore(trial).nMuscles;
end

%% Static optimization
% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization
% NOTE: We do not estimate any parameters here, but these results can serve as
% decent initial guess for the later dynamic optimization
% Extract the muscle-tendon properties
for i=Misc.trials_sel
    % Static optimization using IPOPT solver (used as an initial guess)
    DatStore = SolveStaticOptimization_IPOPT_CasADi(DatStore,Misc,i);
end

%% Input activation and contraction dynamics
% ----------------------------------------------------------------------- %
Misc.tau_act = 0.015;    
Misc.tau_deact = 0.06;   
Misc.b = 0.1;       % tanh coefficient for smooth activation dynamics

% Misc.kT=Misc.kT;
% Misc.shift=Misc.shift;

%% Descretisation

% mesh descretisation
for trial = Misc.trials_sel
    Misc.tauAct{trial} = Misc.tau_act * ones(NMuscles(trial), 1);       % activation time constant (activation dynamics)
    Misc.tauDeact{trial} = Misc.tau_deact * ones(NMuscles(trial),1);  % deactivation time constant (activation dynamics)
    
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    Mesh(trial).N = round((tf-t0)*Misc.Mesh_Frequency);
    Mesh(trial).step = (tf-t0)/Mesh(trial).N;
    Mesh(trial).t = t0:Mesh(trial).step:tf;
end

%% Evaluate splines at Mesh Points
% ----------------------------------------------------------------------- %
% Get IK, ID, muscle analysis and static opt information at mesh points

for trial = Misc.trials_sel
    % Discretization
    N = Mesh(trial).N;
    time_opt = Mesh(trial).t;    
    % Spline approximation of muscle-tendon length (LMT), moment arms (MA) and inverse dynamic torques (ID)
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles(trial)
            DatStore(trial).JointMASpline(dof).Muscle(m) = spline(DatStore(trial).time,squeeze(DatStore(trial).dM(:,dof,m)));
        end
        DatStore(trial).JointIDSpline(dof) = spline(DatStore(trial).time,DatStore(trial).T_exp(:,dof));
    end
    
    for m = 1:NMuscles(trial)
        DatStore(trial).LMTSpline(m) = spline(DatStore(trial).time,DatStore(trial).LMT(:,m));
    end
    
    % Evaluate LMT, VMT, MA and ID at optimization mesh
    DatStore(trial).LMTinterp = zeros(length(time_opt),NMuscles(trial)); % Muscle-tendon length
    for m = 1:NMuscles(trial)
        [DatStore(trial).LMTinterp(:,m),~,~] = SplineEval_ppuval(DatStore(trial).LMTSpline(m),time_opt,1);
    end
    DatStore(trial).MAinterp = zeros(length(time_opt),DatStore(trial).nDOF*NMuscles(trial)); % Moment arm
    DatStore(trial).IDinterp = zeros(length(time_opt),DatStore(trial).nDOF); % Inverse dynamic torque
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles(trial)
            index_sel=(dof-1)*NMuscles(trial)+m;
            DatStore(trial).MAinterp(:,index_sel) = ppval(DatStore(trial).JointMASpline(dof).Muscle(m),time_opt);
        end
        DatStore(trial).IDinterp(:,dof) = ppval(DatStore(trial).JointIDSpline(dof),time_opt);
    end
    
    % Interpolate results of static optimization
    DatStore(trial).SoActInterp = interp1(DatStore(trial).time,DatStore(trial).SoAct,time_opt');
    DatStore(trial).SoRActInterp = interp1(DatStore(trial).time,DatStore(trial).SoRAct,time_opt');
    DatStore(trial).SoForceInterp = interp1(DatStore(trial).time,DatStore(trial).SoForce.*DatStore(trial).cos_alpha./Misc.FMo(:,Misc.idx_allMuscleList{trial}),time_opt);
    [~,DatStore(trial).lMtildeInterp ] = FiberLength_Ftilde(DatStore(trial).SoForceInterp,Misc.params(:,Misc.idx_allMuscleList{trial}),DatStore(trial).LMTinterp,Misc.kT(:,Misc.idx_allMuscleList{trial}),Misc.shift(:,Misc.idx_allMuscleList{trial}));
    DatStore(trial).vMtildeinterp = zeros(size(DatStore(trial).lMtildeInterp));
    for m = 1:NMuscles(trial)
        DatStore(trial).lMtildeSpline = spline(time_opt,DatStore(trial).lMtildeInterp(:,m));
        [~,DatStore(trial).vMtildeinterp_norm,~] = SplineEval_ppuval(DatStore(trial).lMtildeSpline,time_opt,1);
        DatStore(trial).vMtildeinterp(:,m) = DatStore(trial).vMtildeinterp_norm;
    end
end

%% setup options for the solver
% Create an NLP solver
% output.setup.lM_projecteddata = lM_projecteddata;
output.setup.nlp.solver = 'ipopt';
output.setup.nlp.ipoptoptions.linear_solver = 'mumps';
% Set derivativelevel to 'first' for approximating the Hessian
output.setup.derivatives.derivativelevel = 'second';
output.setup.nlp.ipoptoptions.tolerance = 1e-6;
output.setup.nlp.ipoptoptions.maxiterations = 10000;
if strcmp(output.setup.derivatives.derivativelevel, 'first')
    optionssol.ipopt.hessian_approximation = 'limited-memory';
end
% By default, the barrier parameter update strategy is monotone.
% https://www.coin-or.org/Ipopt/documentation/node46.html#SECTION000116020000000000000
% Uncomment the following line to use an adaptive strategy
% optionssol.ipopt.mu_strategy = 'adaptive';
optionssol.ipopt.nlp_scaling_method = 'gradient-based';
optionssol.ipopt.linear_solver = output.setup.nlp.ipoptoptions.linear_solver;
optionssol.ipopt.tol = output.setup.nlp.ipoptoptions.tolerance;
optionssol.ipopt.max_iter = output.setup.nlp.ipoptoptions.maxiterations;


%% Dynamic Optimization - Default parameters
% ----------------------------------------------------------------------- %
% Solve muscle redundancy problem with default parameters
Results = struct;
if Misc.MRSBool == 1
    if Misc.MRS_validation_together
        [Results] = runMRS_allTrialsTogether(Misc,DatStore,Mesh,output,optionssol,Results,NMuscles);
    else
        for trial = Misc.trials_sel
            [Results] = runMRS(Misc,DatStore,Mesh,trial,output,optionssol,Results,NMuscles);
        end
    end
end

%% Dynamic Optimization - Parameter estimation
% ----------------------------------------------------------------------- %

% Parameter optimization selected if EMG information or ultrasound
% information is active
BoolParamOpt = 0;
if Misc.UStracking == 1 || Misc.EMGconstr == 1
    BoolParamOpt = 1;
end

if BoolParamOpt == 1
    [Results,Misc,DatStore,lMo_scaling_param_opt,lTs_scaling_param_opt,...
        kT_scaling_param_opt,EMGscale_opt] = runParameterEstimation(Misc,...
        DatStore,Mesh,output,optionssol,Results,NMuscles);
else
    lMo_scaling_param_opt = ones(NMuscles,1);
    lTs_scaling_param_opt = ones(NMuscles,1);
    kT_scaling_param_opt = ones(NMuscles,1);
end

%% Store Results

% save original and estimated parameters (and the bounds)
Results.Param.lMo_scaling_paramopt  = lMo_scaling_param_opt;
Results.Param.lTs_scaling_paramopt  = lTs_scaling_param_opt;
Results.Param.kT_scaling_paramopt   = kT_scaling_param_opt;
Results.Param.Original.FMo      = Misc.params(1,:);
Results.Param.Original.lMo      = Misc.params(2,:);
Results.Param.Original.lTs   = Misc.params(3,:);
Results.Param.Original.alphao = Misc.params(4,:);
Results.Param.Original.kT   = Misc.kT;
if BoolParamOpt
    Results.Param.Estimated.FMo    = Results.Param.Original.FMo;
    Results.Param.Estimated.lMo    = Results.Param.Original.lMo .* Results.Param.lMo_scaling_paramopt';
    Results.Param.Estimated.lTs    = Results.Param.Original.lTs .* Results.Param.lTs_scaling_paramopt';
    Results.Param.Estimated.alphao = Results.Param.Original.alphao;
    Results.Param.Estimated.kT     = Results.Param.Original.kT .* Results.Param.kT_scaling_paramopt';
    Results.Param.Bound.lMo.lb     = Misc.lb_lMo_scaling;
    Results.Param.Bound.lMo.ub     = Misc.ub_lMo_scaling;
    Results.Param.Bound.lTs.lb     = Misc.lb_lTs_scaling;
    Results.Param.Bound.lTs.ub     = Misc.ub_lTs_scaling;
    Results.Param.Bound.kT.lb      = Misc.lb_kT_scaling;
    Results.Param.Bound.kT.ub      = Misc.ub_kT_scaling;
    Results.Param.Bound.EMG.lb     = Misc.BoundsScaleEMG(1);
    Results.Param.Bound.EMG.ub     = Misc.BoundsScaleEMG(2);
    if Misc.boolEMG
        Results.Param.EMGscale     = EMGscale_opt;
    end
end

%% Validate results parameter estimation 

if Misc.ValidationBool == true && BoolParamOpt
    if Misc.MRS_validation_together
        [Results] = runValidation_allTrialsTogether(Misc,DatStore,Mesh,output,optionssol,Results,NMuscles,lMo_scaling_param_opt,lTs_scaling_param_opt,kT_scaling_param_opt);
    else
        for trial = Misc.trials_sel
            [Results] = runValidation(Misc,DatStore,Mesh,trial,output,optionssol,Results,NMuscles,lMo_scaling_param_opt,lTs_scaling_param_opt,kT_scaling_param_opt);
        end
    end
end
% add selected muscle names to the output structure
Results.MuscleNames = DatStore.MuscleNames;

%% Plot Output
% plot EMG tracking
if Misc.PlotBool && Misc.EMGconstr == 1
    h = PlotEMGTracking(Results,DatStore,Misc);
    if ~isdir(fullfile(Misc.OutPath,'figures'))
        mkdir(fullfile(Misc.OutPath,'figures'));
    end
    saveas(h,fullfile(Misc.OutPath,'figures',[Misc.analysisName '_fig_EMG.fig']));
end

% plot estimated parameters
if Misc.PlotBool == 1 && BoolParamOpt ==1
    h = PlotEstimatedParameters(Results,Misc);
    if ~isdir(fullfile(Misc.OutPath,'figures'))
        mkdir(fullfile(Misc.OutPath,'figures'));
    end
    saveas(h,fullfile(Misc.OutPath,'figures',[Misc.analysisName '_fig_Param.fig']));
end

% plot fiber length
if Misc.PlotBool && Misc.UStracking == 1
    h = PlotFiberLength(Results,DatStore);
    if ~isdir(fullfile(Misc.OutPath,'figures'))
        mkdir(fullfile(Misc.OutPath,'figures'));
    end
    saveas(h,fullfile(Misc.OutPath,'figures',[Misc.analysisName '_fig_FiberLength.fig']));
end

% plot the states of the muscles in the simulation
if Misc.PlotBool
    h = PlotStates(Results,DatStore,Misc);
    if ~isdir(fullfile(Misc.OutPath,'figures'))
        mkdir(fullfile(Misc.OutPath,'figures'));
    end
    saveas(h,fullfile(Misc.OutPath,'figures',[Misc.analysisName '_fig_States.fig']));
end

%% save the results
% plot states and variables from parameter estimation simulation
save(fullfile(Misc.OutPath,[Misc.analysisName 'Results.mat']),'Results','DatStore','Misc');

% write estimated parameters to new duplicate osim model
if BoolParamOpt
    muscleParams = Results.Param.Estimated;
    muscleNames  = Misc.allMuscleList;
    modelPath    = char(Misc.model_path);
    outPath      = Misc.OutPath;
    newModelFile = Misc.newModelFile;
    
    ParamsToOsim(muscleParams,muscleNames,modelPath,outPath,newModelFile); 
end

end


