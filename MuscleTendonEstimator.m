function [Results,Parameters,DatStore] = MuscleTendonEstimator(model_path,IK_path,ID_path,US_path,time,Bounds,OutPath,Misc)
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
close all;
% Default low-pass filter:
%   Butterworth order: 6
%   Cutoff frequency: 6Hz
% Inverse Dynamics
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
% Mes Frequency
if ~isfield(Misc,'Mesh_Frequency') || isempty(Misc.Mesh_Frequency)
    Misc.Mesh_Frequency=100;
end

% ----------------------------------------------------------------------- %
% Muscle analysis ------------------------------------------------------- %
% Perform muscle analysis for the different selected trials
Misc.nTrials = size(IK_path,1);

DatStore = struct;
for i = 1:size(IK_path,1)
    IK_path_trial = IK_path(i,:);
    ID_path_trial = ID_path(i,:);
    
    Misc.time=time;
    MuscleAnalysisPath=fullfile(OutPath,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end
    disp('MuscleAnalysis Running .....');
    OpenSim_Muscle_Analysis(IK_path_trial,model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input)
    disp('MuscleAnalysis Finished');
    Misc.MuscleAnalysisPath=MuscleAnalysisPath;
    
    % ----------------------------------------------------------------------- %
    % Extract muscle information -------------------------------------------- %
    % Get number of degrees of freedom (dofs), muscle-tendon lengths and moment
    % arms for the selected muscles.
    [~,Misc.trialName,~]=fileparts(IK_path_trial);
    if ~isfield(Misc,'MuscleNames_Input') || isempty(Misc.MuscleNames_Input)
        Misc=getMuscles4DOFS(Misc);
    end
    % Shift tendon force-length curve as a function of the tendon stiffness
    Misc.shift = getShift(Misc.Atendon);
    [DatStore] = getMuscleInfo(IK_path_trial,ID_path_trial,Misc,DatStore,i);
    
    DatStore(i).free_lMo = zeros(length(Misc.Estimate_OptFL),1);
    for j = 1:length(DatStore(i).free_lMo)
        DatStore(i).free_lMo(j) = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_OptFL{j}));
    end
    
end


% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization --------- %

% NOTE: We do not estimate any parameters here, but these results can serve as
% decent initial guess for the later dynamic optimization

% The solution of the static optimization is used as initial guess for the
% dynamic optimization
% Extract the muscle-tendon properties
[Misc.params,Misc.lOpt,Misc.L_TendonSlack,Misc.Fiso,Misc.PennationAngle]=ReadMuscleParameters(model_path,DatStore(1).MuscleNames);



% Static optimization using IPOPT solver
for trial = 1:size(IK_path,1)
    DatStore = SolveStaticOptimization_IPOPT_CasADi(DatStore,Misc,trial);
end
% Show difference in fiber lengths from scaled model with rigid tendon vs
% experimental:





% update Tendon stiffness for specific muscles based on input arguments (we
% should replace this function with bounds (we can use it as inspiration
%if isfield(Misc,'Set_ATendon_ByName') && ~isempty(Misc.Set_ATendon_ByName)
%   [Misc,DatStore] = set_ATendon_ByName(Misc,DatStore);
%end

% get the EMG information
% [DatStore] = GetEMGInfo(Misc,DatStore);


% ....... To be continued......


%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART II: OPTIMAL CONTROL PROBLEM FORMULATION -------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% Input arguments
Misc.params = Misc.params;       % Muscle-tendon parameters
Misc.scaling.vMtilde = 10;           % Scaling factor: derivative muscle fiber lengths
Misc.w1 = 1000;                      % Weight objective function
Misc.w2 = 0.01;
Misc.Topt = 150;                     % Scaling factor: reserve actuators

tau_act = 0.015; Misc.tauAct = tau_act * ones(1,DatStore(1).nMuscles);       % activation time constant (activation dynamics)
tau_deact = 0.06; Misc.tauDeact = tau_deact * ones(1,DatStore(1).nMuscles);  % deactivation time constant (activation dynamics)
Misc.b = 0.1;

% Parameters of active muscle force-velocity characteristic
load('ActiveFVParameters.mat','ActiveFVParameters');
Fvparam(1) = 1.475*ActiveFVParameters(1); Fvparam(2) = 0.25*ActiveFVParameters(2);
Fvparam(3) = ActiveFVParameters(3) + 0.75; Fvparam(4) = ActiveFVParameters(4) - 0.027;
Misc.Fvparam = Fvparam;

% Parameters of active muscle force-length characteristic
load('Faparam.mat','Faparam');
Misc.Faparam = Faparam;

% Parameters of passive muscle force-length characteristic
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
Misc.Fpparam = [pp1;pp2];
Misc.Atendon=Misc.Atendon;
Misc.shift=Misc.shift;

% Problem bounds
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;      % bounds on normalized muscle fiber velocity
lMtilde_min = 0.2; lMtilde_max = 1.8;   % bounds on normalized muscle fiber length

for trial = 1:size(IK_path,1)
    DatStore(trial).NMuscles = DatStore(trial).nMuscles;   % number of muscles
    DatStore(trial).Ndof = DatStore(trial).nDOF;           % number of dofs
    % ADiGator works with 2D: convert 3D arrays to 2D structure (moment arms)
    for i = 1:DatStore(trial).Ndof
        DatStore(trial).MA(i).Joint(:,:) = DatStore(trial).dM(:,i,:);  % moment arms
    end
    % Time bounds
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    
    % Discretization
    N = round((tf-t0)*Misc.Mesh_Frequency);
    h = (tf-t0)/N;
    
    % Spline approximation of muscle-tendon length (LMT), moment arms (MA) and inverse dynamic torques (ID)
    for dof = 1:DatStore(trial).Ndof
        for m = 1:DatStore(trial).NMuscles
            DatStore(trial).JointMASpline(dof).Muscle(m) = spline(DatStore(trial).time,DatStore(trial).MA(dof).Joint(:,m));
        end
        DatStore(trial).JointIDSpline(dof) = spline(DatStore(trial).time,DatStore(trial).T_exp(:,dof));
    end
    
    for m = 1:DatStore(trial).NMuscles
        DatStore(trial).LMTSpline(m) = spline(DatStore(trial).time,DatStore(trial).LMT(:,m));
    end
    
    % Evaluate LMT, VMT, MA and ID at optimization mesh
    step = (tf-t0)/(N);
    time_opt = t0:step:tf;
    DatStore(trial).LMTinterp = zeros(length(time_opt),DatStore(trial).NMuscles); % Muscle-tendon length
    for m = 1:DatStore(trial).NMuscles
        [DatStore(trial).LMTinterp(:,m),~,~] = SplineEval_ppuval(DatStore(trial).LMTSpline(m),time_opt,1);
    end
    DatStore(trial).MAinterp = zeros(length(time_opt),DatStore(trial).Ndof*DatStore(trial).NMuscles); % Moment arm
    DatStore(trial).IDinterp = zeros(length(time_opt),DatStore(trial).Ndof); % Inverse dynamic torque
    for dof = 1:DatStore(trial).Ndof
        for m = 1:DatStore(trial).NMuscles
            index_sel=(dof-1)*(DatStore(trial).NMuscles)+m;
            DatStore(trial).MAinterp(:,index_sel) = ppval(DatStore(trial).JointMASpline(dof).Muscle(m),time_opt);
        end
        DatStore(trial).IDinterp(:,dof) = ppval(DatStore(trial).JointIDSpline(dof),time_opt);
    end
    
    % Initial guess static optimization
    DatStore(trial).SoActInterp = interp1(DatStore(trial).time,DatStore(trial).SoAct,time_opt);
    DatStore(trial).SoRActInterp = interp1(DatStore(trial).time,DatStore(trial).SoRAct,time_opt);
    DatStore(trial).SoForceInterp = interp1(DatStore(trial).time,DatStore(trial).SoForce.*DatStore(trial).cos_alpha./Misc.Fiso,time_opt);
    [~,DatStore(trial).lMtildeInterp ] = FiberLength_Ftilde(DatStore(trial).SoForceInterp,Misc.params,DatStore(trial).LMTinterp,Misc.Atendon,Misc.shift);
    DatStore(trial).vMtildeinterp = zeros(size(DatStore(trial).lMtildeInterp));
    for m = 1:DatStore(trial).NMuscles
        DatStore(trial).lMtildeSpline = spline(time_opt,DatStore(trial).lMtildeInterp(:,m));
        [~,DatStore(trial).vMtildeinterp_norm,~] = SplineEval_ppuval(DatStore(trial).lMtildeSpline,time_opt,1);
        DatStore(trial).vMtildeinterp(:,m) = DatStore(trial).vMtildeinterp_norm/Misc.scaling.vMtilde;
    end
end

% CasADi setup
import casadi.*
% Load CasADi function
% CasADiFunctions

opti = casadi.Opti();

% What if different trials involve different nr of muscles?
if Misc.MRSBool == 1
    nTrials = Misc.nTrials;
    N_tot = 0;
    SoActGuess = []; lMtildeGuess = []; vMtildeGuess = []; SoRActGuess = []; SoExcGuess = [];
    for trial = 1:nTrials
        % Time bounds
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        % Discretization
        N = round((tf-t0)*Misc.Mesh_Frequency);
        N_tot = N_tot + N;
        SoActGuess = [SoActGuess DatStore(trial).SoActInterp'];
        SoExcGuess = [SoExcGuess DatStore(trial).SoActInterp(1:end-1,:)'];
        lMtildeGuess = [lMtildeGuess DatStore(trial).lMtildeInterp'];
        vMtildeGuess = [vMtildeGuess DatStore(trial).vMtildeinterp(1:end-1,:)'];
        SoRActGuess = [SoRActGuess DatStore(trial).SoRActInterp(1:end-1,:)'];
        
    end
    
    % Variables - bounds and initial guess
    % States (at mesh and collocation points)
    % Muscle activations
    a = opti.variable(DatStore(trial).NMuscles,N_tot+nTrials);      % Variable at mesh points
    opti.subject_to(a_min < a < a_max);           % Bounds
    opti.set_initial(a,SoActGuess);             % Initial guess (static optimization)
    % Muscle fiber lengths
    lMtilde = opti.variable(DatStore(trial).NMuscles,N_tot+nTrials);
    opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
    opti.set_initial(lMtilde,lMtildeGuess);
    
    % Controls
    e = opti.variable(DatStore(trial).NMuscles,N_tot);
    opti.subject_to(e_min < e < e_max);
    opti.set_initial(e, SoExcGuess);
    % Reserve actuators
    aT = opti.variable(DatStore(trial).Ndof,N_tot);
    opti.subject_to(-1 < aT <1);
    opti.set_initial(aT,SoRActGuess./Misc.Topt);
    % Time derivative of muscle-tendon forces (states)
    vMtilde = Misc.scaling.vMtilde.*opti.variable(DatStore(trial).NMuscles,N_tot);
    opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
    opti.set_initial(vMtilde,vMtildeGuess);
    
    
    % Loop over mesh points formulating NLP
    J = 0; % Initialize cost function
    N_acc = 0;
    for trial = 1:Misc.nTrials
        % Time bounds
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        % Discretization
        N = round((tf-t0)*Misc.Mesh_Frequency);
        h = (tf-t0)/N;
        
        for k=1:N
            % Variables within current mesh interval
            ak = a(:,(N_acc+trial)*(trial-1) + k); lMtildek = lMtilde(:,(N_acc+trial)*(trial-1) + k);
            vMtildek = vMtilde(:,N_acc*(trial-1) + k); aTk = aT(:,N_acc*(trial-1) + k); ek = e(:,N_acc*(trial-1) + k);
            
            % Integration   Uk = (X_(k+1) - X_k)/*dt
            Xk = [ak; lMtildek];
            Zk = [a(:,k+1);lMtilde(:,k+1)];
            Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
            opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
            
            % Get muscle-tendon forces and derive Hill-equilibrium
            [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,DatStore(trial).LMTinterp(k,:)',Misc.params',Misc.Fvparam,Misc.Fpparam,Misc.Faparam,Misc.Atendon',Misc.shift');
            
            % Add path constraints
            % Moment constraints
            for dof = 1:DatStore(trial).Ndof
                T_exp = DatStore(trial).IDinterp(k,dof);
                index_sel = (dof-1)*(DatStore(trial).NMuscles)+1:(dof-1)*(DatStore(trial).NMuscles)+DatStore(trial).NMuscles;
                T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
                opti.subject_to(T_exp - T_sim == 0);
            end
            % Hill-equilibrium constraint
            opti.subject_to(Hilldiffk == 0);
            
        end
        
        J = J + ...
            sumsqr(e)/N/DatStore(trial).NMuscles + ...
            Misc.w1*sumsqr(aT)/N/DatStore(trial).Ndof + ...
            Misc.w2*sumsqr(vMtilde)/N/DatStore(trial).NMuscles;
        N_acc = N_acc + N;
    end
    opti.minimize(J); % Define cost function in opti
    
    % Create an NLP solver
    output.setup.Misc = Misc;
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
    
    opti.solver(output.setup.nlp.solver,optionssol);
    
    % Solve
    diary('GenericMRS.txt');
    sol = opti.solve();
    diary off
    
    %% Extract results
    % Variables at mesh points
    % Muscle activations and muscle-tendon forces
    a_opt = sol.value(a);
    lMtilde_opt = sol.value(lMtilde);
    % Muscle excitations
    e_opt = sol.value(e);
    % Reserve actuators
    aT_opt = sol.value(aT);
    % Time derivatives of muscle-tendon forces
    vMtilde_opt = sol.value(vMtilde);
    
    % Grid
    % Mesh points
    tgrid = linspace(t0,tf,N+1)';
    % Save results
    Time(trial).genericMRS = tgrid;
    MActivation(trial).genericMRS = a_opt;
    lMtildeopt(trial).genericMRS = lMtilde_opt;
    lM(trial).genericMRS = lMtilde_opt.*repmat(DatStore.lOpt',1,length(Time.genericMRS));
    MvMtilde(trial).genericMRS = vMtilde_opt;
    MExcitation(trial).genericMRS = e_opt;
    RActivation(trial).genericMRS = aT_opt*Misc.Topt;
    MuscleNames = DatStore.MuscleNames;
    OptInfo = output;
    % Tendon forces from lMtilde
    lMTinterp(trial).genericMRS = LMTinterp;
    [TForcetilde_,TForce_] = TendonForce_lMtilde(lMtildeopt.genericMRS',Misc.params,lMTinterp.genericMRS,Misc.Atendon,Misc.shift);
    TForcetilde(trial).genericMRS = TForcetilde_';
    TForce(trial).genericMRS = TForce_';
    
    clear opti a lMtilde e vMtilde aT
    
end

opti_MTE = casadi.Opti();
US_data = importdata(US_path);

% Variables - bounds and initial guess
% States (at mesh and collocation points)
% Muscle activations
a = opti_MTE.variable(DatStore.NMuscles,N+1);      % Variable at mesh points
opti_MTE.subject_to(a_min < a < a_max);           % Bounds
% Muscle fiber lengths
lMtilde = opti_MTE.variable(DatStore.NMuscles,N+1);
opti_MTE.subject_to(lMtilde_min < lMtilde < lMtilde_max);

% Controls
e = opti_MTE.variable(DatStore.NMuscles,N);
opti_MTE.subject_to(e_min < e < e_max);
% Reserve actuators
aT = opti_MTE.variable(DatStore.Ndof,N);
opti_MTE.subject_to(-1 < aT <1);
% Time derivative of muscle-tendon forces (states)
vMtilde = Misc.scaling.vMtilde.*opti_MTE.variable(DatStore.NMuscles,N);
opti_MTE.subject_to(vMtilde_min < vMtilde < vMtilde_max);

% Free optimal fiber length
lMo_scaling_param = opti_MTE.variable(DatStore.NMuscles,1);
lb_lMo_scaling = ones(DatStore.NMuscles,1); ub_lMo_scaling = ones(DatStore.NMuscles,1);
lb_lMo_scaling(DatStore.free_lMo(:)) = 0.8*lb_lMo_scaling(DatStore.free_lMo(:));
ub_lMo_scaling(DatStore.free_lMo(:)) = 1.3*ub_lMo_scaling(DatStore.free_lMo(:));
opti_MTE.subject_to(lb_lMo_scaling < lMo_scaling_param < ub_lMo_scaling);

% Free slack length
lTs_scaling_param = opti_MTE.variable(DatStore.NMuscles,1);
lb_lTs_scaling = ones(DatStore.NMuscles,1); ub_lTs_scaling = ones(DatStore.NMuscles,1);
lb_lTs_scaling(DatStore.free_lMo(:)) = 0.8*lb_lTs_scaling(DatStore.free_lMo(:));
ub_lTs_scaling(DatStore.free_lMo(:)) = 1.3*ub_lTs_scaling(DatStore.free_lMo(:));
opti_MTE.subject_to(lb_lTs_scaling < lTs_scaling_param < ub_lTs_scaling);

% Free tendon stifness
kT_scaling_param = opti_MTE.variable(DatStore.NMuscles,1);
lb_kT_scaling_param = ones(DatStore.NMuscles,1); ub_kT_scaling_param = ones(DatStore.NMuscles,1);
lb_kT_scaling_param(DatStore.free_lMo(:)) = 0.2*lb_kT_scaling_param(DatStore.free_lMo(:));
ub_kT_scaling_param(DatStore.free_lMo(:)) = 1.2*ub_kT_scaling_param(DatStore.free_lMo(:));
opti_MTE.subject_to(lb_kT_scaling_param < kT_scaling_param < ub_kT_scaling_param);

% Set initial guess
if Misc.MRSBool == 1
    opti_MTE.set_initial(a,a_opt);             % Initial guess generic MRS
    opti_MTE.set_initial(lMtilde,lMtilde_opt);
    opti_MTE.set_initial(e,e_opt);
    opti_MTE.set_initial(vMtilde,vMtilde_opt);
    opti_MTE.set_initial(aT,aT_opt);
    opti_MTE.set_initial(lMo_scaling_param,1);
    opti_MTE.set_initial(lTs_scaling_param,1);
    opti_MTE.set_initial(kT_scaling_param,1);
else
    opti_MTE.set_initial(a,SoActInterp');             % Initial guess (static optimization)
    opti_MTE.set_initial(lMtilde,lMtildeopt.genericMRS);
    opti_MTE.set_initial(e, SoActInterp(1:N,:)');
    opti_MTE.set_initial(vMtilde,vMtildeinterp(1:N,:)');
    opti_MTE.set_initial(aT,SoRActInterp(1:N,:)'./Misc.Topt);
    opti_MTE.set_initial(lMo_scaling_param,1);
    opti_MTE.set_initial(lTs_scaling_param,1);
    opti_MTE.set_initial(kT_scaling_param,1);
end


% Loop over mesh points formulating NLP
J = 0; % Initialize cost function
for k=1:N
    % Variables within current mesh interval
    ak = a(:,k); lMtildek = lMtilde(:,k);
    vMtildek = vMtilde(:,k); aTk = aT(:,k); ek = e(:,k);
    
    % Integration   Uk = (X_(k+1) - X_k)/*dt
    Xk = [ak; lMtildek];
    Zk = [a(:,k+1);lMtilde(:,k+1)];
    Uk = [f_ActivationDynamics(ek,ak); vMtildek];
    opti_MTE.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
    
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffk,FTk] = f_ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree(ak,lMtildek,vMtildek,LMTinterp(k,:)',[lMo_scaling_param lTs_scaling_param kT_scaling_param]);
    
    % Add path constraints
    % Moment constraints
    for dof = 1:DatStore.Ndof
        T_exp = IDinterp(k,dof);
        index_sel = (dof-1)*(DatStore.NMuscles)+1:(dof-1)*(DatStore.NMuscles)+DatStore.NMuscles;
        T_sim = f_spNMuscles(MAinterp(k,index_sel),FTk) + Misc.Topt*aTk(dof);
        opti_MTE.subject_to(T_exp - T_sim == 0);
    end
    % Hill-equilibrium constraint
    opti_MTE.subject_to(Hilldiffk == 0);
    
end
lMo = lMo_scaling_param(DatStore.free_lMo(:)).*Misc.params(2,DatStore.free_lMo(:))';
lMtilde_tracking = US_data.data(:,2)'./lMo/1000;
lMtilde_simulated = lMtilde(DatStore.free_lMo(:),:);
J = J + ...
    10*sumsqr(lMtilde_simulated-lMtilde_tracking)/N + ...
    0.5*(sumsqr(e)/N/DatStore.NMuscles + sumsqr(a)/N/DatStore.NMuscles) + ...
    Misc.w1*sumsqr(aT)/N/DatStore.Ndof + ...
    Misc.w2*sumsqr(vMtilde)/N/DatStore.NMuscles;

opti_MTE.minimize(J); % Define cost function in opti

% Create an NLP solver
output.setup.auxdata = auxdata;
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

opti_MTE.solver(output.setup.nlp.solver,optionssol);

% Solve
diary('MTE.txt');
sol = opti_MTE.solve();
diary off

%% Extract results
% Variables at mesh points
% Muscle activations and muscle-tendon forces
a_opt = sol.value(a);
lMtilde_opt = sol.value(lMtilde);
lMo_scaling_param_opt = sol.value(lMo_scaling_param);
lTs_scaling_param_opt = sol.value(lTs_scaling_param);
kT_scaling_param_opt = sol.value(kT_scaling_param);

lMo_opt = lMo_scaling_param_opt.*Misc.params(2,:)';
% Muscle excitations
e_opt = sol.value(e);
% Reserve actuators
aT_opt = sol.value(aT);
% Time derivatives of muscle-tendon forces
% vMtilde_opt = sol.value(vMtilde)';

% Grid
% Mesh points
tgrid = linspace(t0,tf,N+1)';

% Plot US tracking
figure(5)
plot(tgrid,lMtilde_opt(DatStore.free_lMo(:),:).*lMo_opt(DatStore.free_lMo(:)),'LineWidth',2); hold on;
plot(tgrid,lMtildeopt.genericMRS(DatStore.free_lMo(:),:).*Misc.params(2,DatStore.free_lMo(:))','LineWidth',2); hold on;
plot(tgrid,US_data.data(:,2)'/1000,'LineWidth',2);
legend('MTE','MRS','USdata');

% Save results
Time.MTE = tgrid;
MActivation.MTE = a_opt;
lMtildeopt.MTE = lMtilde_opt;
lM.MTE = lMtilde_opt.*lMo_opt;
MExcitation.MTE = e_opt;
RActivation.MTE = aT_opt*Misc.Topt;

lMo_scaling_paramopt.MTE = lMo_scaling_param_opt;
lTs_scaling_paramopt.MTE = lTs_scaling_param_opt;
kT_scaling_paramopt.MTE = kT_scaling_param_opt;

MuscleNames = DatStore.MuscleNames;
OptInfo = output;
% Tendon forces from lMtilde
lMTinterp.MTE = LMTinterp;
[TForcetilde_,TForce_] = TendonForce_lMtilde(lMtildeopt.MTE',Misc.params,lMTinterp.MTE,Misc.Atendon,Misc.shift);
TForcetilde.MTE = TForcetilde';
TForce.MTE = TForce';

clear opti a lMtilde e vMtilde aT

Results.Time = Time;
Results.MActivation = MActivation;
Results.lMtildeopt = lMtildeopt;
Results.lM = lM;
Results.MExcitation = MExcitation;
Results.RActivation = RActivation;
Results.lMo_scaling_paramopt = lMo_scaling_paramopt;
Results.lMTinterp = lMTinterp;
Results.TForcetilde = TForcetilde;
Results.TForce = TForce;

Parameters = [];

end

