function [Results,Parameters,DatStore] = MuscleTendonEstimator_vTom(model_path,time,Bounds,OutPath,Misc)
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

%% Default settings
% ----------------------------------------------------------------------- %
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
% Mesh Frequency
if ~isfield(Misc,'Mesh_Frequency') || isempty(Misc.Mesh_Frequency)
    Misc.Mesh_Frequency=100;
end
% Run muscle analysis (default true)
if ~isfield(Misc,'RunAnalysis') || isempty(Misc.RunAnalysis)
    Misc.RunAnalysis = 1;
end

%% Extract muscle information
% ----------------------------------------------------------------------- %
% Perform muscle analysis for the different selected trials
Misc.nTrials = size(Misc.IKfile,1);
DatStore = struct;
for i = 1:size(Misc.IKfile,1)
    IK_path_trial = Misc.IKfile{i};
    ID_path_trial = Misc.IDfile{i};
    
    Misc.time=time;
    MuscleAnalysisPath=fullfile(OutPath,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end    
    if Misc.RunAnalysis
        disp('MuscleAnalysis Running .....');
        OpenSim_Muscle_Analysis(IK_path_trial,model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input)
        disp('MuscleAnalysis Finished');        
    end
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
    
    % get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
    DatStore(i).free_lMo = zeros(length(Misc.Estimate_OptFL),1);
    DatStore(i).free_kT = zeros(length(Misc.Estimate_TendonStifness),1);
    DatStore(i).coupled_kT = zeros(size(Misc.Coupled_TendonStifness))';
    for j = 1:length(DatStore(i).free_lMo)
        DatStore(i).free_lMo(j) = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_OptFL{j}));
    end
    
    for j = 1:length(DatStore(i).free_kT)
        DatStore(i).free_kT(j) = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_TendonStifness{j}));
    end
    
    for k = 1:size((DatStore(i).coupled_kT),1)
        for j = 1:size((DatStore(i).coupled_kT),2)
            DatStore(i).coupled_kT(k,j) = find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_TendonStifness{j,k}));
        end
    end
end

% update Tendon stiffness for specific muscles based on input arguments (we
% should replace this function with bounds (we can use it as inspiration
if isfield(Misc,'Set_ATendon_ByName') && ~isempty(Misc.Set_ATendon_ByName)
   [Misc,DatStore] = set_ATendon_ByName(Misc,DatStore);
end

% get the EMG information
[DatStore] = GetEMGInfo(Misc,DatStore);
[DatStore] = GetUSInfo(Misc,DatStore);


%% Static optimization
% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization
% NOTE: We do not estimate any parameters here, but these results can serve as
% decent initial guess for the later dynamic optimization
% Extract the muscle-tendon properties
[Misc.params,Misc.lOpt,Misc.L_TendonSlack,Misc.Fiso,Misc.PennationAngle]=ReadMuscleParameters(model_path,DatStore(1).MuscleNames);

% Static optimization using IPOPT solver (used as an initial guess)
for trial = 1:Misc.nTrials
    DatStore = SolveStaticOptimization_IPOPT_CasADi(DatStore,Misc,trial);
end


%% Input activation and contraction dynamics
% ----------------------------------------------------------------------- %
Misc.params = Misc.params;          % Muscle-tendon parameters
Misc.scaling.vMtilde = 10;          % Scaling factor: derivative muscle fiber lengths
Misc.w1 = 1000;                     % Weight objective function
Misc.w2 = 0.01;
Misc.Topt = 150;                    % Scaling factor: reserve actuators

tau_act = 0.015;    Misc.tauAct = tau_act * ones(DatStore(1).nMuscles, 1);       % activation time constant (activation dynamics)
tau_deact = 0.06;   Misc.tauDeact = tau_deact * ones(DatStore(1).nMuscles,1);  % deactivation time constant (activation dynamics)
Misc.b = 0.1;       % tanh coefficient for smooth activation dynamics

Misc.Atendon=Misc.Atendon;
Misc.shift=Misc.shift;


%% Descretisation

% mesh descretisation
for trial = 1:Misc.nTrials
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    Mesh(trial).N = round((tf-t0)*Misc.Mesh_Frequency);    
    Mesh(trial).step = (tf-t0)/Mesh(trial).N;
    Mesh(trial).t = t0:Mesh(trial).step:tf;    
end


%% Evaluate splines at Mesh Points
% ----------------------------------------------------------------------- %
% Get IK, ID, muscle analysis and static opt information at mesh points
for trial = 1:size(Misc.IKfile,1)
    
    % Discretization
    N = Mesh(trial).N;
    time_opt = Mesh(trial).t;
    
    % Convert moment arms from 3D to 2D
    DatStore(trial).NMuscles = DatStore(trial).nMuscles;   % number of muscles
    DatStore(trial).Ndof = DatStore(trial).nDOF;           % number of dofs    
    for i = 1:DatStore(trial).Ndof
        DatStore(trial).MA(i).Joint(:,:) = DatStore(trial).dM(:,i,:);  % moment arms
    end
    
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

%% setup options for the solver
% Create an NLP solver
% output.setup.auxdata = auxdata;
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

% Problem bounds
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
lMtilde_min = 0.2; lMtilde_max = 1.8;   % bounds on normalized muscle fiber length

% CasADi setup
import casadi.*

% create opti structure
opti    = casadi.Opti();
nTrials = Misc.nTrials;
N_tot   = 0;
SoActGuess = []; lMtildeGuess = []; vMtildeGuess = []; SoRActGuess = []; SoExcGuess = [];
for trial = 1:nTrials
    % Update guess based on static optimization results
    N = Mesh(trial).N;
    N_tot = N_tot + N;
    SoActGuess = [SoActGuess DatStore(trial).SoActInterp'];
    SoExcGuess = [SoExcGuess DatStore(trial).SoActInterp(1:end-1,:)'];
    lMtildeGuess = [lMtildeGuess DatStore(trial).lMtildeInterp'];
    vMtildeGuess = [vMtildeGuess DatStore(trial).vMtildeinterp(1:end-1,:)'];
    SoRActGuess = [SoRActGuess DatStore(trial).SoRActInterp(1:end-1,:)'];
end

if Misc.MRSBool == 1
    % States 
    %   - muscle activation
    a = opti.variable(DatStore(trial).NMuscles,N_tot+nTrials);      % Variable at mesh points
    opti.subject_to(a_min < a < a_max);           % Bounds
    opti.set_initial(a,SoActGuess);             % Initial guess (static optimization)
    %   - Muscle fiber lengths
    lMtilde = opti.variable(DatStore(trial).NMuscles,N_tot+nTrials);
    opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
    opti.set_initial(lMtilde,lMtildeGuess);    
    %   - Controls
    e = opti.variable(DatStore(trial).NMuscles,N_tot);
    opti.subject_to(e_min < e < e_max);
    opti.set_initial(e, SoExcGuess);
    %   - Reserve actuators
    aT = opti.variable(DatStore(trial).Ndof,N_tot);
    opti.subject_to(-1 < aT <1);
    opti.set_initial(aT,SoRActGuess./Misc.Topt);
    %   - Time derivative of muscle-tendon forces (states)
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
        N = Mesh(trial).N;
        h = Mesh(trial).step;
        
        for k=1:N
            % Variables within current mesh interval
            ak = a(:,(N_acc+trial-1) + k); lMtildek = lMtilde(:,(N_acc+trial-1) + k);
            vMtildek = vMtilde(:,N_acc + k); aTk = aT(:,N_acc + k); ek = e(:,N_acc + k);
            
            % Integration   Uk = (X_(k+1) - X_k)/*dt
            Xk = [ak; lMtildek];
            Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
            Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
            opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
            
            % Get muscle-tendon forces and derive Hill-equilibrium
            [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,DatStore(trial).LMTinterp(k,:)',Misc.params',Misc.Atendon',Misc.shift');
            
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
            0.5*(sumsqr(e)/N/DatStore(trial).NMuscles + sumsqr(a)/N/DatStore(trial).NMuscles) + ...
            Misc.w1*sumsqr(aT)/N/DatStore(trial).Ndof + ...
            Misc.w2*sumsqr(vMtilde)/N/DatStore(trial).NMuscles;
        N_acc = N_acc + N;
    end
    opti.minimize(J); % Define cost function in opti
    
    % Create an NLP solver    
    opti.solver(output.setup.nlp.solver,optionssol);
    
    % Solve
    diary('GenericMRS.txt');
    sol = opti.solve();
    diary off
    
    % Extract results
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
    
    % Save results
    Ntot = 0;
    for trial = 1:nTrials
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        N = round((tf-t0)*Misc.Mesh_Frequency);
        % Time grid
        tgrid = linspace(t0,tf,N+1)';
        % Save results
        Time(trial).genericMRS = tgrid;
        MActivation(trial).genericMRS = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        lMtildeopt(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        lM(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(Misc.lOpt',1,length(tgrid));
        MvMtilde(trial).genericMRS = vMtilde_opt(:,Ntot + 1:Ntot + N);
        MExcitation(trial).genericMRS = e_opt(:,Ntot + 1:Ntot + N);
        RActivation(trial).genericMRS = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
        MuscleNames = DatStore.MuscleNames;
        OptInfo = output;
        % Tendon forces from lMtilde
        lMTinterp(trial).genericMRS = DatStore(trial).LMTinterp;
        [TForcetilde_,TForce_] = TendonForce_lMtilde(lMtildeopt(trial).genericMRS',Misc.params,lMTinterp(trial).genericMRS,Misc.Atendon,Misc.shift);
        TForcetilde(trial).genericMRS = TForcetilde_';
        TForce(trial).genericMRS = TForce_';
        Ntot = Ntot + N;
    end
    clear opti a lMtilde e vMtilde aT
    
end


%% Dynamic Optimization - Parameter estimation
% ----------------------------------------------------------------------- %
% Estimate parameters
opti_MTE = casadi.Opti();
% States
%   - Muscle activations
a = opti_MTE.variable(DatStore(1).NMuscles,N_tot+nTrials);      % Variable at mesh points
opti_MTE.subject_to(a_min < a < a_max);           % Bounds
%   - Muscle fiber lengths
lMtilde = opti_MTE.variable(DatStore(1).NMuscles,N_tot+nTrials);
opti_MTE.subject_to(lMtilde_min < lMtilde < lMtilde_max);

% Controls
%   - Muscle excitations
e = opti_MTE.variable(DatStore(1).NMuscles,N_tot);
opti_MTE.subject_to(e_min < e < e_max);

params = Misc.params';
lMo = params(:,2);
alphao = params(:,4);
w = lMo.*sin(alphao);
e_lMnuttig = (lMo.^2-w.^2).*opti_MTE.variable(DatStore(1).NMuscles,N_tot);
opti_MTE.subject_to(e_lMnuttig(:) > 0);
%   - Reserve actuators
aT = opti_MTE.variable(DatStore(1).Ndof,N_tot);
opti_MTE.subject_to(-1 < aT <1);
%   - Time derivative of muscle-tendon forces (states)
vMtilde = Misc.scaling.vMtilde.*opti_MTE.variable(DatStore(1).NMuscles,N_tot);
opti_MTE.subject_to(vMtilde_min < vMtilde < vMtilde_max);

% Free optimal fiber length
lMo_scaling_param  = opti_MTE.variable(DatStore(1).NMuscles,1);
lb_lMo_scaling     = ones(DatStore(1).NMuscles,1);   % default upper and lower bound is one (equality constraint)
ub_lMo_scaling     = ones(DatStore(1).NMuscles,1);   % default upper and lower bound is one (equality constraint)
iM                 = DatStore(1).free_lMo(:);        % index muscles with parameter estimation
lb_lMo_scaling(iM) = Misc.lb_lMo_scaling;            % update lower bound for these muscles
ub_lMo_scaling(iM) = Misc.ub_lMo_scaling;            % update uppder bound for these muscles
opti_MTE.subject_to(lb_lMo_scaling < lMo_scaling_param < ub_lMo_scaling); % update the upper and lower bounds

% Free slack length
lTs_scaling_param = opti_MTE.variable(DatStore(1).NMuscles,1);
lb_lTs_scaling = ones(DatStore(1).NMuscles,1); ub_lTs_scaling = ones(DatStore(1).NMuscles,1);
lb_lTs_scaling(DatStore(1).free_lMo(:)) = Misc.lb_lTs_scaling*lb_lTs_scaling(DatStore(1).free_lMo(:));
ub_lTs_scaling(DatStore(1).free_lMo(:)) = Misc.ub_lTs_scaling*ub_lTs_scaling(DatStore(1).free_lMo(:));
opti_MTE.subject_to(lb_lTs_scaling < lTs_scaling_param < ub_lTs_scaling);

% Free tendon stifness
kT_scaling_param = opti_MTE.variable(DatStore(1).NMuscles,1);
lb_kT_scaling_param = ones(DatStore(1).NMuscles,1); ub_kT_scaling_param = ones(DatStore(1).NMuscles,1);
lb_kT_scaling_param(DatStore(1).free_kT(:)) = Misc.lb_kT_scaling*lb_kT_scaling_param(DatStore(1).free_kT(:));
ub_kT_scaling_param(DatStore(1).free_kT(:)) = Misc.ub_kT_scaling*ub_kT_scaling_param(DatStore(1).free_kT(:));
opti_MTE.subject_to(lb_kT_scaling_param < kT_scaling_param < ub_kT_scaling_param);
for k = 1:size(DatStore(1).coupled_kT,1)
    for j = 1:size(DatStore(1).coupled_kT,2)-1
        opti_MTE.subject_to(kT_scaling_param(DatStore(1).coupled_kT(k,j)) - kT_scaling_param(DatStore(1).coupled_kT(k,j+1)) == 0);
    end
end

% Scale factor for EMG
if DatStore(1).EMG.boolEMG
    nEMG        = DatStore(1).EMG.nEMG;
    EMGscale    = opti_MTE.variable(nEMG,1);
    opti_MTE.subject_to(0 < EMGscale < Misc.MaxScaleEMG);
end



        
        
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
    opti_MTE.set_initial(e_lMnuttig , (lMtilde_opt(:,1:N).*lMo).^2 - w.^2);

else
    opti_MTE.set_initial(a,SoActGuess);             % Initial guess (static optimization)
    opti_MTE.set_initial(lMtilde,lMtildeGuess);
    opti_MTE.set_initial(e, SoExcGuess);
    opti_MTE.set_initial(vMtilde,vMtildeGuess);
    opti_MTE.set_initial(aT,SoRActGuess./Misc.Topt);
    opti_MTE.set_initial(lMo_scaling_param,1);
    opti_MTE.set_initial(lTs_scaling_param,1);
    opti_MTE.set_initial(kT_scaling_param,1);
    opti_MTE.set_initial(e_lMnuttig , lMtildeGuess(:,1:N).*lMo.^2 - w.^2 );

end

% get EMG at optimization mesh
for trial = 1:Misc.nTrials
    if DatStore(trial).EMG.boolEMG
        EMGTracking(trial).data = interp1(DatStore(trial).EMG.time,DatStore(trial).EMG.EMGsel,Mesh(trial).t(1:end-1));
    end
end

% get US data at optimization mesh
if ~isempty(Misc.USfile)
    for trial = 1:Misc.nTrials
        DatStore(trial).boolUS = 1;
        USdata =  importdata(Misc.USfile{trial});
        USTracking(trial).data = interp1(DatStore(trial).US.time,DatStore(trial).US.USsel,Mesh(trial).t(1:end));
    end
end



% Loop over mesh points formulating NLP
J = 0; % Initialize cost function
N_acc = 0;
for trial = 1:Misc.nTrials
    % Time bounds    
    N = Mesh(trial).N;
    h = Mesh(trial).step;
    for k=1:N
        % Variables within current mesh interval
        ak = a(:,(N_acc+trial-1) + k); lMtildek = lMtilde(:,(N_acc+trial-1) + k);
        vMtildek = vMtilde(:,N_acc + k); aTk = aT(:,N_acc + k); ek = e(:,N_acc + k);
        e_lMnuttigk = e_lMnuttig(:,(N_acc+trial-1) + k);
       
        % Integration   Uk = (X_(k+1) - X_k)/*dt
        Xk = [ak; lMtildek];
        Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
        Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
        opti_MTE.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
        
        % Get muscle-tendon forces and derive Hill-equilibrium        
        [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree_vTom(ak,lMtildek,vMtildek,e_lMnuttigk,DatStore(trial).LMTinterp(k,:)',[lMo_scaling_param lTs_scaling_param kT_scaling_param],Misc.params',Misc.Atendon');
        params = Misc.params';
        lMo = lMo_scaling_param.*params(:,2);
        alphao = params(:,4);
        lMk = lMtildek.*lMo;
        wk = lMo.*sin(alphao);
        opti_MTE.subject_to(e_lMnuttigk == lMk.^2 - wk.^2);
        % Add path constraints
        % Moment constraints
        for dof = 1:DatStore(trial).Ndof
            T_exp = DatStore(trial).IDinterp(k,dof);
            index_sel = (dof-1)*(DatStore(trial).NMuscles)+1:(dof-1)*(DatStore(trial).NMuscles)+DatStore(trial).NMuscles;
            T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
            opti_MTE.subject_to(T_exp - T_sim == 0);
        end
        % Hill-equilibrium constraint
        opti_MTE.subject_to(Hilldiffk == 0);        
    end

    % tracking lMtilde
    if DatStore(trial).boolUS        
        lMo = lMo_scaling_param(DatStore(trial).free_lMo(:)).*Misc.params(2,DatStore(trial).free_lMo(:))';
        lMtilde_tracking = USTracking(trial).data./lMo/1000; % US data expected in mm in the input file.
        lMtilde_simulated = lMtilde(DatStore(trial).US.USindices,(N_acc+trial:N_acc+trial+N)); 
        J = J + Misc.wlM*sumsqr(lMtilde_simulated-lMtilde_tracking)/DatStore(trial).US.nUS/N;
    end

    % tracking Muscle activity
    if DatStore(trial).EMG.boolEMG
        eSim  = e(DatStore(trial).EMG.EMGindices,N_acc:N_acc+N-1);
        eMeas = EMGTracking(trial).data' .* repmat(EMGscale,1,N);
        JEMG  = sumsqr(eSim-eMeas);
        J     = J + Misc.wEMG * JEMG/DatStore(trial).EMG.nEMG/N;
    end
    N_acc = N_acc + N;
end
J = J + 1e-6*sumsqr(e_lMnuttig);
J = J + ...
    0.5*(sumsqr(e)/N_tot/DatStore(trial).NMuscles + sumsqr(a)/N_tot/DatStore(trial).NMuscles) + ...
    Misc.w1*sumsqr(aT)/N_tot/DatStore(trial).Ndof + ...
    Misc.w2*sumsqr(vMtilde)/N_tot/DatStore(trial).NMuscles;

opti_MTE.minimize(J); % Define cost function in opti
opti_MTE.solver(output.setup.nlp.solver,optionssol);

% Solve
diary('MTE.txt');
% result = solve_NLPSOL(opti_MTE,optionssol)
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

lMo_opt_ = lMo_scaling_param_opt.*Misc.params(2,:)';
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
Ntot = 0;
for trial = 1:nTrials
    t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
    N = round((tf-t0)*Misc.Mesh_Frequency);
    % Time grid
    tgrid = linspace(t0,tf,N+1)';
    % Save results
    Time(trial).MTE = tgrid;
    MActivation(trial).MTE = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    lMtildeopt(trial).MTE = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
    lM(trial).MTE = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(Misc.lOpt',1,length(tgrid));
    MvMtilde(trial).MTE = vMtilde_opt(:,Ntot + 1:Ntot + N);
    MExcitation(trial).MTE = e_opt(:,Ntot + 1:Ntot + N);
    RActivation(trial).MTE = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
    lMo_opt(trial).MTE = lMo_opt_;
    MuscleNames = DatStore(trial).MuscleNames;
    OptInfo = output;
    % Tendon forces from lMtilde
    lMTinterp(trial).MTE = DatStore(trial).LMTinterp;
    [TForcetilde_,TForce_] = TendonForce_lMtilde(lMtildeopt(trial).MTE',Misc.params,lMTinterp(trial).MTE,Misc.Atendon,Misc.shift);
    TForcetilde(trial).MTE = TForcetilde_';
    TForce(trial).MTE = TForce_';
    lMo_scaling_paramopt(trial).MTE = lMo_scaling_param_opt;
    lTs_scaling_paramopt(trial).MTE = lTs_scaling_param_opt;
    kT_scaling_paramopt(trial).MTE = kT_scaling_param_opt;
    Ntot = Ntot + N;
end


clear opti_MTE a lMtilde e vMtilde aT


Parameters = [];

% Collect results
Results.Time = Time;
Results.MActivation = MActivation;
Results.lMtildeopt = lMtildeopt;
Results.lM = lM;
Results.MvMtilde = MvMtilde;
Results.MExcitation = MExcitation;
Results.RActivation = RActivation;
Results.lMo_opt = lMo_opt;
Results.lMTinterp = lMTinterp;
Results.lMo_scaling_paramopt = lMo_scaling_paramopt;
Results.lTs_scaling_paramopt = lTs_scaling_paramopt;
Results.kT_scaling_paramopt = kT_scaling_paramopt;

if Misc.ValidationBool == true
    opti_validation = casadi.Opti();
    
    % Variables - bounds and initial guess
    % States (at mesh and collocation points)
    % Muscle activations
    a = opti_validation.variable(DatStore(trial).NMuscles,N_tot+nTrials);      % Variable at mesh points
    opti_validation.subject_to(a_min < a < a_max);           % Bounds
    opti_validation.set_initial(a,a_opt);             % Initial guess (static optimization)
    % Muscle fiber lengths
    lMtilde = opti_validation.variable(DatStore(trial).NMuscles,N_tot+nTrials);
    opti_validation.subject_to(lMtilde_min < lMtilde < lMtilde_max);
    opti_validation.set_initial(lMtilde,lMtilde_opt);
    
    % Controls
    e = opti_validation.variable(DatStore(trial).NMuscles,N_tot);
    opti_validation.subject_to(e_min < e < e_max);
    opti_validation.set_initial(e, e_opt);
    % Reserve actuators
    aT = opti_validation.variable(DatStore(trial).Ndof,N_tot);
    opti_validation.subject_to(-1 < aT <1);
    opti_validation.set_initial(aT,aT_opt./Misc.Topt);
    % Time derivative of muscle-tendon forces (states)
    vMtilde = Misc.scaling.vMtilde.*opti_validation.variable(DatStore(trial).NMuscles,N_tot);
    opti_validation.subject_to(vMtilde_min < vMtilde < vMtilde_max);
    opti_validation.set_initial(vMtilde,vMtilde_opt);
    
    % Generate optimized parameters for specific trial by scaling generic parameters
    optimized_params = Misc.params';
    optimized_params(:,2) = lMo_scaling_paramopt(trial).MTE.*optimized_params(:,2);
    optimized_params(:,3) = lTs_scaling_paramopt(trial).MTE.*optimized_params(:,3);
    optimized_Atendon = kT_scaling_paramopt(trial).MTE.*Misc.Atendon';
    optimized_shift = (exp(optimized_Atendon.*(1 - 0.995)))/5 - (exp(35.*(1 - 0.995)))/5;
    
    % Loop over mesh points formulating NLP
    J = 0; % Initialize cost function
    N_acc = 0;
    for trial = 1:Misc.nTrials
        % Time bounds
        N = Mesh(trial).N;
        h = Mesh(trial).step;
        
        for k=1:N
            % Variables within current mesh interval
            ak = a(:,(N_acc+trial-1) + k); lMtildek = lMtilde(:,(N_acc+trial-1) + k);
            vMtildek = vMtilde(:,N_acc + k); aTk = aT(:,N_acc + k); ek = e(:,N_acc + k);
            
            % Integration   Uk = (X_(k+1) - X_k)/*dt
            Xk = [ak; lMtildek];
            Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
            Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
            opti_validation.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
            
            % Get muscle-tendon forces and derive Hill-equilibrium
            [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,DatStore(trial).LMTinterp(k,:)',optimized_params,optimized_Atendon,optimized_shift);
            
            % Add path constraints
            % Moment constraints
            for dof = 1:DatStore(trial).Ndof
                T_exp = DatStore(trial).IDinterp(k,dof);
                index_sel = (dof-1)*(DatStore(trial).NMuscles)+1:(dof-1)*(DatStore(trial).NMuscles)+DatStore(trial).NMuscles;
                T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
                opti_validation.subject_to(T_exp - T_sim == 0);
            end
            % Hill-equilibrium constraint
            opti_validation.subject_to(Hilldiffk == 0);
            
        end
        
        J = J + ...
            0.5*(sumsqr(e)/N/DatStore(trial).NMuscles + sumsqr(a)/N/DatStore(trial).NMuscles) + ...
            Misc.w1*sumsqr(aT)/N/DatStore(trial).Ndof + ...
            Misc.w2*sumsqr(vMtilde)/N/DatStore(trial).NMuscles;
        N_acc = N_acc + N;
    end
    opti_validation.minimize(J); % Define cost function in opti
    
    % Create an NLP solver
    opti_validation.solver(output.setup.nlp.solver,optionssol);
    
    % Solve
    diary('ValidationMRS.txt');
    sol = opti_validation.solve();
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
    
    % Save results
    Ntot = 0;
    for trial = 1:nTrials
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        N = round((tf-t0)*Misc.Mesh_Frequency);
        % Time grid
        tgrid = linspace(t0,tf,N+1)';
        % Save results
        Time(trial).validationMRS = tgrid;
        MActivation(trial).validationMRS = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        lMtildeopt(trial).validationMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        lM(trial).validationMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(Misc.lOpt',1,length(tgrid));
        MvMtilde(trial).validationMRS = vMtilde_opt(:,Ntot + 1:Ntot + N);
        MExcitation(trial).validationMRS = e_opt(:,Ntot + 1:Ntot + N);
        RActivation(trial).validationMRS = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
        MuscleNames = DatStore.MuscleNames;
        OptInfo = output;
        % Tendon forces from lMtilde
        lMTinterp(trial).validationMRS = DatStore(trial).LMTinterp;
        [TForcetilde_,TForce_] = TendonForce_lMtilde(lMtildeopt(trial).validationMRS',optimized_params',lMTinterp(trial).validationMRS,optimized_Atendon',optimized_shift');
        TForcetilde(trial).validationMRS = TForcetilde_';
        TForce(trial).validationMRS = TForce_';
        Ntot = Ntot + N;
    end
    clear opti_validation a lMtilde e vMtilde aT
end

% Collect results
Results.Time = Time;
Results.MActivation = MActivation;
Results.lMtildeopt = lMtildeopt;
Results.lM = lM;
Results.MvMtilde = MvMtilde;
Results.MExcitation = MExcitation;
Results.RActivation = RActivation;
Results.lMTinterp = lMTinterp;

% Plot US tracking
figure(5)
for trial = 1:nTrials
    subplot(nTrials,1,trial)
    plot(Time(trial).MTE,lMtildeopt(trial).MTE(DatStore(trial).free_lMo(:),:).*lMo_opt(trial).MTE(DatStore(trial).free_lMo(:)),'LineWidth',2); hold on;
    if Misc.MRSBool == 1
        plot(Time(trial).MTE,lMtildeopt(trial).genericMRS(DatStore(trial).free_lMo(:),:).*Misc.params(2,DatStore(trial).free_lMo(:))','LineWidth',2); hold on;
    end
    if Misc.ValidationBool == 1
        plot(Time(trial).validationMRS,lMtildeopt(trial).validationMRS(DatStore(trial).free_lMo(:),:).*optimized_params(DatStore(trial).free_lMo(:),2),'LineWidth',2); hold on;
    end
    plot(Time(trial).MTE,USTracking(trial).data/1000,'LineWidth',2);
    if Misc.MRSBool == 1 && Misc.ValidationBool == 1        
        legend('MTE','Generic MRS','Optimized MRS','USdata');
    elseif Misc.MRSBool == 1
        legend('MTE','MRS','USdata');
    elseif Misc.ValidationBool == 1
        legend('MTE','Optimized MRS','USdata');
    else
        legend('MTE','USdata');
    end
end


end

