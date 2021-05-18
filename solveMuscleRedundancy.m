function [Results,DatStore,Misc] = solveMuscleRedundancy(model_path,time,OutPath,Misc)
% -----------------------------------------------------------------------%
% INPUTS:
%           model_path: path to the .osim model        
%           time: time window
%           OutPath: path to folder where results will be saved
%           Misc: structure with general input data (see manual for more details)
%
% OUTPUTS:
%           Results:    structure with outputs (states, controls, ...)
%           DatStore:   structure with data used for solving the optimal control problem
%           Misc:       structure with general input data (see manual for more details)
% -----------------------------------------------------------------------%

% add the model name to the misc structure (mainly for post processing)
Misc.model_path = model_path;

% update default settings
Misc = DefaultSettings(Misc);

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%check if we have to adapt the start and end time so that it corresponds to
%time frames in the IK solution
[time] = Check_TimeIndices(Misc,time);
Misc.time=time;



%% Extract muscle information
% ----------------------------------------------------------------------- %
% Perform muscle analysis for the different selected trials
DatStore = struct;
for i = 1:Misc.nTrials
    % select the IK and ID file
    IK_path_trial = Misc.IKfile{i};
    ID_path_trial = Misc.IDfile{i};
    % Run muscle analysis    
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
    if ~isfield(Misc,'kT') || isempty(Misc.kT)
        Misc.kT =ones(1,length(Misc.MuscleNames_Input)).*35;
    end
    % Shift tendon force-length curve as a function of the tendon stiffness
    Misc.shift = getShift(Misc.kT);
    [DatStore] = getMuscleInfo(IK_path_trial,ID_path_trial,Misc,DatStore,i);
    
    % display warnings in muscle selection
    Warnings_MuscleNames(DatStore,Misc,i);
    
    % get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
    [DatStore] = GetIndices_US(DatStore,Misc,i);
end

% set the default value of the tendon stiffness
if isfield(Misc,'Set_kT_ByName') && ~isempty(Misc.Set_kT_ByName)
    [Misc,DatStore] = set_kT_ByName(Misc,DatStore);
end

% get the EMG information
[DatStore] = GetEMGInfo(Misc,DatStore);
[DatStore] = GetUSInfo(Misc,DatStore);

% get the number of muscles
NMuscles = length(DatStore(1).MuscleNames);

%% Static optimization
% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization
% NOTE: We do not estimate any parameters here, but these results can serve as
% decent initial guess for the later dynamic optimization
% Extract the muscle-tendon properties
[Misc.params,Misc.lMo,Misc.lTs,Misc.FMo,Misc.alphao]=ReadMuscleParameters(model_path,DatStore(1).MuscleNames);

% Static optimization using IPOPT solver (used as an initial guess)
for trial = 1:Misc.nTrials
    DatStore = SolveStaticOptimization_IPOPT_CasADi(DatStore,Misc,trial);
end

%% Input activation and contraction dynamics
% ----------------------------------------------------------------------- %
tau_act = 0.015;    Misc.tauAct = tau_act * ones(NMuscles, 1);       % activation time constant (activation dynamics)
tau_deact = 0.06;   Misc.tauDeact = tau_deact * ones(NMuscles,1);  % deactivation time constant (activation dynamics)
Misc.b = 0.1;       % tanh coefficient for smooth activation dynamics

Misc.kT=Misc.kT;
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

for trial = 1:Misc.nTrials
    % Discretization
    N = Mesh(trial).N;
    time_opt = Mesh(trial).t;    
    % Spline approximation of muscle-tendon length (LMT), moment arms (MA) and inverse dynamic torques (ID)
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles
            DatStore(trial).JointMASpline(dof).Muscle(m) = spline(DatStore(trial).time,squeeze(DatStore(trial).dM(:,dof,m)));
        end
        DatStore(trial).JointIDSpline(dof) = spline(DatStore(trial).time,DatStore(trial).T_exp(:,dof));
    end
    
    for m = 1:NMuscles
        DatStore(trial).LMTSpline(m) = spline(DatStore(trial).time,DatStore(trial).LMT(:,m));
    end
    
    % Evaluate LMT, VMT, MA and ID at optimization mesh
    DatStore(trial).LMTinterp = zeros(length(time_opt),NMuscles); % Muscle-tendon length
    for m = 1:NMuscles
        [DatStore(trial).LMTinterp(:,m),~,~] = SplineEval_ppuval(DatStore(trial).LMTSpline(m),time_opt,1);
    end
    DatStore(trial).MAinterp = zeros(length(time_opt),DatStore(trial).nDOF*NMuscles); % Moment arm
    DatStore(trial).IDinterp = zeros(length(time_opt),DatStore(trial).nDOF); % Inverse dynamic torque
    for dof = 1:DatStore(trial).nDOF
        for m = 1:NMuscles
            index_sel=(dof-1)*NMuscles+m;
            DatStore(trial).MAinterp(:,index_sel) = ppval(DatStore(trial).JointMASpline(dof).Muscle(m),time_opt);
        end
        DatStore(trial).IDinterp(:,dof) = ppval(DatStore(trial).JointIDSpline(dof),time_opt);
    end
    
    % Interpolate results of static optimization
    DatStore(trial).SoActInterp = interp1(DatStore(trial).time,DatStore(trial).SoAct,time_opt');
    DatStore(trial).SoRActInterp = interp1(DatStore(trial).time,DatStore(trial).SoRAct,time_opt');
    DatStore(trial).SoForceInterp = interp1(DatStore(trial).time,DatStore(trial).SoForce.*DatStore(trial).cos_alpha./Misc.FMo,time_opt);
    [~,DatStore(trial).lMtildeInterp ] = FiberLength_Ftilde(DatStore(trial).SoForceInterp,Misc.params,DatStore(trial).LMTinterp,Misc.kT,Misc.shift);
    DatStore(trial).vMtildeinterp = zeros(size(DatStore(trial).lMtildeInterp));
    for m = 1:NMuscles
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

% Problem bounds
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

% CasADi setup
import casadi.*
opti    = casadi.Opti();    % create opti structure

% get total number of mesh points
nTrials = Misc.nTrials;
N_tot = sum([Mesh().N]);

% set intial guess based on static opt data
SoActGuess = zeros(NMuscles,N_tot);
SoExcGuess = zeros(NMuscles,N_tot-nTrials);
lMtildeGuess = zeros(NMuscles,N_tot);
vMtildeGuess = zeros(NMuscles,N_tot-nTrials);
SoRActGuess = zeros(DatStore(1).nDOF,N_tot-nTrials);
ctx = 1;  ctu= 1;
for trial = 1:nTrials
    ctx_e = ctx+Mesh(trial).N;      % counter for states
    ctu_e = ctu+Mesh(trial).N-1;    % counter for controls
    SoActGuess(:,ctx:ctx_e) = DatStore(trial).SoActInterp';
    SoExcGuess(:,ctu:ctu_e) = DatStore(trial).SoActInterp(1:end-1,:)';
    lMtildeGuess(:,ctx:ctx_e) = DatStore(trial).lMtildeInterp';
    vMtildeGuess(:,ctu:ctu_e) = DatStore(trial).vMtildeinterp(1:end-1,:)';
    SoRActGuess(:,ctu:ctu_e) = DatStore(trial).SoRActInterp(1:end-1,:)';
    ctx = ctx_e+1;
    ctu = ctu_e+1;
end

if Misc.MRSBool == 1
    % States
    %   - muscle activation
    a = opti.variable(NMuscles,N_tot+nTrials);      % Variable at mesh points
    opti.subject_to(a_min < a < a_max);             % Bounds
    opti.set_initial(a,SoActGuess);                 % Initial guess (static optimization)
    %   - Muscle fiber lengths
    lMtilde = opti.variable(NMuscles,N_tot+nTrials);
    opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
    opti.set_initial(lMtilde,lMtildeGuess);
    %   - Controls
    e = opti.variable(NMuscles,N_tot);
    opti.subject_to(e_min < e < e_max);
    opti.set_initial(e, SoExcGuess);
    %   - Reserve actuators
    aT = opti.variable(DatStore(trial).nDOF,N_tot);
    opti.subject_to(-1 < aT <1);
    %   - Time derivative of muscle-tendon forces (states)
    vMtilde = opti.variable(NMuscles,N_tot);
    opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
    opti.set_initial(vMtilde,vMtildeGuess);    
    %   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
    lM_projected = opti.variable(NMuscles,N_tot + nTrials);
    opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
    % Initial guess for this variable is retrieved from lMtilde guess
    % and geometric relationship between pennation angle, muscle length
    % and width
    lMo = Misc.params(2,:)';
    alphao = Misc.params(4,:)';
    lMGuess = lMtildeGuess.*lMo;
    w = lMo.*sin(alphao);
    lM_projectedGuess = sqrt((lMGuess.^2 - w.^2));
    opti.set_initial(lM_projected,lM_projectedGuess);
    
    % constraint on projected fiber length
    w = lMo.*sin(alphao);  
    lM = lMtilde.*lMo;
    opti.subject_to(lM.^2 - w.^2 == lM_projected.^2);
    
    N_acc = 0; % Index that keeps track of trials that are accumulated         
    % Loop over trials --> one simulation for each trial
    for trial = 1:Misc.nTrials
        % Time bounds
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        % Discretization
        N = Mesh(trial).N;
        h = Mesh(trial).step;
        
        % Loop over mesh points formulating NLP
        for k=1:N
            % Variables within current mesh interval
            ak = a(:,(N_acc+trial-1) + k); lMtildek = lMtilde(:,(N_acc+trial-1) + k);
            vMtildek = vMtilde(:,N_acc + k); aTk = aT(:,N_acc + k); ek = e(:,N_acc + k);
            lM_projectedk = lM_projected(:,(N_acc+trial-1) + k);
            
            % Euler integration  Uk = (X_(k+1) - X_k)/*dt
            Xk = [ak; lMtildek];
            Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
            Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
            opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
            
            % Get muscle-tendon forces and derive Hill-equilibrium
            [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,lM_projectedk,...
                DatStore(trial).LMTinterp(k,:)',Misc.params',Misc.kT',Misc.shift');
            
            % Add path constraints
            % Moment constraints
            for dof = 1:DatStore(trial).nDOF
                T_exp = DatStore(trial).IDinterp(k,dof);
                index_sel = (dof-1)*(NMuscles)+1:(dof*NMuscles); % moment is a vector with the different dofs "below" each other
                T_sim = DatStore(trial).MAinterp(k,index_sel)*FTk + Misc.Topt*aTk(dof);
                opti.subject_to(T_exp - T_sim == 0);
            end
            % Hill-equilibrium constraint
            opti.subject_to(Hilldiffk == 0);
        end
        N_acc = N_acc + N;
    end
    % Cost function
    J = 0.5*(sumsqr(e)/N/NMuscles + sumsqr(a)/N/NMuscles) + ...
        Misc.wTres*sumsqr(aT)/N/DatStore(trial).nDOF + ...
        Misc.wVm*sumsqr(vMtilde)/N/NMuscles;
    
    opti.minimize(J); % Define cost function in opti
    
    % Create an NLP solver
    opti.solver(output.setup.nlp.solver,optionssol);
    
    % Solve
    diary(fullfile(OutPath,[Misc.OutName 'GenericMRS.txt']));
    tic
    sol = opti.solve();
    dt = toc;
    disp(['Computation time solving OCP: ' num2str(dt) ' s'])
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
    % Optimal lM_projectedilary variable
    lM_projected_opt = sol.value(lM_projected);
    
    % Append results to output structures
    Ntot = 0;
    for trial = 1:nTrials
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        N = round((tf-t0)*Misc.Mesh_Frequency);
        % Time grid
        tgrid = linspace(t0,tf,N+1)';
        % Save results
        Results.Time(trial).genericMRS = tgrid;
        Results.MActivation(trial).genericMRS = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        Results.lMtildeopt(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        Results.lM(trial).genericMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(Misc.lMo',1,length(tgrid));
        Results.vMtilde(trial).genericMRS = vMtilde_opt(:,Ntot + 1:Ntot + N);
        Results.MExcitation(trial).genericMRS = e_opt(:,Ntot + 1:Ntot + N);
        Results.RActivation(trial).genericMRS = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
        Results.MuscleNames = DatStore.MuscleNames;
        Results.OptInfo = output;
        % Tendon force
        Results.lMTinterp(trial).genericMRS = DatStore(trial).LMTinterp';
        [TForcetilde_,TForce_] = TendonForce_lMtilde(Results.lMtildeopt(trial).genericMRS',Misc.params,Results.lMTinterp(trial).genericMRS',Misc.kT,Misc.shift);    
        Results.TForcetilde(trial).genericMRS = TForcetilde_';
        Results.TForce(trial).genericMRS = TForce_';
        % get information F/l and F/v properties
        [Fpe_,FMltilde_,FMvtilde_] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).genericMRS',Results.vMtilde(trial).genericMRS',Misc.params(5,:));
        FMo = ones(N+1,1)*Misc.params(1,:);
        Results.Fpe(trial).genericMRS = (Fpe_.*FMo)';
        Results.FMltilde(trial).genericMRS = FMltilde_';
        Results.FMvtilde(trial).genericMRS = FMvtilde_';
        Ntot = Ntot + N;
    end
    clear opti a lMtilde e vMtilde aT
    
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
    % Estimate parameters
    opti_MTE = casadi.Opti();
    % States
    %   - Muscle activations
    a = opti_MTE.variable(NMuscles,N_tot+nTrials);      % Variable at mesh points
    opti_MTE.subject_to(a_min < a < a_max);           % Bounds
    %   - Muscle fiber lengths
    lMtilde = opti_MTE.variable(NMuscles,N_tot+nTrials);
    opti_MTE.subject_to(lMtilde_min < lMtilde < lMtilde_max);
    
    % Controls
    %   - Muscle excitations
    e = opti_MTE.variable(NMuscles,N_tot);
    opti_MTE.subject_to(e_min < e < e_max);
    %   - Reserve actuators
    aT = opti_MTE.variable(DatStore(1).nDOF,N_tot);
    opti_MTE.subject_to(-1 < aT <1);
    %   - Time derivative of muscle-tendon forces (states)
    vMtilde = opti_MTE.variable(NMuscles,N_tot);
    opti_MTE.subject_to(vMtilde_min < vMtilde < vMtilde_max);
    
    %   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
    lM_projected = opti_MTE.variable(NMuscles,N_tot+nTrials);
    opti_MTE.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber has strict positive length
    
    % Free optimal fiber length
    lMo_scaling_param  = opti_MTE.variable(NMuscles,1);
    lb_lMo_scaling     = ones(NMuscles,1);   % default upper and lower bound is one (equality constraint)
    ub_lMo_scaling     = ones(NMuscles,1);   % default upper and lower bound is one (equality constraint)
    iM                 = DatStore(1).free_lMo(:);        % index muscles with parameter estimation
    if ~isempty(iM)
        lb_lMo_scaling(iM) = Misc.lb_lMo_scaling;            % update lower bound for these muscles
        ub_lMo_scaling(iM) = Misc.ub_lMo_scaling;            % update uppder bound for these muscles
    end
    opti_MTE.subject_to(lb_lMo_scaling < lMo_scaling_param < ub_lMo_scaling); % update the upper and lower bounds
    
    % Free slack length
    lTs_scaling_param = opti_MTE.variable(NMuscles,1);    
    lb_lTs_scaling = ones(NMuscles,1); 
    ub_lTs_scaling = ones(NMuscles,1);
    iM             = DatStore(1).free_lMo(:);        % index muscles with parameter estimation
    if ~isempty(iM)
        lb_lTs_scaling(iM) =  Misc.lb_lTs_scaling;
        ub_lTs_scaling(iM) =  Misc.ub_lTs_scaling;
    end
    opti_MTE.subject_to(lb_lTs_scaling < lTs_scaling_param < ub_lTs_scaling);
    
    % Free tendon stiffness
    kT_scaling_param = opti_MTE.variable(NMuscles,1);
    lb_kT_scaling_param = ones(NMuscles,1); 
    ub_kT_scaling_param = ones(NMuscles,1);
    iM = DatStore(1).free_kT(:);
    if ~isempty(iM)
        lb_kT_scaling_param(iM) = Misc.lb_kT_scaling;
        ub_kT_scaling_param(iM) = Misc.ub_kT_scaling;
    end
    opti_MTE.subject_to(lb_kT_scaling_param < kT_scaling_param < ub_kT_scaling_param);
    
    % added tendon stiffness coupling
    for k = 1:size(DatStore(1).coupled_kT,1)
        for j = 1:size(DatStore(1).coupled_kT,2)-1
            opti_MTE.subject_to(kT_scaling_param(DatStore(1).coupled_kT(k,j)) - kT_scaling_param(DatStore(1).coupled_kT(k,j+1)) == 0);
        end
    end    
    % added fibre length coupling
    for k = 1:size(DatStore(1).coupled_lMo,1)
        for j = 1:size(DatStore(1).coupled_lMo,2)-1
            opti_MTE.subject_to(lMo_scaling_param(DatStore(1).coupled_lMo(k,j)) - lMo_scaling_param(DatStore(1).coupled_lMo(k,j+1)) == 0);
        end
    end
    % added tendon slack length coupling
    for k = 1:size(DatStore(1).coupled_lTs,1)
        for j = 1:size(DatStore(1).coupled_lTs,2)-1
            opti_MTE.subject_to(lTs_scaling_param(DatStore(1).coupled_lTs(k,j)) - lTs_scaling_param(DatStore(1).coupled_lTs(k,j+1)) == 0);
        end
    end    
    % Scale factor for EMG
    if DatStore(1).EMG.boolEMG
        nEMG        = DatStore(1).EMG.nEMG;
        EMGscale    = opti_MTE.variable(nEMG,1);
        opti_MTE.subject_to(Misc.BoundsScaleEMG(1) < EMGscale < Misc.BoundsScaleEMG(2));
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
        opti_MTE.set_initial(lM_projected,lM_projected_opt);
    else
        opti_MTE.set_initial(a,SoActGuess);             % Initial guess (static optimization)
        opti_MTE.set_initial(lMtilde,lMtildeGuess);
        opti_MTE.set_initial(e, SoExcGuess);
        opti_MTE.set_initial(vMtilde,vMtildeGuess);
        opti_MTE.set_initial(aT,SoRActGuess./Misc.Topt);
        opti_MTE.set_initial(lMo_scaling_param,1);
        opti_MTE.set_initial(lTs_scaling_param,1);
        opti_MTE.set_initial(kT_scaling_param,1);
        % Hill-type muscle model: geometric relationships
        lMo = Misc.params(2,:)';
        alphao = Misc.params(4,:)';
        lMGuess = lMtildeGuess.*lMo;
        w = lMo.*sin(alphao);
        lM_projectedGuess = sqrt((lMGuess.^2 - w.^2));
        opti_MTE.set_initial(lM_projected,lM_projectedGuess);
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
            DatStore(trial).USTracking = interp1(DatStore(trial).US.time,DatStore(trial).US.USsel,Mesh(trial).t(1:end));
        end
    end
    
    % constraint on projected fiber length
    lMo = lMo_scaling_param.*Misc.params(2,:)';
    alphao = Misc.params(4,:)';
    w = lMo.*sin(alphao);  
    lM = lMtilde.*lMo;
    opti_MTE.subject_to(lM.^2 - w.^2 == lM_projected.^2);
    
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
            lM_projectedk = lM_projected(:,(N_acc+trial-1) + k);
            
            % Euler integration  Uk = (X_(k+1) - X_k)/*dt
            Xk = [ak; lMtildek];
            Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
            Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
            opti_MTE.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
            
            % Get muscle-tendon forces and derive Hill-equilibrium
            [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree(ak,lMtildek,vMtildek,lM_projectedk,...
                DatStore(trial).LMTinterp(k,:)',[lMo_scaling_param lTs_scaling_param kT_scaling_param],Misc.params',Misc.kT');
                                   
            % Add path constraints
            % Moment constraints
            for dof = 1:DatStore(trial).nDOF
                T_exp = DatStore(trial).IDinterp(k,dof);
                index_sel = (dof-1)*(NMuscles)+1:dof*NMuscles;
                T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
                opti_MTE.subject_to(T_exp - T_sim == 0);
            end
            % Hill-equilibrium constraint
            opti_MTE.subject_to(Hilldiffk == 0);
        end        
        % tracking lMtilde
        if DatStore(trial).US.boolUS
            lMo = lMo_scaling_param(DatStore(trial).USsel)'.*Misc.params(2,DatStore(trial).USsel(:));
            lMo = ones(size(DatStore(trial).USTracking,1),1)*lMo;
            lMtilde_tracking = DatStore(trial).USTracking./lMo/1000; % US data expected in mm in the input file.
            lMtilde_simulated = lMtilde(DatStore(trial).US.USindices,(N_acc+trial:N_acc+trial+N));
            if size(lMtilde_simulated,1) ~= size(lMtilde_tracking,1)
                lMtilde_simulated = lMtilde_simulated';
            end
            J = J + Misc.wlM*sumsqr(lMtilde_simulated-lMtilde_tracking)/DatStore(trial).US.nUS/N;
        end        
        % tracking Muscle activity
        if DatStore(trial).EMG.boolEMG
            eSim  = e(DatStore(trial).EMG.EMGindices,N_acc:N_acc+N-1);
            eMeas = EMGTracking(trial).data' .* repmat(EMGscale,1,N);
            JEMG  = sumsqr(eSim-eMeas);
            J     = J + Misc.wEMG * JEMG/DatStore(trial).EMG.nEMG/N;
        end        
        % Bounds on difference EMG - activations
        if DatStore(trial).EMG.boolEMG && ~isempty(Misc.EMGbounds)
            eSim  = e(DatStore(trial).EMG.EMGindices,N_acc:N_acc+N-1);
            eMeas = EMGTracking(trial).data' .* repmat(EMGscale,1,N);
            if Misc.EMGbounds(1) == 0
                opti_MTE.subject_to(eSim == eMeas); % "EMG driven"
            else
                opti_MTE.subject_to(Misc.EMGbounds(1) < eSim-eMeas < Misc.EMGbounds(2)); % "EMG constrained"
            end
            
        end
        N_acc = N_acc + N;
    end
    
    J = J + ...
        Misc.wAct*0.5*(sumsqr(e)/N_tot/NMuscles + sumsqr(a)/N_tot/NMuscles) + ...
        Misc.wTres*sumsqr(aT)/N_tot/DatStore(trial).nDOF + ...
        Misc.wVm*sumsqr(vMtilde)/N_tot/NMuscles;
    
    opti_MTE.minimize(J); % Define cost function in opti
    opti_MTE.solver(output.setup.nlp.solver,optionssol);
    diary(fullfile(OutPath,[Misc.OutName 'MTE.txt']));
    tic
    
    % Note: we don't use opti.solve() here because opti does not
    % recognisize a constraint as an equality constraint if the lower and
    % upper bounds are zero. Threfore, we convert is first to a "standard
    % NLP in casadi. This makes mainly the extraction of the results a bit 
    % more difficult (see below).
    
    %sol = opti_MTE.solve();
    [w_opt] = solve_NLPSOL(opti_MTE,optionssol);
    dt = toc;
    disp(['Computation time solving OCP: ' num2str(dt) ' s'])
    diary off
    
    %% Extract results
    % Variables at mesh points
    % Muscle activations and muscle-tendon forces
    
    % get the results from w_opt vector (this is vector with all the
    % variables (constructed with opti.variable)
    ct = 0; % help for extracting results
    iSel = ct+1:ct+NMuscles*(N_tot+nTrials);
    a_opt = reshape(w_opt(iSel),[NMuscles,(N_tot+nTrials)]);
    ct = ct + NMuscles*(N_tot+nTrials);
    
    iSel = ct+1:ct+NMuscles*(N_tot+nTrials);
    lMtilde_opt = reshape(w_opt(iSel),[NMuscles,(N_tot+nTrials)]);
    ct = ct + NMuscles*(N_tot+nTrials);
    
    iSel = ct+1:ct+NMuscles*N_tot;
    e_opt = reshape(w_opt(iSel),[NMuscles,(N_tot)]);
    ct = ct + NMuscles*(N_tot);
    
    iSel = ct+1:ct+DatStore(1).nDOF*N_tot;
    aT_opt = reshape(w_opt(iSel),[DatStore(1).nDOF,(N_tot)]);
    ct = ct + DatStore(1).nDOF*(N_tot);
    
    iSel = ct+1:ct+NMuscles*N_tot;
    vMtilde_opt = reshape(w_opt(iSel),[NMuscles,(N_tot)]);
    ct = ct + NMuscles*(N_tot);
    
    iSel = ct+1:ct+NMuscles*(N_tot+nTrials);
    lM_projected_opt = reshape(w_opt(iSel),[NMuscles,(N_tot+nTrials)]);
    ct = ct + NMuscles*(N_tot+nTrials); 
    
    iSel = ct+1:ct+NMuscles;
    lMo_scaling_param_opt = reshape(w_opt(iSel),[NMuscles 1]);
    ct = ct + NMuscles; 
    
    iSel = ct+1:ct+NMuscles;
    lTs_scaling_param_opt = reshape(w_opt(iSel),[NMuscles 1]);
    ct = ct + NMuscles; 
    
    iSel = ct+1:ct+NMuscles;
    kT_scaling_param_opt = reshape(w_opt(iSel),[NMuscles 1]);
    ct = ct + NMuscles; 
    
    if DatStore(1).EMG.boolEMG
        iSel = ct+1:ct+nEMG;
        EMGscale_opt = reshape(w_opt(iSel),[nEMG 1]);
        ct = ct + nEMG; 
    end
    
    % adapted optimal fiber length
    lMo_opt_ = lMo_scaling_param_opt.*Misc.params(2,:)';
    
    if ct ~= length(w_opt)
        error('Problem with extracting results from NLP solution in parameter estimation');
    end
    
    Ntot = 0;
    for trial = 1:nTrials
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        N = round((tf-t0)*Misc.Mesh_Frequency);
        % Time grid
        tgrid = linspace(t0,tf,N+1)';
        % Save results
        Results.Time(trial).MTE = tgrid;
        Results.MActivation(trial).MTE = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        Results.lMtildeopt(trial).MTE = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        Results.lM(trial).MTE = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(lMo_opt_,1,length(tgrid));
        Results.vMtilde(trial).MTE = vMtilde_opt(:,Ntot + 1:Ntot + N);
        Results.MExcitation(trial).MTE = e_opt(:,Ntot + 1:Ntot + N);
        Results.RActivation(trial).MTE = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
        % Update muscle params      
        kTOpt = Misc.kT .* kT_scaling_param_opt';
        shiftOpt = getShift(kTOpt);        
        paramsOpt =  Misc.params;
        paramsOpt(2,:) = paramsOpt(2,:).*lMo_scaling_param_opt'; % updated optimal fiber length
        paramsOpt(3,:) = paramsOpt(3,:).*lTs_scaling_param_opt';  % updated tendon slack length       
        Results.lMTinterp(trial).MTE = DatStore(trial).LMTinterp';
        % get force from muscle state
        [TForcetilde_,TForce_] = TendonForce_lMtilde(Results.lMtildeopt(trial).MTE',paramsOpt,Results.lMTinterp(trial).MTE',kTOpt,shiftOpt);
        Results.TForcetilde(trial).MTE = TForcetilde_';
        Results.TForce(trial).MTE = TForce_';
        [Fpe_,FMltilde_,FMvtilde_] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).MTE',Results.vMtilde(trial).MTE',Misc.params(5,:));
        FMo = ones(N+1,1)*Misc.params(1,:);
        Results.Fpe(trial).MTE = (Fpe_.*FMo)';
        Results.FMltilde(trial).MTE = FMltilde_';
        Results.FMvtilde(trial).MTE = FMvtilde_';
        Ntot = Ntot + N;
    end
else
    lMo_scaling_param_opt = ones(NMuscles,1);
    lTs_scaling_param_opt = ones(NMuscles,1);
    kT_scaling_param_opt = ones(NMuscles,1);
end

clear opti_MTE a lMtilde e vMtilde aT


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
    Results.Param.Estimated.FMo     = Results.Param.Original.FMo;
    Results.Param.Estimated.lMo     = Results.Param.Original.lMo .* Results.Param.lMo_scaling_paramopt';
    Results.Param.Estimated.lTs  = Results.Param.Original.lTs .* Results.Param.lTs_scaling_paramopt';
    Results.Param.Estimated.alphao = Results.Param.Original.alphao;
    Results.Param.Estimated.kT  = Results.Param.Original.kT .* Results.Param.kT_scaling_paramopt';
    Results.Param.Bound.lMo.lb       = Misc.lb_lMo_scaling;
    Results.Param.Bound.lMo.ub       = Misc.ub_lMo_scaling;
    Results.Param.Bound.lTs.lb       = Misc.lb_lTs_scaling;
    Results.Param.Bound.lTs.ub       = Misc.ub_lTs_scaling;
    Results.Param.Bound.kT.lb        = Misc.lb_kT_scaling;
    Results.Param.Bound.kT.ub        = Misc.ub_kT_scaling;
    Results.Param.Bound.EMG.lb       = Misc.BoundsScaleEMG(1);
    Results.Param.Bound.EMG.ub       = Misc.BoundsScaleEMG(2);
    if DatStore(1).EMG.boolEMG
        Results.Param.EMGscale       = EMGscale_opt;
    end
end

% store the Misc structure as well in the results
Results.Misc = Misc;

%% Validate results parameter estimation 
if Misc.ValidationBool == true && BoolParamOpt
    opti_validation = casadi.Opti();
    
    % Variables - bounds and initial guess
    % States (at mesh and collocation points)
    % Muscle activations
    a = opti_validation.variable(NMuscles,N_tot+nTrials);      % Variable at mesh points
    opti_validation.subject_to(a_min < a < a_max);           % Bounds
    opti_validation.set_initial(a,a_opt);             % Initial guess (static optimization)
    % Muscle fiber lengths
    lMtilde = opti_validation.variable(NMuscles,N_tot+nTrials);
    opti_validation.subject_to(lMtilde_min < lMtilde < lMtilde_max);
    opti_validation.set_initial(lMtilde,lMtilde_opt);
    
    % Controls
    e = opti_validation.variable(NMuscles,N_tot);
    opti_validation.subject_to(e_min < e < e_max);
    opti_validation.set_initial(e, e_opt);
    % Reserve actuators
    aT = opti_validation.variable(DatStore(trial).nDOF,N_tot);
    opti_validation.subject_to(-1 < aT <1);
    opti_validation.set_initial(aT,aT_opt);
    % Time derivative of muscle-tendon forces (states)
    vMtilde = opti_validation.variable(NMuscles,N_tot);
    opti_validation.subject_to(vMtilde_min < vMtilde < vMtilde_max);
    opti_validation.set_initial(vMtilde,vMtilde_opt);
    
    % Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
    lM_projected = opti_validation.variable(NMuscles,N_tot+nTrials);
    opti_validation.subject_to(1e-4 < lM_projected(:));
    opti_validation.set_initial(lM_projected,lM_projected_opt);    
    
    % Generate optimized parameters for specific trial by scaling generic parameters
    optimized_params = Misc.params';
    optimized_params(:,2) = lMo_scaling_param_opt.*optimized_params(:,2);
    optimized_params(:,3) = lTs_scaling_param_opt.*optimized_params(:,3);
    optimized_kT = kT_scaling_param_opt.*Misc.kT';
    optimized_shift = getShift(optimized_kT);
    
    % constraint on projected fiber length
    lMo = optimized_params(:,2);
    alphao = optimized_params(:,4);
    lM = lMtilde.*lMo;
    w = lMo.*sin(alphao);
    opti_validation.subject_to(lM.^2 - w.^2 == lM_projected.^2);
    
    
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
            lM_projectedk = lM_projected(:,(N_acc+trial-1) + k);
            
            % Integration   Uk = (X_(k+1) - X_k)/*dt
            Xk = [ak; lMtildek];
            Zk = [a(:,(N_acc+trial-1) + k + 1);lMtilde(:,(N_acc+trial-1) + k + 1)];
            Uk = [ActivationDynamics(ek,ak,Misc.tauAct,Misc.tauDeact,Misc.b); vMtildek];
            opti_validation.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
            
            % Impose that auxilary variable lM_projected behaves as defined
            [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,lM_projectedk,DatStore(trial).LMTinterp(k,:)',optimized_params,optimized_kT,optimized_shift);                  
            
            % Add path constraints
            % Moment constraints
            for dof = 1:DatStore(trial).nDOF
                T_exp = DatStore(trial).IDinterp(k,dof);
                index_sel = (dof-1)*(NMuscles)+1:dof*NMuscles;
                T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
                opti_validation.subject_to(T_exp - T_sim == 0);
            end
            % Hill-equilibrium constraint
            opti_validation.subject_to(Hilldiffk == 0);
            
        end
        % Cost function
        J = J + ...
            0.5*(sumsqr(e)/N/NMuscles + sumsqr(a)/N/NMuscles) + ...
            Misc.wTres*sumsqr(aT)/N/DatStore(trial).nDOF + ...
            Misc.wVm*sumsqr(vMtilde)/N/NMuscles;
        N_acc = N_acc + N;
    end
    opti_validation.minimize(J); % Define cost function in opti
    
    % Create an NLP solver
    opti_validation.solver(output.setup.nlp.solver,optionssol);
    
    % Solve
    diary(fullfile(OutPath,[Misc.OutName 'ValidationMRS.txt']));
    tic
    sol = opti_validation.solve();
    dt = toc;
    disp(['Computation time solving OCP: ' num2str(dt) ' s'])
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
        Results.Time(trial).validationMRS = tgrid;
        Results.MActivation(trial).validationMRS = a_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        Results.lMtildeopt(trial).validationMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1);
        Results.lM(trial).validationMRS = lMtilde_opt(:,(Ntot + trial - 1) + 1:(Ntot + trial - 1) + N + 1).*repmat(optimized_params(:,2),1,length(tgrid));
        Results.vMtilde(trial).validationMRS = vMtilde_opt(:,Ntot + 1:Ntot + N);
        Results.MExcitation(trial).validationMRS = e_opt(:,Ntot + 1:Ntot + N);
        Results.RActivation(trial).validationMRS = aT_opt(:,Ntot + 1:Ntot + N)*Misc.Topt;
        % Tendon forces from lMtilde
        Results.lMTinterp(trial).validationMRS = DatStore(trial).LMTinterp';
        [TForcetilde_,TForce_] = TendonForce_lMtilde(Results.lMtildeopt(trial).validationMRS',optimized_params',Results.lMTinterp(trial).validationMRS',optimized_kT',optimized_shift');
        Results.TForcetilde(trial).validationMRS = TForcetilde_';
        Results.TForce(trial).validationMRS = TForce_';
        [Fpe_,FMltilde_,FMvtilde_] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).MTE',Results.vMtilde(trial).MTE',Misc.params(5,:));
        FMo = ones(N+1,1)*Misc.params(1,:);
        Results.Fpe(trial).validationMRS = (Fpe_.*FMo)';
        Results.FMltilde(trial).validationMRS = FMltilde_';
        Results.FMvtilde(trial).validationMRS = FMvtilde_';
        Ntot = Ntot + N;
    end
    clear opti_validation a lMtilde e vMtilde aT
end
% add selected muscle names to the output structure
Results.MuscleNames = DatStore.MuscleNames;

%% Plot Output
% plot EMG tracking
if Misc.PlotBool && Misc.EMGconstr == 1
    h = PlotEMGTracking(Results,DatStore);
    if ~isdir(fullfile(OutPath,'figures'))
        mkdir(fullfile(OutPath,'figures'));
    end
    saveas(h,fullfile(OutPath,'figures',[Misc.OutName '_fig_EMG.fig']));
end

% plot estimated parameters
if Misc.PlotBool == 1 && BoolParamOpt ==1
    h = PlotEstimatedParameters(Results,DatStore,Misc);
    if ~isdir(fullfile(OutPath,'figures'))
        mkdir(fullfile(OutPath,'figures'));
    end
    saveas(h,fullfile(OutPath,'figures',[Misc.OutName '_fig_Param.fig']));
end

% plot fiber length
if Misc.PlotBool && Misc.UStracking == 1
    h = PlotFiberLength(Results,DatStore);
    if ~isdir(fullfile(OutPath,'figures'))
        mkdir(fullfile(OutPath,'figures'));
    end
    saveas(h,fullfile(OutPath,'figures',[Misc.OutName '_fig_FiberLength.fig']));
end

% plot the states of the muscles in the simulation
if Misc.PlotBool
    h = PlotStates(Results,DatStore,Misc);
    if ~isdir(fullfile(OutPath,'figures'))
        mkdir(fullfile(OutPath,'figures'));
    end
    saveas(h,fullfile(OutPath,'figures',[Misc.OutName '_fig_States.fig']));
end

%% save the results
% plot states and variables from parameter estimation simulation
save(fullfile(OutPath,[Misc.OutName 'Results.mat']),'Results','DatStore','Misc');

% write estimated parameters to new duplicate osim model
if BoolParamOpt
    muscleParams = Results.Param.Estimated;
    muscleNames  = DatStore.MuscleNames;
    modelPath    = char(model_path);
    outPath      = OutPath;
    newModelFile = Misc.newModelFile;
    
    ParamsToOsim(muscleParams,muscleNames,modelPath,outPath,newModelFile); 
end

end


