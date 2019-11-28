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


% ----------------------------------------------------------------------- %
% Solve the muscle redundancy problem using static optimization --------- %

% NOTE: We do not estimate any parameters here, but these results can serve as
% decent initial guess for the later dynamic optimization

% The solution of the static optimization is used as initial guess for the
% dynamic optimization
% Extract the muscle-tendon properties
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=ReadMuscleParameters(model_path,DatStore.MuscleNames);
DatStore.free_lMo = zeros(length(Misc.Estimate_OptFL),1);
for i = 1:length(DatStore.free_lMo)
DatStore.free_lMo(i) = find(strcmp(DatStore.MuscleNames,Misc.Estimate_OptFL{i}));
end
% Static optimization using IPOPT solver
DatStore = SolveStaticOptimization_IPOPT_CasADi(DatStore);

% Show difference in fiber lengths from scaled model with rigid tendon vs
% experimental:

lM_FileName=fullfile(Misc.MuscleAnalysisPath,[Misc.trialName '_MuscleAnalysis_FiberLength.sto']);
lMRigidTendon = importdata(lM_FileName);
lM_index = find(strcmp(lMRigidTendon.colheaders,Misc.FL_expdata));
US_data = importdata(US_path);
figure(2)
plot(lMRigidTendon.data(:,1),1000*lMRigidTendon.data(:,lM_index)); hold on;
plot(US_data.data(:,1),US_data.data(:,2));
legend('GenericMuscleAnalysis','Experimental')


% update Tendon stiffness for specific muscles based on input arguments (we
% should replace this function with bounds (we can use it as inspiration
%if isfield(Misc,'Set_ATendon_ByName') && ~isempty(Misc.Set_ATendon_ByName)
%   [Misc,DatStore] = set_ATendon_ByName(Misc,DatStore);
%end

% get the EMG information
% [DatStore] = GetEMGInfo(Misc,DatStore);

% Extract the muscle-tendon properties
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=...
    ReadMuscleParameters(model_path,DatStore.MuscleNames);

% ....... To be continued......


%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART II: OPTIMAL CONTROL PROBLEM FORMULATION -------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% Input arguments
auxdata.NMuscles = DatStore.nMuscles;   % number of muscles
auxdata.Ndof = DatStore.nDOF;           % number of dofs
auxdata.ID = DatStore.T_exp;            % inverse dynamics
auxdata.params = DatStore.params;       % Muscle-tendon parameters
auxdata.scaling.vMtilde = 10;           % Scaling factor: derivative muscle fiber lengths
auxdata.w1 = 1000;                      % Weight objective function
auxdata.w2 = 0.01;
auxdata.Topt = 150;                     % Scaling factor: reserve actuators

% ADiGator works with 2D: convert 3D arrays to 2D structure (moment arms)
for i = 1:auxdata.Ndof
    auxdata.MA(i).Joint(:,:) = DatStore.dM(:,i,:);  % moment arms
end
auxdata.DOFNames = DatStore.DOFNames;   % names of dofs

tau_act = 0.015; auxdata.tauAct = tau_act * ones(1,auxdata.NMuscles);       % activation time constant (activation dynamics)
tau_deact = 0.06; auxdata.tauDeact = tau_deact * ones(1,auxdata.NMuscles);  % deactivation time constant (activation dynamics)
auxdata.b = 0.1;

% Parameters of active muscle force-velocity characteristic
load('ActiveFVParameters.mat','ActiveFVParameters');
Fvparam(1) = 1.475*ActiveFVParameters(1); Fvparam(2) = 0.25*ActiveFVParameters(2);
Fvparam(3) = ActiveFVParameters(3) + 0.75; Fvparam(4) = ActiveFVParameters(4) - 0.027;
auxdata.Fvparam = Fvparam;

% Parameters of active muscle force-length characteristic
load('Faparam.mat','Faparam');
auxdata.Faparam = Faparam;

% Parameters of passive muscle force-length characteristic
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
auxdata.Fpparam = [pp1;pp2];
auxdata.Atendon=Misc.Atendon;
auxdata.shift=Misc.shift;

% Problem bounds
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;      % bounds on normalized muscle fiber velocity
lMtilde_min = 0.2; lMtilde_max = 1.8;   % bounds on normalized muscle fiber length
% Time bounds
t0 = DatStore.time(1); tf = DatStore.time(end);

% CasADi setup
import casadi.*
opti = casadi.Opti();

% Load CasADi function
CasADiFunctions

% Discretization
N = round((tf-t0)*Misc.Mesh_Frequency);
h = (tf-t0)/N;

% Spline approximation of muscle-tendon length (LMT), moment arms (MA) and inverse dynamic torques (ID)
for dof = 1:auxdata.Ndof
    for m = 1:auxdata.NMuscles
        auxdata.JointMASpline(dof).Muscle(m) = spline(DatStore.time,auxdata.MA(dof).Joint(:,m));
    end
    auxdata.JointIDSpline(dof) = spline(DatStore.time,DatStore.T_exp(:,dof));
end

for m = 1:auxdata.NMuscles
    auxdata.LMTSpline(m) = spline(DatStore.time,DatStore.LMT(:,m));
end

% Evaluate LMT, VMT, MA and ID at optimization mesh
step = (tf-t0)/(N);
time_opt = t0:step:tf;
LMTinterp = zeros(length(time_opt),auxdata.NMuscles); % Muscle-tendon length
for m = 1:auxdata.NMuscles
    [LMTinterp(:,m),~,~] = SplineEval_ppuval(auxdata.LMTSpline(m),time_opt,1);
end
MAinterp = zeros(length(time_opt),auxdata.Ndof*auxdata.NMuscles); % Moment arm
IDinterp = zeros(length(time_opt),auxdata.Ndof); % Inverse dynamic torque
for dof = 1:auxdata.Ndof
    for m = 1:auxdata.NMuscles
        index_sel=(dof-1)*(auxdata.NMuscles)+m;
        MAinterp(:,index_sel) = ppval(auxdata.JointMASpline(dof).Muscle(m),time_opt);
    end
    IDinterp(:,dof) = ppval(auxdata.JointIDSpline(dof),time_opt);
end

% Initial guess static optimization
SoActInterp = interp1(DatStore.time,DatStore.SoAct,time_opt);
SoRActInterp = interp1(DatStore.time,DatStore.SoRAct,time_opt);
SoForceInterp = interp1(DatStore.time,DatStore.SoForce.*DatStore.cos_alpha./DatStore.Fiso,time_opt);
[~,lMtildeInterp ] = FiberLength_Ftilde(SoForceInterp,DatStore.params,LMTinterp,Misc.Atendon,Misc.shift);
vMtildeinterp = zeros(size(lMtildeInterp));
for m = 1:auxdata.NMuscles
    lMtildeSpline = spline(time_opt,lMtildeInterp(:,m));
    [~,vMtildeinterp_norm,~] = SplineEval_ppuval(lMtildeSpline,time_opt,1);
    vMtildeinterp(:,m) = vMtildeinterp_norm/auxdata.scaling.vMtilde;
end


% Variables - bounds and initial guess
% States (at mesh and collocation points)
% Muscle activations
a = opti.variable(auxdata.NMuscles,N+1);      % Variable at mesh points
opti.subject_to(a_min < a < a_max);           % Bounds
opti.set_initial(a,SoActInterp');             % Initial guess (static optimization)
% Muscle fiber lengths
lMtilde = opti.variable(auxdata.NMuscles,N+1);
opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
opti.set_initial(lMtilde,lMtildeInterp');

% Controls
e = opti.variable(auxdata.NMuscles,N);
opti.subject_to(e_min < e < e_max);
opti.set_initial(e, SoActInterp(1:N,:)');
% Reserve actuators
aT = opti.variable(auxdata.Ndof,N);
opti.subject_to(-1 < aT <1);
opti.set_initial(aT,SoRActInterp(1:N,:)'./auxdata.Topt);
% Time derivative of muscle-tendon forces (states)
vMtilde = auxdata.scaling.vMtilde.*opti.variable(auxdata.NMuscles,N);
opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
opti.set_initial(vMtilde,vMtildeinterp(1:N,:)');


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
    opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
    
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffk,FTk] = f_forceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,LMTinterp(k,:)');
    
    % Add path constraints
    % Moment constraints
    for dof = 1:auxdata.Ndof
        T_exp = IDinterp(k,dof);
        index_sel = (dof-1)*(auxdata.NMuscles)+1:(dof-1)*(auxdata.NMuscles)+auxdata.NMuscles;
        T_sim = f_spNMuscles(MAinterp(k,index_sel),FTk) + auxdata.Topt*aTk(dof);
        opti.subject_to(T_exp - T_sim == 0);
    end
    % Hill-equilibrium constraint
    opti.subject_to(Hilldiffk == 0);

end

J = J + ...
    sumsqr(e)/N/auxdata.NMuscles + ...
    auxdata.w1*sumsqr(aT)/N/auxdata.Ndof + ...
    auxdata.w2*sumsqr(vMtilde)/N/auxdata.NMuscles;
    
opti.minimize(J); % Define cost function in opti

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
% vMtilde_opt = sol.value(vMtilde)';

% Grid
% Mesh points
tgrid = linspace(t0,tf,N+1)';
% Save results
Time.genericMRS = tgrid;
MActivation.genericMRS = a_opt;
lMtildeopt.genericMRS = lMtilde_opt;
lM.genericMRS = lMtilde_opt.*repmat(DatStore.lOpt',1,length(Time.genericMRS));
MExcitation.genericMRS = e_opt;
RActivation.genericMRS = aT_opt*auxdata.Topt;
MuscleNames = DatStore.MuscleNames;
OptInfo = output;
% Tendon forces from lMtilde
lMTinterp.genericMRS = LMTinterp;
[TForcetilde_,TForce_] = TendonForce_lMtilde(lMtildeopt.genericMRS',auxdata.params,lMTinterp.genericMRS,auxdata.Atendon,auxdata.shift);
TForcetilde.genericMRS = TForcetilde_';
TForce.genericMRS = TForce_';

clear opti a lMtilde e vMtilde aT

opti = casadi.Opti();

% Variables - bounds and initial guess
% States (at mesh and collocation points)
% Muscle activations
a = opti.variable(auxdata.NMuscles,N+1);      % Variable at mesh points
opti.subject_to(a_min < a < a_max);           % Bounds
opti.set_initial(a,SoActInterp');             % Initial guess (static optimization)
% Muscle fiber lengths
lMtilde = opti.variable(auxdata.NMuscles,N+1);
opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
opti.set_initial(lMtilde,lMtildeInterp');

% Controls
e = opti.variable(auxdata.NMuscles,N);
opti.subject_to(e_min < e < e_max);
opti.set_initial(e, SoActInterp(1:N,:)');
% Reserve actuators
aT = opti.variable(auxdata.Ndof,N);
opti.subject_to(-1 < aT <1);
opti.set_initial(aT,SoRActInterp(1:N,:)'./auxdata.Topt);
% Time derivative of muscle-tendon forces (states)
vMtilde = auxdata.scaling.vMtilde.*opti.variable(auxdata.NMuscles,N);
opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
opti.set_initial(vMtilde,vMtildeinterp(1:N,:)');

% Free optimal fiber length
lMo_scaling_param = opti.variable(auxdata.NMuscles,1);
lb_lMo_scaling = ones(auxdata.NMuscles,1); ub_lMo_scaling = ones(auxdata.NMuscles,1);
lb_lMo_scaling(DatStore.free_lMo(:)) = 0.9*lb_lMo_scaling(DatStore.free_lMo(:));
ub_lMo_scaling(DatStore.free_lMo(:)) = 1.3*ub_lMo_scaling(DatStore.free_lMo(:));
opti.subject_to(lb_lMo_scaling < lMo_scaling_param < ub_lMo_scaling);
opti.set_initial(lMo_scaling_param,1);


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
    opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);
    
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffk,FTk] = f_ForceEquilibrium_lMtildeState_lMoFree(ak,lMtildek,vMtildek,LMTinterp(k,:)',lMo_scaling_param);
    
    % Add path constraints
    % Moment constraints
    for dof = 1:auxdata.Ndof
        T_exp = IDinterp(k,dof);
        index_sel = (dof-1)*(auxdata.NMuscles)+1:(dof-1)*(auxdata.NMuscles)+auxdata.NMuscles;
        T_sim = f_spNMuscles(MAinterp(k,index_sel),FTk) + auxdata.Topt*aTk(dof);
        opti.subject_to(T_exp - T_sim == 0);
    end
    % Hill-equilibrium constraint
    opti.subject_to(Hilldiffk == 0);

end
lMo = lMo_scaling_param(DatStore.free_lMo(:)).*auxdata.params(2,DatStore.free_lMo(:))';
lMtilde_tracking = US_data.data(:,2)'./lMo/1000;
lMtilde_simulated = lMtilde(DatStore.free_lMo(:),:);
J = J + ...
    10*sumsqr(lMtilde_simulated-lMtilde_tracking)/N + ...
    0.5*(sumsqr(e)/N/auxdata.NMuscles + sumsqr(a)/N/auxdata.NMuscles) + ...
    auxdata.w1*sumsqr(aT)/N/auxdata.Ndof + ...
    auxdata.w2*sumsqr(vMtilde)/N/auxdata.NMuscles;
    
opti.minimize(J); % Define cost function in opti

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
lMo_scaling_param_opt = sol.value(lMo_scaling_param);
lMo_opt = lMo_scaling_param_opt.*auxdata.params(2,:)';
% Muscle excitations
e_opt = sol.value(e);
% Reserve actuators
aT_opt = sol.value(aT);
% Time derivatives of muscle-tendon forces
% vMtilde_opt = sol.value(vMtilde)';

% Grid
% Mesh points
tgrid = linspace(t0,tf,N+1)';
figure(5)
plot(tgrid,lMtilde_opt(DatStore.free_lMo(:),:).*lMo_opt(DatStore.free_lMo(:))); hold on;
plot(tgrid,lMtildeopt.genericMRS(DatStore.free_lMo(:),:).*auxdata.params(2,DatStore.free_lMo(:))'); hold on;
plot(tgrid,US_data.data(:,2)'/1000);
legend('MTE','MRS','USdata')

% Save results
Time.MTE = tgrid;
MActivation.MTE = a_opt;
lMtildeopt.MTE = lMtilde_opt;
lM.MTE = lMtilde_opt.*lMo_opt;
MExcitation.MTE = e_opt;
RActivation.MTE = aT_opt*auxdata.Topt;
lMo_scaling_param_opt.MTE = lMo_scaling_param_opt;
MuscleNames = DatStore.MuscleNames;
OptInfo = output;
% Tendon forces from lMtilde
lMTinterp.MTE = LMTinterp;
[TForcetilde,TForce] = TendonForce_lMtilde(lMtildeopt.MTE',auxdata.params,lMTinterp.MTE,auxdata.Atendon,auxdata.shift);
TForcetilde.MTE = TForcetilde';
TForce.MTE = TForce';

clear opti a lMtilde e vMtilde aT

Results.Time = Time;
Results.MActivation = MActivation;
Results.lMtildeopt = lMtildeopt;
Results.lM = lM;
Results.MExcitation = MExcitation;
Results.RActivation = RActivation;
Results.lMo_scaling_param_opt = lMo_scaling_param_opt;
Results.lMTinterp = lMTinterp;
Results.TForcetilde = TForcetilde;
Results.TForce = TForce;

Parameters = [];

end

