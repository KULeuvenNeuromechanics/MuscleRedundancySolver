% SolveMuscleRedundancy_lMtildeState, version 2.1 (November 2018)
%
% This function solves the muscle redundancy problem in the leg using the
% direct collocation optimal control software GPOPS-II as described in De
% Groote F, Kinney AL, Rao AV, Fregly BJ. Evaluation of direct
% collocation optimal control problem formulations for solving the muscle
% redundancy problem. Annals of Biomedical Engineering (2016).
%
% Change with regards to version 0.1: Activation dynamics as described in
% De Goote F, Pipeleers G, Jonkers I, Demeulenaere B, Patten C, Swevers J,
% De Schutter J. A physiology based inverse dynamic analysis of human gait:
% potential and perspectives. Computer Methods in Biomechanics and
% Biomedical Engineering (2009).
%
% Authors:  F. De Groote, M. Afschrift, A. Falisse
% Emails:   friedl.degroote@kuleuven.be
%           maarten.afschrift@kuleuven.be
%           antoine.falisse@kuleuven.be
%
% ----------------------------------------------------------------------- %
% This function uses the normalized muscle fiber length as a state (see
% aforementioned publication for more details)
%
% INPUTS:
%           model_path: path to the .osim model
%           IK_path: path to the inverse kinematics results
%           ID_path: path to the inverse dynamics results
%           time: time window
%           OutPath: path to folder where results will be saved
%           Misc: structure of input data (see manual for more details)
%
% OUTPUTS:
%           Time: time window (as used when solving the optimal control
%           problem)
%           MExcitation: muscle excitation
%           MActivation: muscle activation
%           RActivation: activation of the reserve actuators
%           TForce_tilde: normalized tendon force
%           TForce: tendon force
%           lMtilde: normalized muscle fiber length
%           lM: muscle fiber length
%           MuscleNames: names of muscles
%           OptInfo: output of GPOPS-II
%           DatStore: structure with data used for solving the optimal
%           control problem
%
% ----------------------------------------------------------------------- %
%%

function [Time,MExcitation,MActivation,RActivation,TForcetilde,TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=SolveMuscleRedundancy_lMtildeState_actdyn_GPOPS(model_path,IK_path,ID_path,time,OutPath,Misc)

%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART I: INPUTS FOR OPTIMAL CONTROL PROBLEM ---------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Check for optional input arguments, see manual for details------------- %

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

% Round time window to 2 decimals
time=round(time,2);
if time(1)==time(2)
   warning('Time window should be at least 0.01s'); 
end

% ------------------------------------------------------------------------%
% Compute ID -------------------------------------------------------------%
if isempty(ID_path) || ~exist(ID_path,'file')
    disp('ID path was not specified or the file does not exist, computation ID started');
    if ~isfield(Misc,'Loads_path') || isempty(Misc.Loads_path) || ~exist(Misc.Loads_path,'file')
        error('External loads file was not specified or does not exist, please add the path to the external loads file: Misc.Loads_path');
    else
        %check the output path for the ID results
        if isfield(Misc,'ID_ResultsPath')
            [idpath,~]=fileparts(Misc.ID_ResultsPath);
            if ~isdir(idpath); mkdir(idpath); end
        else 
            % save results in the directory of the external loads
            [Lpath,name,~]=fileparts(Misc.Loads_path);
            Misc.ID_ResultsPath=fullfile(Lpath,name);
        end
        [ID_outPath,ID_outName,ext]=fileparts(Misc.ID_ResultsPath);
        output_settings=fullfile(ID_outPath,[ID_outName '_settings.xml']);
        Opensim_ID(model_path,[time(1)-0.1 time(2)+0.1],Misc.Loads_path,IK_path,ID_outPath,[ID_outName ext],output_settings);
        ID_path=Misc.ID_ResultsPath;
    end    
end

% ----------------------------------------------------------------------- %
% Muscle analysis ------------------------------------------------------- %

Misc.time=time;
MuscleAnalysisPath=fullfile(OutPath,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end
disp('MuscleAnalysis Running .....');
OpenSim_Muscle_Analysis(IK_path,model_path,MuscleAnalysisPath,[time(1) time(end)])
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

% The solution of the static optimization is used as initial guess for the
% dynamic optimization
% Extract the muscle-tendon properties
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=ReadMuscleParameters(model_path,DatStore.MuscleNames);
% Static optimization using IPOPT solver
DatStore = SolveStaticOptimization_IPOPT_GPOPS(DatStore);

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
auxdata.scaling.vMtilde = 10;

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
a_min = 0; a_max = 1;                   % bounds on muscle activation
vA_min = -1/100; vA_max = 1/100;        % bounds on derivative of muscle activation (scaled)
vMtilde_min = -1; vMtilde_max = 1;      % bounds on normalized muscle fiber velocity
lMtilde_min = 0.2; lMtilde_max = 1.8;   % bounds on normalized muscle fiber length

% Time bounds
t0 = DatStore.time(1); tf = DatStore.time(end);
bounds.phase.initialtime.lower = t0; bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; bounds.phase.finaltime.upper = tf;
% Controls bounds
vAmin = vA_min./auxdata.tauDeact; vAmax = vA_max./auxdata.tauAct;
vMtildemin = vMtilde_min*ones(1,auxdata.NMuscles); vMtildemax = vMtilde_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof); aTmax = 1*ones(1,auxdata.Ndof);
bounds.phase.control.lower = [vAmin aTmin vMtildemin]; bounds.phase.control.upper = [vAmax aTmax vMtildemax];
% States bunds
actMin = a_min*ones(1,auxdata.NMuscles); actMax = a_max*ones(1,auxdata.NMuscles);
lMtildemin = lMtilde_min*ones(1,auxdata.NMuscles); lMtildemax = lMtilde_max*ones(1,auxdata.NMuscles);
bounds.phase.initialstate.lower = [actMin, lMtildemin]; bounds.phase.initialstate.upper = [actMax, lMtildemax];
bounds.phase.state.lower = [actMin, lMtildemin]; bounds.phase.state.upper = [actMax, lMtildemax];
bounds.phase.finalstate.lower = [actMin, lMtildemin]; bounds.phase.finalstate.upper = [actMax, lMtildemax];
% Integral bounds
bounds.phase.integral.lower = 0; bounds.phase.integral.upper = 10000*(tf-t0);

% Path constraints
HillEquil = zeros(1, auxdata.NMuscles);
ID_bounds = zeros(1, auxdata.Ndof);
act1_lower = zeros(1, auxdata.NMuscles);
act1_upper = inf * ones(1, auxdata.NMuscles);
act2_lower = -inf * ones(1, auxdata.NMuscles);
act2_upper =  ones(1, auxdata.NMuscles)./auxdata.tauAct;
bounds.phase.path.lower = [ID_bounds,HillEquil,act1_lower,act2_lower]; bounds.phase.path.upper = [ID_bounds,HillEquil,act1_upper, act2_upper];

% Eventgroup
% Impose mild periodicity
pera_lower = -1 * ones(1, auxdata.NMuscles); pera_upper = 1 * ones(1, auxdata.NMuscles);
perlMtilde_lower = -1*ones(1,auxdata.NMuscles); perlMtilde_upper = 1*ones(1,auxdata.NMuscles);
bounds.eventgroup.lower = [pera_lower perlMtilde_lower]; bounds.eventgroup.upper = [pera_upper perlMtilde_upper];

% Initial guess static optimization 
time_opt = DatStore.time;
SoActInterp = interp1(DatStore.time,DatStore.SoAct,time_opt);
SoRActInterp = interp1(DatStore.time,DatStore.SoRAct,time_opt);
SoForceInterp = interp1(DatStore.time,DatStore.SoForce.*DatStore.cos_alpha./DatStore.Fiso,time_opt);
[~,lMtildeInterp ] = FiberLength_Ftilde(SoForceInterp,DatStore.params,DatStore.LMT,Misc.Atendon,Misc.shift);
vMtildeinterp = zeros(size(lMtildeInterp));
for m = 1:auxdata.NMuscles
lMtildeSpline = spline(time_opt,lMtildeInterp(:,m));
[~,vMtildeinterp_norm,~] = SplineEval_ppuval(lMtildeSpline,time_opt,1);
vMtildeinterp(:,m) = vMtildeinterp_norm/auxdata.scaling.vMtilde;
end

% Initial guess
N = length(DatStore.time);
guess.phase.time = DatStore.time;
% Based on SO result
guess.phase.control = [zeros(N,auxdata.NMuscles) SoRActInterp./150 vMtildeinterp];
guess.phase.state =  [SoActInterp lMtildeInterp];
% Random
% guess.phase.control = [zeros(N,auxdata.NMuscles) zeros(N,auxdata.Ndof) 0.01*ones(N,auxdata.NMuscles)];
% guess.phase.state =  [0.2*ones(N,auxdata.NMuscles) ones(N,auxdata.NMuscles)];
guess.phase.integral = 0;

% Spline structures
for dof = 1:auxdata.Ndof
    for m = 1:auxdata.NMuscles       
        auxdata.JointMASpline(dof).Muscle(m) = spline(DatStore.time,auxdata.MA(dof).Joint(:,m));       
    end
    auxdata.JointIDSpline(dof) = spline(DatStore.time,DatStore.T_exp(:,dof));
end

for m = 1:auxdata.NMuscles
    auxdata.LMTSpline(m) = spline(DatStore.time,DatStore.LMT(:,m));
end

% GPOPS setup
setup.name = 'DynamicOptimization_lMtildeState_actdyn_GPOPS_';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'mumps';
setup.derivatives.derivativelevel = 'second';
setup.nlp.ipoptoptions.tolerance = 1e-6;
setup.nlp.ipoptoptions.maxiterations = 10000;
setup.derivatives.supplier = 'adigator';
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 0;
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
NMeshIntervals = round((tf-t0)*Misc.Mesh_Frequency);
setup.mesh.phase.colpoints = 3*ones(1,NMeshIntervals);
setup.mesh.phase.fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous = @musdynContinous_lMtildeState_vA;
setup.functions.endpoint = @musdynEndpoint_lMtildeState;
    
% ADiGator setup
persistent splinestruct
input.auxdata = auxdata;
tdummy = guess.phase.time;
splinestruct = SplineInputData(tdummy,input);
splinenames = fieldnames(splinestruct);
for Scount = 1:length(splinenames)
  secdim = size(splinestruct.(splinenames{Scount}),2);
  splinestructad.(splinenames{Scount}) = adigatorCreateAuxInput([Inf,secdim]);
  splinestruct.(splinenames{Scount}) = zeros(0,secdim);
end
setup.auxdata.splinestruct = splinestructad;
adigatorGenFiles4gpops2(setup)

setup.functions.continuous = @Wrap4musdynContinous_lMtildeState_vA;
setup.adigatorgrd.continuous = @musdynContinous_lMtildeState_vAGrdWrap;
setup.adigatorgrd.endpoint   = @musdynEndpoint_lMtildeStateADiGatorGrd;
setup.adigatorhes.continuous = @musdynContinous_lMtildeState_vAHesWrap;
setup.adigatorhes.endpoint   = @musdynEndpoint_lMtildeStateADiGatorHes;

%% ---------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
% PART III: SOLVE OPTIMAL CONTROL PROBLEM ------------------------------- %
% ----------------------------------------------------------------------- %
% ----------------------------------------------------------------------- %
output = gpops2(setup);

% Delete output files from ADiGator
delete musdynEndpoint_lMtildeStateADiGatorGrd.mat
delete musdynEndpoint_lMtildeStateADiGatorGrd.m
delete musdynEndpoint_lMtildeStateADiGatorHes.mat
delete musdynEndpoint_lMtildeStateADiGatorHes.m
delete musdynContinous_lMtildeState_vAADiGatorHes.mat
delete musdynContinous_lMtildeState_vAADiGatorHes.m
delete musdynContinous_lMtildeState_vAADiGatorGrd.mat
delete musdynContinous_lMtildeState_vAADiGatorGrd.m

res=output.result.solution.phase(1);
Time=res.time;
MActivation=res.state(:,1:auxdata.NMuscles);
lMtilde=res.state(:,auxdata.NMuscles+1:auxdata.NMuscles*2);
lM=lMtilde.*(ones(length(Time),1)*DatStore.lOpt);
vA=100*res.control(:,1:auxdata.NMuscles);
MExcitation=computeExcitationRaasch(MActivation,vA, auxdata.tauDeact, auxdata.tauAct);
RActivation=res.control(:,auxdata.NMuscles+1:auxdata.NMuscles+auxdata.Ndof)*150;
MuscleNames=DatStore.MuscleNames;
OptInfo=output;

% Tendon force from lMtilde
% Interpolation lMT
lMTinterp = interp1(DatStore.time,DatStore.LMT,Time);
[TForcetilde,TForce] = TendonForce_lMtilde(lMtilde,auxdata.params,lMTinterp,auxdata.Atendon,auxdata.shift);
end
