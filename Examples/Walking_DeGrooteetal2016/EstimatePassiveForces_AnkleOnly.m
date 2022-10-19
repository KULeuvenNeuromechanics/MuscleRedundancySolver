%% Default example solve muscle redundancy
% (as in DeGroote2016)

clear all; clc; close all;
%% Input information

% select datafolder
ExamplePath = pwd;
DataPath = [pwd '\WalkingData'];

% Add here the paths of IK, ID , US and EMG data trials you want to work with
Misc.model_path  = fullfile(DataPath,'subject1.osim');
Misc.OutPath     = fullfile(ExamplePath,'PassiveEst');                    % folder to store results
Misc.IKfile = {fullfile(ExamplePath,'silder_IK.mot')};
Misc.IDfile = {fullfile(ExamplePath,'silder_ID.sto')};
% IK_file = fullfile(ExamplePath,'silder_IK.mot');
% ID_file = fullfile(ExamplePath,'silder_ID.sto');

% name output
Misc.OutName = 'Walking3_';

% input dofs
Misc.DofNames_Input={'ankle_angle_r'};    % select the DOFs you want to include in the optimization

% adapt the stiffness of the achilles tendon (optional input argument)
Misc.Set_kT_ByName = {'soleus_r',20;
    'med_gas_r',20;
    'lat_gas_r',20};

time = [0 10];

%% Remove the extreme values from the IK file

IK = ReadMotFile(Misc.IKfile{1});
ID = ReadMotFile(Misc.IDfile{1});

% select some specific frames
i0 = find(ID.data(:,2) == 0);
iFrames = 1:length(IK.data(:,1));
iFrames(i0) = [];

% I think we have to change the direction of the knee angle
IK_header = {'time'	'ankle_angle_r','knee_angle_r','hip_flexion_r'};
IK.data(:,3) = -IK.data(:,3);
IK.data = IK.data(iFrames,:);
IKOutNew = fullfile(ExamplePath,'silder_IK_adapt.mot');
generateMotFile(IK.data, IK_header, IKOutNew);

ID_header = {'time','ankle_angle_r_moment','knee_angle_r_moment','hip_flexion_r_moment'};
ID.data(:,3) = -ID.data(:,3);
ID.data = ID.data(iFrames,:);
IDOutNew = fullfile(ExamplePath,'silder_ID_adapt.sto');
generateMotFile(ID.data, ID_header, IDOutNew);

Misc.IKfile = {IKOutNew};
Misc.IDfile = {IDOutNew};

time = [ID.data(1,1) ID.data(end,1)];

%% Run muscle tendon estimator:

% update default settings
Misc = DefaultSettings(Misc);

% read the muscle properties
[Misc] = getMuscleProperties(Misc.model_path,Misc);

% number of motion trials
Misc.nTrials = length(Misc.IKfile);

%check if we have to adapt the start and end time so that it corresponds to
%time frames in the IK solution
[time] = Check_TimeIndices(Misc,time);
Misc.time=time;

%% Perform muscle analysis for all trials
DatStore = struct;
MuscleAnalysisPath=fullfile(Misc.OutPath,'MuscleAnalysis');
Misc.MuscleAnalysisPath=MuscleAnalysisPath;
if ~exist(MuscleAnalysisPath,'dir')
    mkdir(MuscleAnalysisPath);
end
for i = 1:Misc.nTrials
    % select the IK and ID file
    IK_path_trial = Misc.IKfile{i};
    % Run muscle analysis
    if Misc.RunAnalysis
        disp('MuscleAnalysis Running .....');
        OpenSim_Muscle_Analysis(IK_path_trial,Misc.model_path,MuscleAnalysisPath,[time(i,1) time(i,end)],Misc.DofNames_Input{i,1})
        disp('MuscleAnalysis Finished');
    end
end

%% Extract muscle information
% Get number of degrees of freedom (dofs), muscle-tendon lengths, moment
% arms, stiffness and shift for the selected muscles.
for trial = 1:Misc.nTrials
    [~,Misc.MAtrialName{trial},~]=fileparts(Misc.IKfile{trial});
end

% select muscles with a moment arm for the input dofs
Misc = getMuscles4DOFS(Misc);

% get IK, ID, muscle analysis data
[Misc,DatStore] = getMuscleInfo(Misc,DatStore);

% display warnings in muscle selection
[Misc] = Warnings_MuscleNames(DatStore,Misc);

% get indexes of the muscles for which optimal fiber length, tendon stiffness are estimated
[DatStore,Misc] = GetIndices(DatStore,Misc);

% get the EMG and ultrasound information
[Misc,DatStore] = GetEMGInfo(Misc,DatStore);
[Misc,DatStore] = GetUSInfo(Misc,DatStore);

% get the number of muscles in a vector
NMuscles = zeros(Misc.nTrials,1);
for trial = 1:Misc.nTrials
    NMuscles(trial) = DatStore(trial).nMuscles;
end

%% set solver

% Create an NLP solver
SolverSetup.nlp.solver = 'ipopt';
% Set derivativelevel to 'first' for approximating the Hessian
SolverSetup.derivatives.derivativelevel = 'second';
% if strcmp(SolverSetup.derivatives.derivativelevel, 'first')
% SolverSetup.optionssol.ipopt.hessian_approximation = 'limited-memory';
% end
SolverSetup.optionssol.ipopt.nlp_scaling_method = 'gradient-based';
SolverSetup.optionssol.ipopt.linear_solver = 'mumps';
SolverSetup.optionssol.ipopt.tol = 1e-6;
SolverSetup.optionssol.ipopt.max_iter = 10000;

%% initial guess
nMuscles = length(DatStore.MuscleNames);
trial = 1;
Nfr = length(DatStore.time);
lTs = Misc.lTs(Misc.idx_allMuscleList{trial})';
IG.lM_projected = zeros(nMuscles,Nfr);
IG.lMtilde = zeros(nMuscles,Nfr);
for k=1:Nfr
    lMT = DatStore.LMT(k,:)';
    lMGuess = lMT-lTs;
    lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
    alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
    w = lMo.*sin(alphao);
    IG.lM_projected(:,k) = sqrt((lMGuess.^2 - w.^2));
    IG.lMtilde(:,k) = lMGuess./lMo;
end

figure();
plot(IG.lMtilde');
xlabel('frames');
ylabel('rigid tendon estimate fiber length');


%% computer fiber lengths for given passive F/l properties

% thi sjust solves for equilibrium between tendon and muscle force

%
import casadi.*
opti    = casadi.Opti();    % create opti structure

% projected fiber length as optimization variable
Nfr = length(DatStore.time);
lM_projected = opti.variable(nMuscles,Nfr);
opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
opti.set_initial(lM_projected,IG.lM_projected);

% lMtilde as optimization variable
lMtilde = opti.variable(nMuscles,Nfr);
opti.subject_to(0.05 < lMtilde < 2);
opti.set_initial(lMtilde,1);

% trial number 1
trial = 1;

% select muscle properties
MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
MuscProperties.kT = Misc.kT(:,Misc.idx_allMuscleList{trial}')';
MuscProperties.shift = Misc.shift(:,Misc.idx_allMuscleList{trial}')';

% constraint on projected fiber length
lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
opti.subject_to(lM.^2 - w.^2 == lM_projected.^2);

% impose hill equilibrium
TsimVect = MX(Nfr,3);
for k=1:Nfr
    % solve muscle equilibrium
    [err, FTk] = ForceEquilibrium_lMtildeState_optPassive(0,lMtilde(:,k),0,lM_projected(:,k),...
        DatStore.LMT(k,:)',MuscProperties.params,MuscProperties.kT,MuscProperties.shift,4,1);
    % impose equilibrium
    opti.subject_to(err == 0);
    % comput moment generated by muscles
    for dof = 1:DatStore(trial).nDOF
        T_sim = FTk'*squeeze(DatStore(trial).dM(k,dof,:));
        TsimVect(k,dof) = T_sim;
    end
end

opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);
SolGuess = opti.solve();

% store values for initial guess
IGSol.lMProjected = SolGuess.value(lM_projected);
IGSol.lMtilde = SolGuess.value(lMtilde);
IGSol.Tsim = SolGuess.value(TsimVect);

%% formulate optimization problem

opti    = casadi.Opti();    % create opti structure

% optimization variables
kpe = opti.variable(nMuscles,1);
ksf = opti.variable(nMuscles,1);

% bounds
opti.subject_to( 1 < kpe <10);
opti.subject_to( 0.1 < ksf < 3 );

% set initial gues
opti.set_initial(kpe, 4);
opti.set_initial(ksf, 1);

% projected fiber length as optimization variable
Nfr = length(DatStore.time);
lM_projected = opti.variable(nMuscles,Nfr);
opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
opti.set_initial(lM_projected,IGSol.lMProjected);

% create reserve actuators
aTk = opti.variable(3,Nfr);
opti.subject_to(-1< aTk < 1)

% lMtilde as optimization variable
lMtilde = opti.variable(nMuscles,Nfr);
opti.subject_to(0.05 < lMtilde < 2);
opti.set_initial(lMtilde,IGSol.lMtilde);

% trial number 1
trial = 1;

% select muscle properties
MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
MuscProperties.kT = Misc.kT(:,Misc.idx_allMuscleList{trial}')';
MuscProperties.shift = Misc.shift(:,Misc.idx_allMuscleList{trial}')';

% constraint on projected fiber length
lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
opti.subject_to(lM.^2 - w.^2 == lM_projected.^2);

% find passive
TsimVect = MX(Nfr,3);
for k=1:Nfr
    % solve muscle equilibrium
        [err, FTk] = ForceEquilibrium_lMtildeState_optPassive(0.05,lMtilde(:,k),0,lM_projected(:,k),...
            DatStore.LMT(k,:)',MuscProperties.params,MuscProperties.kT,MuscProperties.shift,kpe,ksf);
    % compute joint moment
    % Add path constraints
    % Moment constraints
    for dof = 1:DatStore(trial).nDOF
        T_exp = ID.data(k,dof+1); % + 1 due to the time col
        T_sim = FTk'*squeeze(DatStore(trial).dM(k,dof,:)) + Misc.Topt*aTk(dof,k);
        opti.subject_to(T_exp - T_sim == 0);
        TsimVect(k,dof) = FTk'*squeeze(DatStore(trial).dM(k,dof,:));
    end
    % impose equilibrium
    opti.subject_to(err == 0);
end

% J = sumsqr(aTk) + 0.01.*sumsqr(kpe/100) + 0.01.*sumsqr(ksf/100);
J = sumsqr(aTk) + 0.01.*sumsqr(kpe/100) + 0.01.*sumsqr(ksf/100);
opti.minimize(J);
opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);
Sol = opti.solve();
% wOut = solve_NLPSOL(opti,SolverSetup.optionssol)

%% plot solution

Out.kpe = Sol.value(kpe);
Out.ksf = Sol.value(ksf);
Out.aTk = Sol.value(aTk).*Misc.Topt;
Out.Tsim = Sol.value(TsimVect);

% plot torque error on each frame
figure();
for i=1:DatStore.nDOF
    subplot(1,DatStore.nDOF,i)
    plot(Out.aTk(i,:)); hold on;
    xlabel('frames');
    ylabel('torque error');
    title(DatStore.DOFNames{i});
end

% plot torque as a function of the joint angle
figure();
for i=1:DatStore.nDOF
    subplot(1,DatStore.nDOF,i)
    plot(IK.data(:,i+1),Out.Tsim(:,i),'or'); hold on;
    plot(IK.data(:,i+1),IGSol.Tsim(:,i),'ob'); hold on;
    plot(IK.data(:,i+1),ID.data(:,i+1),'ok'); hold on;
    xlabel('joint angle');
    ylabel('torque');
    title(DatStore.DOFNames{i});
end
legend('estimated','defaultsim','exp');


% plot torque on each frame for default model, estimated model and ID
figure();
for i=1:DatStore.nDOF
    subplot(1,3,DatStore.nDOF)
    plot(Out.Tsim(:,i),'r'); hold on;
    plot(IGSol.Tsim(:,i),'b'); hold on;
    plot(ID.data(:,i+1),'--k');
end
legend('estimated','defaultsim','exp');

%% compute error
JSol = sumsqr((Out.Tsim-ID.data(:,2:end))/Misc.Topt);
JSol_Guess = sumsqr((IGSol.Tsim-ID.data(:,2:end))/Misc.Topt);

disp('Objective value ');
disp([' Parameter optimization ', num2str(JSol)]);
disp([' default values ', num2str(JSol_Guess)]);


%% plot

figure()
for i=1:3
    subplot(1,3,i)
    plot(IK.data(:,i+1),ID.data(:,i+1),'ok'); hold on;
%     xlabel('ankle angle');
%     ylabel('ankle moment (ID)');
end