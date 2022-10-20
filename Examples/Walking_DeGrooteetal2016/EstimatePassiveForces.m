%% Default example solve muscle redundancy
% (as in DeGroote2016)

clear all; clc; close all;
%% Input information

% select datafolder
ExamplePath = pwd;
DataPath = [pwd '\DataPassiveMotion'];

% Add here the paths of IK, ID , US and EMG data trials you want to work with
Misc.model_path  = fullfile(DataPath,'subject1.osim');

Misc.OutPath     = fullfile(ExamplePath,'PassiveEst');                    % folder to store results

Misc.IKfile = {fullfile(DataPath,'trialAnkle_knee-15hip0_knee-60hip0.mot'),...
    fullfile(DataPath,'trialHip_ankle0knee0_ankle0knee-15_ankle0knee-60.mot'),...
    fullfile(DataPath,'trialKnee_ankle0hip-20_ankle0hip15_ankle15hip-20.mot')};

Misc.IDfile = {fullfile(DataPath,'trialAnkle_knee-15hip0_knee-60hip0.sto'),...
    fullfile(DataPath,'trialHip_ankle0knee0_ankle0knee-15_ankle0knee-60.sto'),...
    fullfile(DataPath,'trialKnee_ankle0hip-20_ankle0hip15_ankle15hip-20.sto')};


% name output
Misc.OutName = 'ParmEst';

% input dofs
Misc.DofNames_Input={{'ankle_angle_r'};{'hip_flexion_r'};{'knee_angle_r'}};    % select the DOFs you want to include in the optimization

% adapt the stiffness of the achilles tendon (optional input argument)
Misc.Set_kT_ByName = {'soleus_r',20;
    'med_gas_r',20;
    'lat_gas_r',20};

time = [0 10;...
    0 10;...
    0 10];

% run muscle analysis
Misc.RunAnalysis = false;


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

for trial = 1:Misc.nTrials
    nMuscles = length(DatStore(trial).MuscleNames);
    Nfr = length(DatStore(trial).time);
    lTs = Misc.lTs(Misc.idx_allMuscleList{trial})';
    IG(trial).lM_projected = zeros(nMuscles,Nfr);
    IG(trial).lMtilde = zeros(nMuscles,Nfr);
    for k=1:Nfr
        lMT = DatStore(trial).LMT(k,:)';
        lMGuess = lMT-lTs;
        lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
        alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
        w = lMo.*sin(alphao);
        IG(trial).lM_projected(:,k) = sqrt((lMGuess.^2 - w.^2));
        IG(trial).lMtilde(:,k) = lMGuess./lMo;
    end

    figure();
    plot(IG(trial).lMtilde');
    xlabel('frames');
    ylabel('rigid tendon estimate fiber length');
end

%% formulate optimization problem
import casadi.*
opti    = casadi.Opti();    % create opti structure
J = 0;
for trial=1:Misc.nTrials
    % get number of muscles and frames
    nMuscles = length(DatStore(trial).MuscleNames);
    Nfr = length(DatStore(trial).time);

    % optimization variables
    TR(trial).kpe = opti.variable(nMuscles,1);
    TR(trial).ksf = opti.variable(nMuscles,1);

    % bounds
    opti.subject_to( 2 < TR(trial).kpe <6);
    opti.subject_to( 0.4 < TR(trial).ksf < 3 );

    % set initial gues
    opti.set_initial(TR(trial).kpe, 4);
    opti.set_initial(TR(trial).ksf, 1);

    % projected fiber length as optimization variable
    TR(trial).lM_projected = opti.variable(nMuscles,Nfr);
    opti.subject_to(1e-4 < TR(trial).lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
    opti.set_initial(TR(trial).lM_projected,IG(trial).lM_projected);

    % create reserve actuators
    TR(trial).aTk = opti.variable(DatStore(trial).nDOF,Nfr);
    opti.subject_to(-1< TR(trial).aTk < 1)

    % lMtilde as optimization variable
    TR(trial).lMtilde = opti.variable(nMuscles,Nfr);
    opti.subject_to(0.05 < TR(trial).lMtilde < 2);
    opti.set_initial(TR(trial).lMtilde,IG(trial).lMtilde);

    % select muscle properties
    MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
    MuscProperties.kT = Misc.kT(:,Misc.idx_allMuscleList{trial}')';
    MuscProperties.shift = Misc.shift(:,Misc.idx_allMuscleList{trial}')';

    % constraint on projected fiber length
    lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
    alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
    lM = TR(trial).lMtilde.*lMo;
    w = lMo.*sin(alphao);
    opti.subject_to(lM.^2 - w.^2 == TR(trial).lM_projected.^2);

    % find passive
    TR(trial).TsimVect = MX(Nfr,DatStore(trial).nDOF);
    for k=1:Nfr
        %         disp(['nVar ' num2str(opti.nx) '  nConstr ' num2str(opti.ng)]);
        % solve muscle equilibrium
        [err, FTk] = ForceEquilibrium_lMtildeState_optPassive(0.05,TR(trial).lMtilde(:,k),...
            0,TR(trial).lM_projected(:,k),DatStore(trial).LMT(k,:)',MuscProperties.params,...
            MuscProperties.kT,MuscProperties.shift,TR(trial).kpe,TR(trial).ksf);

        % Moment constraints
        for dof = 1:DatStore(trial).nDOF
            T_exp = DatStore(trial).T_exp(k,dof); % + 1 due to the time col
            T_sim = FTk'*squeeze(DatStore(trial).dM(k,dof,:)) + Misc.Topt*TR(trial).aTk(dof,k);
            opti.subject_to(T_exp - T_sim == 0);
            TR(trial).TsimVect(k,dof) = FTk'*squeeze(DatStore(trial).dM(k,dof,:));
        end
        % impose equilibrium
        opti.subject_to(err == 0);
    end
    J = J + sumsqr(TR(trial).aTk) + 0.01.*sumsqr(TR(trial).kpe/100) + 0.01.*sumsqr(TR(trial).ksf/100);

    % display contstrains
    %     disp(['nVar ' num2str(opti.nx) '  nConstr ' num2str(opti.ng)]);
end

% J = sumsqr(aTk) + 0.01.*sumsqr(kpe/100) + 0.01.*sumsqr(ksf/100);
opti.minimize(J);
opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);
optiDef = opti.copy(); % copy for solution with default params

% impose contraint that values should be equal for each muscle in the
% ankle, knee and hip trials (this is important for bi-articular muscles)
% find all muscles of trial 1 that appear in trial 2 or trial 3
Trial1Copies = nan(10,3);
ct = 1;
for iM = 1:length(DatStore(1).MuscleNames)
    iTrial2 = find(strcmp(DatStore(2).MuscleNames,DatStore(1).MuscleNames{iM}));
    iTrial3 = find(strcmp(DatStore(3).MuscleNames,DatStore(1).MuscleNames{iM}));
    if ~isempty(iTrial2)
        Trial1Copies(ct,1) = iM;
        Trial1Copies(ct,2) = iTrial2;
        Trial1Copies(ct,3) = 2;
        ct = ct+1;
    end
    if ~isempty(iTrial3)
        Trial1Copies(ct,1) = iM;
        Trial1Copies(ct,2) = iTrial3;
        Trial1Copies(ct,3) = 3;
        ct = ct+1;
    end
end
Trial1Copies(ct:end,:) = [];

Trial2Copies = nan(10,3);
ct = 1;
for iM = 1:length(DatStore(2).MuscleNames)
    iTrial3 = find(strcmp(DatStore(3).MuscleNames,DatStore(2).MuscleNames{iM}));   
    if ~isempty(iTrial3)
        Trial2Copies(ct,1) = iM;
        Trial2Copies(ct,2) = iTrial3;
        Trial2Copies(ct,3) = 3;
        ct = ct+1;
    end
end
Trial2Copies(ct:end,:) = [];

% set equality constraints
for i=1:length(Trial1Copies(:,1))
    % equality constraint
    iM_tr1 = Trial1Copies(i,1);
    iM_trx = Trial1Copies(i,2);
    trx = Trial1Copies(i,3);
    opti.subject_to(TR(1).kpe(iM_tr1,1) == TR(trx).kpe(iM_trx));
    opti.subject_to(TR(1).ksf(iM_tr1,1) == TR(trx).ksf(iM_trx));
end
for i=1:length(Trial2Copies(:,1))
    % equality constraint
    iM_tr2 = Trial2Copies(i,1);
    iM_trx = Trial2Copies(i,2);
    trx = Trial2Copies(i,3);
    opti.subject_to(TR(2).kpe(iM_tr2,1) == TR(trx).kpe(iM_trx));
    opti.subject_to(TR(2).ksf(iM_tr2,1) == TR(trx).ksf(iM_trx));
end

Sol = opti.solve();
disp(['Number of bi-articular muscles ' num2str(length(Trial2Copies(:,1)) + length(Trial1Copies(:,1)))])
% wOut = solve_NLPSOL(opti,SolverSetup.optionssol)

% solve without default values
for trial=1:Misc.nTrials
    % equality constraints for parameters of passive force length
    optiDef.subject_to( TR(trial).kpe == 4);
    optiDef.subject_to( TR(trial).ksf == 1 );
    nMuscles = length(DatStore(trial).MuscleNames);
    Nfr = length(DatStore(trial).time);
    % Loop over frames just to compute the torque generqted by the muscles
    % (this is actually postprocessing to create the matrix TsimVect)
    MuscProperties.params = Misc.params(:,Misc.idx_allMuscleList{trial})';
    MuscProperties.kT = Misc.kT(:,Misc.idx_allMuscleList{trial}')';
    MuscProperties.shift = Misc.shift(:,Misc.idx_allMuscleList{trial}')';
    TRDef(trial).TsimVect = MX(Nfr,DatStore(trial).nDOF);
    for k=1:Nfr
        [err, FTk] = ForceEquilibrium_lMtildeState_optPassive(0.05,TR(trial).lMtilde(:,k),...
            0,TR(trial).lM_projected(:,k),DatStore(trial).LMT(k,:)',MuscProperties.params,...
            MuscProperties.kT,MuscProperties.shift,TR(trial).kpe,TR(trial).ksf);
       
        for dof = 1:DatStore(trial).nDOF
            T_sim = FTk'*squeeze(DatStore(trial).dM(k,dof,:)) + Misc.Topt*TR(trial).aTk(dof,k);
            TRDef(trial).TsimVect(k,dof) = FTk'*squeeze(DatStore(trial).dM(k,dof,:));
        end
    end
end
SolDefault = optiDef.solve();


%% Get solutions

for trial=1:Misc.nTrials
    % solution with parameter optimization
    Out(trial).kpe = Sol.value(TR(trial).kpe);
    Out(trial).ksf = Sol.value(TR(trial).ksf);
    Out(trial).aTk = Sol.value(TR(trial).aTk).*Misc.Topt;
    Out(trial).Tsim = Sol.value(TR(trial).TsimVect);
    % solution with default parameters
    OutDef(trial).kpe = SolDefault.value(TR(trial).kpe);
    OutDef(trial).ksf = SolDefault.value(TR(trial).ksf);
    OutDef(trial).aTk = SolDefault.value(TR(trial).aTk).*Misc.Topt;
    OutDef(trial).Tsim = SolDefault.value(TRDef(trial).TsimVect);
end

%% Plot outcome
% plot torque on each frame for default model, estimated model and ID
JointLeg = {'ankle','hip','knee'};
figure();
for trial=1:Misc.nTrials
    subplot(1,3,trial)
    plot(Out(trial).Tsim,'r'); hold on;
    plot(OutDef(trial).Tsim,'b'); hold on;
    plot(DatStore(trial).T_exp,'--k');
    xlabel('frame');
    ylabel('Moment [Nm]')
    title(JointLeg{trial})
end
legend('estimated','defaultsim','exp');

figure();
mk = 4;
for trial=1:Misc.nTrials
    subplot(1,3,trial)
    plot(DatStore(trial).q_exp,Out(trial).Tsim,'or','MarkerFaceColor',[1 0 0],'MarkerSize',mk); hold on;
    plot(DatStore(trial).q_exp,OutDef(trial).Tsim,'ob','MarkerFaceColor',[0 0 1],'MarkerSize',mk); hold on;
    plot(DatStore(trial).q_exp,DatStore(trial).T_exp,'ok','MarkerFaceColor',[0 0 0],'MarkerSize',mk);
    xlabel('angle (deg)');
    ylabel('Moment [Nm]')
    title(JointLeg{trial})
end
legend('estimated','defaultsim','exp');
