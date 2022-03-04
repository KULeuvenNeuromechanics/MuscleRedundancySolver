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
    for i = 19:Misc.nTrials
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
%         [Results] = runMRS_allTrialsTogether(Misc,DatStore,Mesh,output,optionssol,Results,NMuscles);

        % Problem bounds
        e_min = 0; e_max = 1;                   % bounds on muscle excitation
        a_min = 0; a_max = 1;                   % bounds on muscle activation
        vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
        lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

        % CasADi setup
        import casadi.*
        opti    = casadi.Opti();    % create opti structure

        for trial = Misc.trials_sel
            % set intial guess based on static opt data
            SoActGuess{trial} = zeros(NMuscles(trial),Mesh(trial).N+1);
            SoExcGuess{trial} = zeros(NMuscles(trial),Mesh(trial).N);
            lMtildeGuess{trial} = zeros(NMuscles(trial),Mesh(trial).N+1);
            vMtildeGuess{trial} = zeros(NMuscles(trial),Mesh(trial).N);
            SoRActGuess{trial} = zeros(DatStore(trial).nDOF,Mesh(trial).N);

            SoActGuess{trial} = DatStore(trial).SoActInterp';
            SoExcGuess{trial} = DatStore(trial).SoActInterp(1:end-1,:)';
            lMtildeGuess{trial} = DatStore(trial).lMtildeInterp';
            vMtildeGuess{trial} = DatStore(trial).vMtildeinterp(1:end-1,:)';
            SoRActGuess{trial} = DatStore(trial).SoRActInterp(1:end-1,:)';
        end

        J = 0; 

        % Loop over trials --> one simulation for each trial
        ct = 0;
        for trial = Misc.trials_sel
        %     for trial = 5:6
            ct = ct + 1;
            % States
            %   - muscle activation
            a{ct} = opti.variable(NMuscles(trial),Mesh(trial).N+1);      % Variable at mesh points
            opti.subject_to(a_min < a{ct} < a_max);             % Bounds
            opti.set_initial(a{ct},SoActGuess{trial});                 % Initial guess (static optimization)
            %   - Muscle fiber lengths
            lMtilde{ct} = opti.variable(NMuscles(trial),Mesh(trial).N+1);
            opti.subject_to(lMtilde_min < lMtilde{ct} < lMtilde_max);
            opti.set_initial(lMtilde{ct},lMtildeGuess{trial});
            %   - Controls
            e{ct} = opti.variable(NMuscles(trial),Mesh(trial).N);
            opti.subject_to(e_min < e{ct} < e_max);
            opti.set_initial(e{ct}, SoExcGuess{trial});
            %   - Reserve actuators
            aT{ct} = opti.variable(DatStore(trial).nDOF,Mesh(trial).N);
            opti.subject_to(-1 < aT{ct} <1);
            %   - Time derivative of muscle-tendon forces (states)
            vMtilde{ct} = opti.variable(NMuscles(trial),Mesh(trial).N);
            opti.subject_to(vMtilde_min < vMtilde{ct} < vMtilde_max);
            opti.set_initial(vMtilde{ct},vMtildeGuess{trial});    
            %   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
            lM_projected{ct} = opti.variable(NMuscles(trial),Mesh(trial).N+1);
            opti.subject_to(1e-4 < lM_projected{ct}(:)); % We impose that projected muscle fiber length has strict positive length
            % Initial guess for this variable is retrieved from lMtilde guess
            % and geometric relationship between pennation angle, muscle length
            % and width
            lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
            alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
            lMGuess = lMtildeGuess{trial}.*lMo;
            w = lMo.*sin(alphao);
            lM_projectedGuess{trial} = sqrt((lMGuess.^2 - w.^2));
            opti.set_initial(lM_projected{ct},lM_projectedGuess{trial});

            % constraint on projected fiber length
            w = lMo.*sin(alphao);  
            lM = lMtilde{ct}.*lMo;
            opti.subject_to(lM.^2 - w.^2 == lM_projected{ct}.^2);

            % Discretization
            h = Mesh(trial).step;

            % Loop over mesh points formulating NLP
            for k=1:Mesh(trial).N
                % Variables within current mesh interval
                ak = a{ct}(:,k); lMtildek = lMtilde{ct}(:,k);
                vMtildek = vMtilde{ct}(:,k); aTk = aT{ct}(:,k);
                ek = e{ct}(:,k); lM_projectedk = lM_projected{ct}(:,k);

                % Euler integration  Uk = (X_(k+1) - X_k)/*dt
                Xk = [ak; lMtildek];
                Zk = [a{ct}(:,k + 1);lMtilde{ct}(:,k + 1)];
                Uk = [ActivationDynamics(ek,ak,Misc.tau_act,Misc.tau_deact,Misc.b); vMtildek];
                opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);

                % Get muscle-tendon forces and derive Hill-equilibrium
                [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,...
                    vMtildek,lM_projectedk,DatStore(trial).LMTinterp(k,:)',...
                    Misc.params(:,Misc.idx_allMuscleList{trial})',...
                    Misc.kT(1,Misc.idx_allMuscleList{trial})',...
                    Misc.shift(1,Misc.idx_allMuscleList{trial})');
                % Hill-equilibrium constraint
                opti.subject_to(Hilldiffk == 0);

                % Add path constraints
                % Moment constraints
                for dof = 1:DatStore(trial).nDOF
                    T_exp = DatStore(trial).IDinterp(k,dof);
                    index_sel = (dof-1)*(NMuscles(trial))+1:(dof*NMuscles(trial)); % moment is a vector with the different dofs "below" each other
                    T_sim = DatStore(trial).MAinterp(k,index_sel)*FTk + Misc.Topt*aTk(dof);
                    opti.subject_to(T_exp - T_sim == 0);
                end
            end
            % Cost function
            J = J + Misc.wAct*0.5*(sumsqr(e{ct})/Mesh(trial).N/NMuscles(trial) + sumsqr(a{ct})/Mesh(trial).N/NMuscles(trial)) + ...
                Misc.wTres*sumsqr(aT{ct})/Mesh(trial).N/DatStore(trial).nDOF + ...
                Misc.wVm*sumsqr(vMtilde{ct})/Mesh(trial).N/NMuscles(trial);
        end    
        opti.minimize(J); % Define cost function in opti

        % Create an NLP solver
        opti.solver(output.setup.nlp.solver,optionssol);

        % Solve
        diary(fullfile(Misc.OutPath,[Misc.analysisName '_GenericMRS.txt']));
        tic
        sol = opti.solve();
        dt = toc;
        disp(['Computation time solving OCP: ' num2str(dt) ' s'])
        diary off

        ct = 0;
        for trial = Misc.trials_sel
            ct = ct + 1;
            % Extract results
            % Variables at mesh points
            % Muscle activations and muscle-tendon forces
            a_opt{trial} = sol.value(a{ct});
            lMtilde_opt{trial} = sol.value(lMtilde{ct});
            % Muscle excitations
            e_opt{trial} = sol.value(e{ct});
            % Reserve actuators
            aT_opt{trial} = sol.value(aT{ct});
            % Time derivatives of muscle-tendon forces
            vMtilde_opt{trial} = sol.value(vMtilde{ct});
            % Optimal lM_projectedilary variable
            lM_projected_opt{trial} = sol.value(lM_projected{ct});

            t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
            N = round((tf-t0)*Misc.Mesh_Frequency);
            % Time grid
            tgrid = linspace(t0,tf,N+1)';
            % Save results
            Results.Time(trial).genericMRS = tgrid;
            Results.lM_projected_opt(trial).genericMRS = lM_projected_opt{trial};
            Results.MActivation(trial).genericMRS = a_opt{trial};
            Results.lMtildeopt(trial).genericMRS = lMtilde_opt{trial};
            Results.lM(trial).genericMRS = lMtilde_opt{trial}.*repmat(Misc.lMo(:,Misc.idx_allMuscleList{trial})',1,length(tgrid));
            Results.vMtilde(trial).genericMRS = vMtilde_opt{trial};
            Results.MExcitation(trial).genericMRS = e_opt{trial};
            Results.RActivation(trial).genericMRS = aT_opt{trial}*Misc.Topt;
            Results.MuscleNames{trial} = DatStore(trial).MuscleNames;
            Results.OptInfo = output;
            % Tendon force
            Results.lMTinterp(trial).genericMRS = DatStore(trial).LMTinterp';
            [TForcetilde,TForce] = TendonForce_lMtilde(Results.lMtildeopt(trial).genericMRS',Misc.params(:,Misc.idx_allMuscleList{trial}),Results.lMTinterp(trial).genericMRS',Misc.kT(:,Misc.idx_allMuscleList{trial}),Misc.shift(:,Misc.idx_allMuscleList{trial}));    
            Results.TForcetilde(trial).genericMRS = TForcetilde';
            Results.TForce(trial).genericMRS = TForce';
            % get information F/l and F/v properties
            [Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).genericMRS',Results.vMtilde(trial).genericMRS',Misc.params(5,Misc.idx_allMuscleList{trial}));
            FMo = ones(N+1,1)*Misc.params(1,Misc.idx_allMuscleList{trial});
            Results.Fpe(trial).genericMRS = (Fpe.*FMo)';
            Results.FMltilde(trial).genericMRS = FMltilde';
            Results.FMvtilde(trial).genericMRS = FMvtilde';
        end
        
    
    else
        for trial = Misc.trials_sel
%             [Results] = runMRS(Misc,DatStore,Mesh,trial,output,optionssol,Results,NMuscles);

            % Problem bounds
            e_min = 0; e_max = 1;                   % bounds on muscle excitation
            a_min = 0; a_max = 1;                   % bounds on muscle activation
            vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
            lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

            SoActGuess = zeros(NMuscles(trial),Mesh(trial).N+1);
            SoExcGuess = zeros(NMuscles(trial),Mesh(trial).N);
            lMtildeGuess = zeros(NMuscles(trial),Mesh(trial).N+1);
            vMtildeGuess = zeros(NMuscles(trial),Mesh(trial).N);
            SoRActGuess = zeros(DatStore(trial).nDOF,Mesh(trial).N);

            SoActGuess = DatStore(trial).SoActInterp';
            SoExcGuess = DatStore(trial).SoActInterp(1:end-1,:)';
            lMtildeGuess = DatStore(trial).lMtildeInterp';
            vMtildeGuess = DatStore(trial).vMtildeinterp(1:end-1,:)';
            SoRActGuess = DatStore(trial).SoRActInterp(1:end-1,:)';

            % CasADi setup
            import casadi.*
            opti    = casadi.Opti();    % create opti structure

            a = opti.variable(NMuscles(trial),Mesh(trial).N+1);      % Variable at mesh points
            opti.subject_to(a_min < a < a_max);             % Bounds
            opti.set_initial(a,SoActGuess);                 % Initial guess (static optimization)
            %   - Muscle fiber lengths
            lMtilde = opti.variable(NMuscles(trial),Mesh(trial).N+1);
            opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
            opti.set_initial(lMtilde,lMtildeGuess);
            %   - Controls
            e = opti.variable(NMuscles(trial),Mesh(trial).N);
            opti.subject_to(e_min < e < e_max);
            opti.set_initial(e, SoExcGuess);
            %   - Reserve actuators
            aT = opti.variable(DatStore(trial).nDOF,Mesh(trial).N);
            opti.subject_to(-1 < aT <1);
            %   - Time derivative of muscle-tendon forces (states)
            vMtilde = opti.variable(NMuscles(trial),Mesh(trial).N);
            opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
            opti.set_initial(vMtilde,vMtildeGuess);    
            %   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
            lM_projected = opti.variable(NMuscles(trial),Mesh(trial).N+1);
            opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
            % Initial guess for this variable is retrieved from lMtilde guess
            % and geometric relationship between pennation angle, muscle length
            % and width
            lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
            alphao = Misc.params(4,Misc.idx_allMuscleList{trial})';
            lMGuess = lMtildeGuess.*lMo;
            w = lMo.*sin(alphao);
            lM_projectedGuess = sqrt((lMGuess.^2 - w.^2));
            opti.set_initial(lM_projected,lM_projectedGuess);

            % constraint on projected fiber length
            w = lMo.*sin(alphao);  
            lM = lMtilde.*lMo;
            opti.subject_to(lM.^2 - w.^2 == lM_projected.^2);

            % Discretization
            N = Mesh(trial).N;
            h = Mesh(trial).step;

            % Loop over mesh points formulating NLP
            for k=1:N
                % Variables within current mesh interval
                ak = a(:,k); lMtildek = lMtilde(:,k);
                vMtildek = vMtilde(:,k); aTk = aT(:,k);
                ek = e(:,k); lM_projectedk = lM_projected(:,k);

                % Euler integration  Uk = (X_(k+1) - X_k)/*dt
                Xk = [ak; lMtildek];
                Zk = [a(:,k + 1);lMtilde(:,k + 1)];
                Uk = [ActivationDynamics(ek,ak,Misc.tau_act,Misc.tau_deact,Misc.b); vMtildek];
                opti.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);

                % Get muscle-tendon forces and derive Hill-equilibrium
                [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,...
                    vMtildek,lM_projectedk,DatStore(trial).LMTinterp(k,:)',...
                    Misc.params(:,Misc.idx_allMuscleList{trial})',...
                    Misc.kT(:,Misc.idx_allMuscleList{trial})',...
                    Misc.shift(:,Misc.idx_allMuscleList{trial})');
                % Hill-equilibrium constraint
                opti.subject_to(Hilldiffk == 0);

                % Add path constraints
                % Moment constraints
                for dof = 1:DatStore(trial).nDOF
                    T_exp = DatStore(trial).IDinterp(k,dof);
                    index_sel = (dof-1)*(NMuscles(trial))+1:(dof*NMuscles(trial)); % moment is a vector with the different dofs "below" each other
                    T_sim = DatStore(trial).MAinterp(k,index_sel)*FTk + Misc.Topt*aTk(dof);
                    opti.subject_to(T_exp - T_sim == 0);
                end
            end
            % Cost function
            J = Misc.wAct*0.5*(sumsqr(e)/N/NMuscles(trial) + sumsqr(a)/N/NMuscles(trial)) + ...
                Misc.wTres*sumsqr(aT)/N/DatStore(trial).nDOF + ...
                Misc.wVm*sumsqr(vMtilde)/N/NMuscles(trial);

            opti.minimize(J); % Define cost function in opti

            % Create an NLP solver
            opti.solver(output.setup.nlp.solver,optionssol);

            % Solve
            diary(fullfile(Misc.OutPath,[Misc.OutName{trial} 'GenericMRS.txt']));
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

            t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
            N = round((tf-t0)*Misc.Mesh_Frequency);
            % Time grid
            tgrid = linspace(t0,tf,N+1)';
            % Save results
            Results.Time(trial).genericMRS = tgrid;
            Results.lM_projected_opt(trial).genericMRS = lM_projected_opt;
            Results.MActivation(trial).genericMRS = a_opt;
            Results.lMtildeopt(trial).genericMRS = lMtilde_opt;
            Results.lM(trial).genericMRS = lMtilde_opt.*repmat(Misc.lMo(:,Misc.idx_allMuscleList{trial})',1,length(tgrid));
            Results.vMtilde(trial).genericMRS = vMtilde_opt;
            Results.MExcitation(trial).genericMRS = e_opt;
            Results.RActivation(trial).genericMRS = aT_opt*Misc.Topt;
            Results.MuscleNames{trial} = DatStore(trial).MuscleNames;
            Results.OptInfo = output;
            % Tendon force
            Results.lMTinterp(trial).genericMRS = DatStore(trial).LMTinterp';
            [TForcetilde,TForce] = TendonForce_lMtilde(Results.lMtildeopt(trial).genericMRS',Misc.params(:,Misc.idx_allMuscleList{trial}),Results.lMTinterp(trial).genericMRS',Misc.kT(:,Misc.idx_allMuscleList{trial}),Misc.shift(:,Misc.idx_allMuscleList{trial}));    
            Results.TForcetilde(trial).genericMRS = TForcetilde';
            Results.TForce(trial).genericMRS = TForce';
            % get information F/l and F/v properties
            [Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).genericMRS',Results.vMtilde(trial).genericMRS',Misc.params(5,Misc.idx_allMuscleList{trial}));
            FMo = ones(N+1,1)*Misc.params(1,Misc.idx_allMuscleList{trial});
            Results.Fpe(trial).genericMRS = (Fpe.*FMo)';
            Results.FMltilde(trial).genericMRS = FMltilde';
            Results.FMvtilde(trial).genericMRS = FMvtilde';

        end
    end
end

%% Normalize EMG
if Misc.boolEMG
    maxEMG = nan(Misc.nAllMuscList,1);
    for t = Misc.trials_sel
        for m=1:DatStore(t).EMG.nEMG
            idx_m = Misc.idx_EMGsel{t}(m,1);
            maxEMG(idx_m,1) = max(maxEMG(idx_m,1),max(DatStore(t).EMG.EMGsel(:,m)));
        end
    end
    if Misc.normalizeToMRS
        maxMRS = nan(Misc.nAllMuscList,1);
        for t = Misc.trials_sel
            for m=1:NMuscles(t)
                idx_m = Misc.idx_allMuscleList{t}(m);
                maxMRS(idx_m,1) = max(maxMRS(idx_m,1),max(Results.MActivation(t).genericMRS(m,:)));
            end
        end
        for t = Misc.trials_sel
            for m=1:DatStore(t).EMG.nEMG
                idx_m = Misc.idx_EMGsel{t}(m,1);
                DatStore(t).EMG.EMGsel(:,m) = DatStore(t).EMG.EMGsel(:,m)./(maxEMG(idx_m,1)/maxMRS(idx_m,1));
            end
        end
    else
        for t = Misc.trials_sel
            for m=1:DatStore(t).EMG.nEMG
                idx_m = Misc.idx_EMGsel{t}(m,1);
                DatStore(t).EMG.EMGsel(:,m) = DatStore(t).EMG.EMGsel(:,m)./maxEMG(idx_m,1);
            end
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
%     [Results,Misc,DatStore,lMo_scaling_param_opt,lTs_scaling_param_opt,...
%         kT_scaling_param_opt,EMGscale_opt] = runParameterEstimation(Misc,...
%         DatStore,Mesh,output,optionssol,Results,NMuscles);

    % Problem bounds
    e_min = 0; e_max = 1;                   % bounds on muscle excitation
    a_min = 0; a_max = 1;                   % bounds on muscle activation
    vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
    lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

    [free_lMo,free_lTs,free_kT] = getFreeIndecies(Misc,DatStore,Misc.trials_sel);

    % Estimate parameters
    % CasADi setup
    import casadi.*
    opti_MTE = casadi.Opti();
    J = 0; % Initialize cost function

    % Free optimal fiber length
    lMo_scaling_param  = opti_MTE.variable(Misc.nAllMuscList,1);
    lb_lMo_scaling     = ones(Misc.nAllMuscList,1);   % default upper and lower bound is one (equality constraint)
    ub_lMo_scaling     = ones(Misc.nAllMuscList,1);   % default upper and lower bound is one (equality constraint)
    iM                 = free_lMo;        % index muscles with parameter estimation
    if ~isempty(iM)
        lb_lMo_scaling(iM) = Misc.lb_lMo_scaling;            % update lower bound for these muscles
        ub_lMo_scaling(iM) = Misc.ub_lMo_scaling;            % update uppder bound for these muscles
    end
    opti_MTE.subject_to(lb_lMo_scaling < lMo_scaling_param < ub_lMo_scaling); % update the upper and lower bounds

    % Free slack length
    lTs_scaling_param  = opti_MTE.variable(Misc.nAllMuscList,1);    
    lb_lTs_scaling     = ones(Misc.nAllMuscList,1); 
    ub_lTs_scaling     = ones(Misc.nAllMuscList,1);
    iM                 = free_lTs;        % index muscles with parameter estimation
    if ~isempty(iM)
        lb_lTs_scaling(iM) =  Misc.lb_lTs_scaling;
        ub_lTs_scaling(iM) =  Misc.ub_lTs_scaling;
    end
    opti_MTE.subject_to(lb_lTs_scaling < lTs_scaling_param < ub_lTs_scaling);

    % Free tendon stiffness
    kT_scaling_param   = opti_MTE.variable(Misc.nAllMuscList,1);
    lb_kT_scaling      = ones(Misc.nAllMuscList,1); 
    ub_kT_scaling      = ones(Misc.nAllMuscList,1);
    iM = free_kT;
    if ~isempty(iM)
        lb_kT_scaling(iM)  = Misc.lb_kT_scaling;
        ub_kT_scaling(iM)  = Misc.ub_kT_scaling;
    end
    opti_MTE.subject_to(lb_kT_scaling < kT_scaling_param < ub_kT_scaling);

    % added tendon stiffness coupling
    for k = 1:size(Misc.coupled_kT,1)
        opti_MTE.subject_to(kT_scaling_param(Misc.coupled_kT(k,1)) - kT_scaling_param(Misc.coupled_kT(k,2)) == 0);
    end    
    % added fibre length coupling
    for k = 1:size(Misc.coupled_lMo,1)
        opti_MTE.subject_to(lMo_scaling_param(Misc.coupled_lMo(k,1)) - lMo_scaling_param(Misc.coupled_lMo(k,2)) == 0);
    end
    % added tendon slack length coupling
    for k = 1:size(Misc.coupled_lTs,1)
        opti_MTE.subject_to(lTs_scaling_param(Misc.coupled_lTs(k,1)) - lTs_scaling_param(Misc.coupled_lTs(k,2)) == 0);
    end

    % Scale factor for EMG
    if Misc.boolEMG
        [DatStore,scaledEMGmusc] = getEMG_scaleIndecies(DatStore,Misc.trials_sel);
        EMGscale    = opti_MTE.variable(length(scaledEMGmusc),1);
        opti_MTE.subject_to(Misc.BoundsScaleEMG(1) < EMGscale < Misc.BoundsScaleEMG(2));
    end
    Misc.scaledEMGmusc = scaledEMGmusc;

    % Set initial guess
    opti_MTE.set_initial(lMo_scaling_param,1);
    opti_MTE.set_initial(lTs_scaling_param,1);
    opti_MTE.set_initial(kT_scaling_param,1);
    opti_MTE.set_initial(EMGscale,1);
    %%
    ct = 0;
    for trial = Misc.trials_sel
        ct = ct + 1;
        % States
        %   - Muscle activations
        a{ct} = opti_MTE.variable(NMuscles(trial),Mesh(trial).N+1);      % Variable at mesh points
        opti_MTE.subject_to(a_min < a{ct} < a_max);           % Bounds
        %   - Muscle fiber lengths
        lMtilde{ct} = opti_MTE.variable(NMuscles(trial),Mesh(trial).N+1);
        opti_MTE.subject_to(lMtilde_min < lMtilde{ct} < lMtilde_max);

        % Controls
        %   - Muscle excitations
        e{ct} = opti_MTE.variable(NMuscles(trial),Mesh(trial).N);
        opti_MTE.subject_to(e_min < e{ct} < e_max);
        %   - Reserve actuators
        aT{ct} = opti_MTE.variable(DatStore(trial).nDOF,Mesh(trial).N);
        opti_MTE.subject_to(-1 < aT{ct} <1);
        %   - Time derivative of muscle-tendon forces (states)
        vMtilde{ct} = opti_MTE.variable(NMuscles(trial),Mesh(trial).N);
        opti_MTE.subject_to(vMtilde_min < vMtilde{ct} < vMtilde_max);

        %   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
        lM_projected{ct} = opti_MTE.variable(NMuscles(trial),Mesh(trial).N+1);
        opti_MTE.subject_to(1e-4 < lM_projected{ct}(:)); % We impose that projected muscle fiber has strict positive length

        % Set initial guess
        if Misc.MRSBool == 1
            opti_MTE.set_initial(a{ct},Results.MActivation(trial).genericMRS);             % Initial guess generic MRS
            opti_MTE.set_initial(lMtilde{ct},Results.lMtildeopt(trial).genericMRS);
            opti_MTE.set_initial(e{ct},Results.MExcitation(trial).genericMRS);
            opti_MTE.set_initial(vMtilde{ct},Results.vMtilde(trial).genericMRS);
            opti_MTE.set_initial(aT{ct},Results.RActivation(trial).genericMRS./Misc.Topt);
            opti_MTE.set_initial(lM_projected{ct},Results.lM_projected_opt(trial).genericMRS);
        else
            SoActGuess = DatStore(trial).SoActInterp';
            SoExcGuess = DatStore(trial).SoActInterp(1:end-1,:)';
            lMtildeGuess = DatStore(trial).lMtildeInterp';
            vMtildeGuess = DatStore(trial).vMtildeinterp(1:end-1,:)';
            SoRActGuess = DatStore(trial).SoRActInterp(1:end-1,:)';
            opti_MTE.set_initial(a{ct},SoActGuess);             % Initial guess (static optimization)
            opti_MTE.set_initial(lMtilde{ct},lMtildeGuess);
            opti_MTE.set_initial(e{ct}, SoExcGuess);
            opti_MTE.set_initial(vMtilde{ct},vMtildeGuess);
            opti_MTE.set_initial(aT{ct},SoRActGuess./Misc.Topt);
            % Hill-type muscle model: geometric relationships
            lMo = Misc.params(2,Misc.idx_allMuscleList{trial})';
            alphao = Misc.params{trial}(4,Misc.idx_allMuscleList{trial})';
            w = lMo.*sin(alphao);
            lMGuess = lMtildeGuess.*lMo;
            lM_projectedGuess = sqrt((lMGuess.^2 - w.^2));
            opti_MTE.set_initial(lM_projected{ct},lM_projectedGuess);
        end

        % get EMG at optimization mesh
        if DatStore(trial).EMG.boolEMG
            EMGTracking(trial).data = interp1(DatStore(trial).EMG.time,DatStore(trial).EMG.EMGsel,Mesh(trial).t(1:end-1));
        end

        % get US data at optimization mesh
        if ~isempty(Misc.USfile)
            DatStore(trial).boolUS = 1;
            DatStore(trial).USTracking = interp1(DatStore(trial).US.time,DatStore(trial).US.USsel,Mesh(trial).t(1:end));
        end

        % constraint on projected fiber length
        lMo_all = lMo_scaling_param.*Misc.params(2,:)';
        alphao_all = Misc.params(4,:)';
        lMo = lMo_all(Misc.idx_allMuscleList{trial},1);
        alphao = alphao_all(Misc.idx_allMuscleList{trial},1);
        w = lMo.*sin(alphao);  
        lM = lMtilde{ct}.*lMo;
        opti_MTE.subject_to(lM.^2 - w.^2 == lM_projected{ct}.^2);

        % Time bounds
        N = Mesh(trial).N;
        h = Mesh(trial).step;
        for k=1:N
            % Variables within current mesh interval
            ak = a{ct}(:, k); lMtildek = lMtilde{ct}(:, k);
            vMtildek = vMtilde{ct}(:, k); aTk = aT{ct}(:, k); ek = e{ct}(:, k);
            lM_projectedk = lM_projected{ct}(:, k);

            % Euler integration  Uk = (X_(k+1) - X_k)/*dt
            Xk = [ak; lMtildek];
            Zk = [a{ct}(:,k + 1);lMtilde{ct}(:,k + 1)];
            Uk = [ActivationDynamics(ek,ak,Misc.tau_act,Misc.tau_deact,Misc.b); vMtildek];
            opti_MTE.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);

            % Get muscle-tendon forces and derive Hill-equilibrium
            [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree(ak,lMtildek,...
                vMtildek,lM_projectedk,DatStore(trial).LMTinterp(k,:)',...
                [lMo_scaling_param(Misc.idx_allMuscleList{trial},1) lTs_scaling_param(Misc.idx_allMuscleList{trial},1) kT_scaling_param(Misc.idx_allMuscleList{trial},1)],...
                Misc.params(:,Misc.idx_allMuscleList{trial})',Misc.kT(1,Misc.idx_allMuscleList{trial})');
            % Hill-equilibrium constraint
            opti_MTE.subject_to(Hilldiffk == 0);

            % Add path constraints
            % Moment constraints
            for dof = 1:DatStore(trial).nDOF
                T_exp = DatStore(trial).IDinterp(k,dof);
                index_sel = (dof-1)*(NMuscles(trial))+1:dof*NMuscles(trial);
                T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
                opti_MTE.subject_to(T_exp - T_sim == 0);
            end
        end

        % tracking lMtilde
        if DatStore(trial).US.boolUS
            lMo = lMo_scaling_param(DatStore(trial).US.idx_USsel(:,1))'.*Misc.params(2,DatStore(trial).US.idx_USsel(:,1));
            lMo = ones(size(DatStore(trial).USTracking,1),1)*lMo;
            lMtilde_tracking = DatStore(trial).USTracking./lMo/1000; % US data expected in mm in the input file.
            lMtilde_simulated = lMtilde{ct}(DatStore(trial).US.idx_USsel(:,3),:);
            if size(lMtilde_simulated,1) ~= size(lMtilde_tracking,1)
                lMtilde_simulated = lMtilde_simulated';
            end
            J = J + Misc.wlM*sumsqr(lMtilde_simulated-lMtilde_tracking)/DatStore(trial).US.nUS/N;
        end        
        % tracking Muscle activity
        if DatStore(trial).EMG.boolEMG
            eSim  = e{ct}(DatStore(trial).EMG.idx_EMGsel(:,3),:);
            eMeas = EMGTracking(trial).data' .* repmat(EMGscale(DatStore(trial).EMG.idx_EMGsel(:,4),1),1,N);
            JEMG  = sumsqr(eSim-eMeas);
            J     = J + Misc.wEMG * JEMG/DatStore(trial).EMG.nEMG/N;
        end        
        % Bounds on difference EMG - activations
        if DatStore(trial).EMG.boolEMG && ~isempty(Misc.EMGbounds)
            eSim  = e{ct}(DatStore(trial).EMG.idx_EMGsel(:,3),:);
            eMeas = EMGTracking(trial).data' .* repmat(EMGscale(DatStore(trial).EMG.idx_EMGsel(:,4),1),1,N);
            if Misc.EMGbounds(1) == 0
                opti_MTE.subject_to(eSim == eMeas); % "EMG driven"
            else
                opti_MTE.subject_to(Misc.EMGbounds(1) < eSim-eMeas < Misc.EMGbounds(2)); % "EMG constrained"
            end
        end

        J = J + ...
            Misc.wAct*0.5*(sumsqr(e{ct})/Mesh(trial).N/NMuscles(trial) + sumsqr(a{ct})/Mesh(trial).N/NMuscles(trial)) + ...
            Misc.wTres*sumsqr(aT{ct})/Mesh(trial).N/DatStore(trial).nDOF + ...
            Misc.wVm*sumsqr(vMtilde{ct})/Mesh(trial).N/NMuscles(trial);
    end
    %%
    opti_MTE.minimize(J); % Define cost function in opti
    opti_MTE.solver(output.setup.nlp.solver,optionssol);
    diary(fullfile(Misc.OutPath,[Misc.subjectName '_MTE.txt']));
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
    % get the results from w_opt vector (this is vector with all the
    % variables (constructed with opti.variable)
    iStart = 0; % help for extracting results
    lMo_scaling_param_opt = w_opt(iStart+1:iStart+Misc.nAllMuscList,1);
    iStart = iStart + Misc.nAllMuscList;
    % adapted optimal fiber length
    lMo_opt = lMo_scaling_param_opt.*Misc.params(2,:)';

    lTs_scaling_param_opt = w_opt(iStart+1:iStart+Misc.nAllMuscList,1);
    iStart = iStart + Misc.nAllMuscList;

    % Update muscle params
    paramsOpt =  Misc.params;
    paramsOpt(2,:) = paramsOpt(2,:).*lMo_scaling_param_opt'; % updated optimal fiber length
    paramsOpt(3,:) = paramsOpt(3,:).*lTs_scaling_param_opt';  % updated tendon slack length       
    Misc.paramsOpt = paramsOpt;

    kT_scaling_param_opt = w_opt(iStart+1:iStart+Misc.nAllMuscList,1);
    iStart = iStart + Misc.nAllMuscList;
    % Update muscle params
    kTOpt = Misc.kT .* kT_scaling_param_opt';
    shiftOpt = getShift(kTOpt);

    if Misc.boolEMG
        EMGscale_opt = w_opt(iStart+1:iStart+length(scaledEMGmusc),1);
        iStart = iStart + length(scaledEMGmusc);
    else
        EMGscale_opt = [];
    end

    for trial = Misc.trials_sel
        a_opt{trial} = reshape(w_opt(iStart+1:iStart+(NMuscles(trial)*(Mesh(trial).N+1)),1), NMuscles(trial),Mesh(trial).N+1);
        iStart = iStart + (NMuscles(trial)*(Mesh(trial).N+1));

        lMtilde_opt{trial} = reshape(w_opt(iStart+1:iStart+(NMuscles(trial)*(Mesh(trial).N+1)),1), NMuscles(trial),Mesh(trial).N+1);
        iStart = iStart + (NMuscles(trial)*(Mesh(trial).N+1));

        e_opt{trial} = reshape(w_opt(iStart+1:iStart+(NMuscles(trial)*Mesh(trial).N),1), NMuscles(trial),Mesh(trial).N);
        iStart = iStart + (NMuscles(trial)*Mesh(trial).N);

        aT_opt{trial} = reshape(w_opt(iStart+1:iStart+(DatStore(trial).nDOF*Mesh(trial).N),1), DatStore(trial).nDOF,Mesh(trial).N);
        iStart = iStart + (DatStore(trial).nDOF*Mesh(trial).N);

        vMtilde_opt{trial} = reshape(w_opt(iStart+1:iStart+(NMuscles(trial)*Mesh(trial).N),1), NMuscles(trial),Mesh(trial).N);
        iStart = iStart + (NMuscles(trial)*Mesh(trial).N);

        lM_projected_opt{trial} = reshape(w_opt(iStart+1:iStart+(NMuscles(trial)*(Mesh(trial).N+1)),1), NMuscles(trial),Mesh(trial).N+1);
        iStart = iStart + (NMuscles(trial)*(Mesh(trial).N+1));

        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        N = round((tf-t0)*Misc.Mesh_Frequency);
        % Time grid
        tgrid = linspace(t0,tf,N+1)';
        % Save results
        Results.Time(trial).MTE = tgrid;
        Results.MActivation(trial).MTE = a_opt{trial};
        Results.lMtildeopt(trial).MTE = lMtilde_opt{trial};
        Results.lM(trial).MTE = lMtilde_opt{trial}.*repmat(lMo_opt(Misc.idx_allMuscleList{trial},1),1,length(tgrid));
        Results.vMtilde(trial).MTE = vMtilde_opt{trial};
        Results.MExcitation(trial).MTE = e_opt{trial};
        Results.RActivation(trial).MTE = aT_opt{trial}*Misc.Topt;
        Results.lM_projected_opt(trial).MTE = lM_projected_opt{trial};

        Results.lMTinterp(trial).MTE = DatStore(trial).LMTinterp';
        % get force from muscle state
        [TForcetilde,TForce] = TendonForce_lMtilde(Results.lMtildeopt(trial).MTE',paramsOpt(:,Misc.idx_allMuscleList{trial}),Results.lMTinterp(trial).MTE',kTOpt(1,Misc.idx_allMuscleList{trial}),shiftOpt(1,Misc.idx_allMuscleList{trial}));
        Results.TForcetilde(trial).MTE = TForcetilde';
        Results.TForce(trial).MTE = TForce';
        [Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).MTE',Results.vMtilde(trial).MTE',Misc.params(5,Misc.idx_allMuscleList{trial}));
        FMo = ones(N+1,1)*Misc.params(1,Misc.idx_allMuscleList{trial});
        Results.Fpe(trial).MTE = (Fpe.*FMo)';
        Results.FMltilde(trial).MTE = FMltilde';
        Results.FMvtilde(trial).MTE = FMvtilde';
    end
    if iStart ~= length(w_opt)
        error('Problem with extracting results from NLP solution in parameter estimation');
    end


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
%         [Results] = runValidation_allTrialsTogether(Misc,DatStore,Mesh,output,optionssol,Results,NMuscles,lMo_scaling_param_opt,lTs_scaling_param_opt,kT_scaling_param_opt);

        % Problem bounds
        e_min = 0; e_max = 1;                   % bounds on muscle excitation
        a_min = 0; a_max = 1;                   % bounds on muscle activation
        vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
        lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

        % CasADi setup
        import casadi.*
        opti_validation = casadi.Opti();

        % Generate optimized parameters for specific trial by scaling generic parameters
        optimized_params = Misc.params';
        optimized_params(:,2) = lMo_scaling_param_opt.*optimized_params(:,2);
        optimized_params(:,3) = lTs_scaling_param_opt.*optimized_params(:,3);
        optimized_kT = kT_scaling_param_opt.*Misc.kT';
        optimized_shift = getShift(optimized_kT);

        % Loop over mesh points formulating NLP
        J = 0; % Initialize cost function

        ct = 0;
        for trial = Misc.trials_sel
            ct = ct + 1;
            % Variables - bounds and initial guess
            % States (at mesh and collocation points)
            % States
            %   - Muscle activations
            a{ct} = opti_validation.variable(NMuscles(trial),Mesh(trial).N+1);      % Variable at mesh points
            opti_validation.subject_to(a_min < a{ct} < a_max);           % Bounds
            opti_validation.set_initial(a{ct},Results.MActivation(trial).MTE);             % Initial guess

            %   - Muscle fiber lengths
            lMtilde{ct} = opti_validation.variable(NMuscles(trial),Mesh(trial).N+1);
            opti_validation.subject_to(lMtilde_min < lMtilde{ct} < lMtilde_max);
            opti_validation.set_initial(lMtilde{ct},Results.lMtildeopt(trial).MTE);

            % Controls
            %   - Muscle excitations
            e{ct} = opti_validation.variable(NMuscles(trial),Mesh(trial).N);
            opti_validation.subject_to(e_min < e{ct} < e_max);
            opti_validation.set_initial(e{ct}, Results.MExcitation(trial).MTE);
            %   - Reserve actuators
            aT{ct} = opti_validation.variable(DatStore(trial).nDOF,Mesh(trial).N);
            opti_validation.subject_to(-1 < aT{ct} <1);
            opti_validation.set_initial(aT{ct},Results.RActivation(trial).MTE/Misc.Topt);
            %   - Time derivative of muscle-tendon forces (states)
            vMtilde{ct} = opti_validation.variable(NMuscles(trial),Mesh(trial).N);
            opti_validation.subject_to(vMtilde_min < vMtilde{ct} < vMtilde_max);
            opti_validation.set_initial(vMtilde{ct},Results.vMtilde(trial).MTE);

            %   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
            lM_projected{ct} = opti_validation.variable(NMuscles(trial),Mesh(trial).N+1);
            opti_validation.subject_to(1e-4 < lM_projected{ct}(:)); % We impose that projected muscle fiber has strict positive length
            opti_validation.set_initial(lM_projected{ct},Results.lM_projected_opt(trial).MTE);    

            % constraint on projected fiber length
            lMo = optimized_params(Misc.idx_allMuscleList{trial},2);
            alphao = optimized_params(Misc.idx_allMuscleList{trial},4);
            lM = lMtilde{ct}.*lMo;
            w = lMo.*sin(alphao);
            opti_validation.subject_to(lM.^2 - w.^2 == lM_projected{ct}.^2);

            % Time bounds
            N = Mesh(trial).N;
            h = Mesh(trial).step;
            for k=1:N
                % Variables within current mesh interval
                ak = a{ct}(:, k); lMtildek = lMtilde{ct}(:, k);
                vMtildek = vMtilde{ct}(:, k); aTk = aT{ct}(:, k); ek = e{ct}(:, k);
                lM_projectedk = lM_projected{ct}(:, k);

                % Integration   Uk = (X_(k+1) - X_k)/*dt
                Xk = [ak; lMtildek];
                Zk = [a{ct}(:,k + 1);lMtilde{ct}(:,k + 1)];
                Uk = [ActivationDynamics(ek,ak,Misc.tau_act,Misc.tau_deact,Misc.b); vMtildek];
                opti_validation.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);

                % Impose that auxilary variable lM_projected behaves as defined
                [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,lM_projectedk,DatStore(trial).LMTinterp(k,:)',optimized_params(Misc.idx_allMuscleList{trial},:),optimized_kT(Misc.idx_allMuscleList{trial},1),optimized_shift(Misc.idx_allMuscleList{trial},1));
                % Hill-equilibrium constraint
                opti_validation.subject_to(Hilldiffk == 0);

                % Add path constraints
                % Moment constraints
                for dof = 1:DatStore(trial).nDOF
                    T_exp = DatStore(trial).IDinterp(k,dof);
                    index_sel = (dof-1)*(NMuscles(trial))+1:dof*NMuscles(trial);
                    T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
                    opti_validation.subject_to(T_exp - T_sim == 0);
                end
            end
            % Cost function
            J = J + ...
                0.5*(sumsqr(e{ct})/Mesh(trial).N/NMuscles(trial) + sumsqr(a{ct})/Mesh(trial).N/NMuscles(trial)) + ...
                Misc.wTres*sumsqr(aT{ct})/Mesh(trial).N/DatStore(trial).nDOF + ...
                Misc.wVm*sumsqr(vMtilde{ct})/Mesh(trial).N/NMuscles(trial);
        end
        opti_validation.minimize(J); % Define cost function in opti

        % Create an NLP solver
        opti_validation.solver(output.setup.nlp.solver,optionssol);

        % Solve
        diary(fullfile(Misc.OutPath,[Misc.subjectName '_ValidationMRS.txt']));
        tic
        sol = opti_validation.solve();
        dt = toc;
        disp(['Computation time solving OCP: ' num2str(dt) ' s'])
        diary off

        %% Extract results
        ct = 0;
        for trial = Misc.trials_sel
            ct = ct + 1;
            % Variables at mesh points
            % Muscle activations and muscle-tendon forces
            a_opt{trial} = sol.value(a{ct});
            lMtilde_opt{trial} = sol.value(lMtilde{ct});
            % Muscle excitations
            e_opt{trial} = sol.value(e{ct});
            % Reserve actuators
            aT_opt{trial} = sol.value(aT{ct});
            % Time derivatives of muscle-tendon forces
            vMtilde_opt{trial} = sol.value(vMtilde{ct});

            lM_projected_opt{trial} = sol.value(lM_projected{ct});

            % Save results
            t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
            N = round((tf-t0)*Misc.Mesh_Frequency);
            % Time grid
            tgrid = linspace(t0,tf,N+1)';
            % Save results
            Results.Time(trial).validationMRS = tgrid;
            Results.MActivation(trial).validationMRS = a_opt{trial};
            Results.lMtildeopt(trial).validationMRS = lMtilde_opt{trial};
            Results.lM(trial).validationMRS = lMtilde_opt{trial}.*repmat(optimized_params(Misc.idx_allMuscleList{trial},2),1,length(tgrid));
            Results.vMtilde(trial).validationMRS = vMtilde_opt{trial};
            Results.MExcitation(trial).validationMRS = e_opt{trial};
            Results.RActivation(trial).validationMRS = aT_opt{trial}*Misc.Topt;
            Results.lM_projected_opt(trial).validationMRS = lM_projected_opt{trial};
            % Tendon forces from lMtilde
            Results.lMTinterp(trial).validationMRS = DatStore(trial).LMTinterp';
            [TForcetilde,TForce] = TendonForce_lMtilde(Results.lMtildeopt(trial).validationMRS',optimized_params(Misc.idx_allMuscleList{trial},:)',Results.lMTinterp(trial).validationMRS',optimized_kT(Misc.idx_allMuscleList{trial},1)',optimized_shift(Misc.idx_allMuscleList{trial},1)');
            Results.TForcetilde(trial).validationMRS = TForcetilde';
            Results.TForce(trial).validationMRS = TForce';
            [Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).MTE',Results.vMtilde(trial).MTE',Misc.params(5,Misc.idx_allMuscleList{trial}));
            FMo = ones(N+1,1)*Misc.params(1,Misc.idx_allMuscleList{trial});
            Results.Fpe(trial).validationMRS = (Fpe.*FMo)';
            Results.FMltilde(trial).validationMRS = FMltilde';
            Results.FMvtilde(trial).validationMRS = FMvtilde';
        end

    
    else
        for trial = Misc.trials_sel
%             [Results] = runValidation(Misc,DatStore,Mesh,trial,output,optionssol,Results,NMuscles,lMo_scaling_param_opt,lTs_scaling_param_opt,kT_scaling_param_opt);

            % Problem bounds
            e_min = 0; e_max = 1;                   % bounds on muscle excitation
            a_min = 0; a_max = 1;                   % bounds on muscle activation
            vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
            lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

            % CasADi setup
            import casadi.*
            opti_validation = casadi.Opti();

            % Generate optimized parameters for specific trial by scaling generic parameters
            optimized_params = Misc.params';
            optimized_params(:,2) = lMo_scaling_param_opt.*optimized_params(:,2);
            optimized_params(:,3) = lTs_scaling_param_opt.*optimized_params(:,3);
            optimized_kT = kT_scaling_param_opt.*Misc.kT';
            optimized_shift = getShift(optimized_kT);

            % Variables - bounds and initial guess
            % States (at mesh and collocation points)
            % States
            %   - Muscle activations
            a = opti_validation.variable(NMuscles(trial),Mesh(trial).N+1);      % Variable at mesh points
            opti_validation.subject_to(a_min < a < a_max);           % Bounds
            opti_validation.set_initial(a,Results.MActivation(trial).MTE);             % Initial guess

            %   - Muscle fiber lengths
            lMtilde = opti_validation.variable(NMuscles(trial),Mesh(trial).N+1);
            opti_validation.subject_to(lMtilde_min < lMtilde < lMtilde_max);
            opti_validation.set_initial(lMtilde,Results.lMtildeopt(trial).MTE);

            % Controls
            %   - Muscle excitations
            e = opti_validation.variable(NMuscles(trial),Mesh(trial).N);
            opti_validation.subject_to(e_min < e < e_max);
            opti_validation.set_initial(e, Results.MExcitation(trial).MTE);
            %   - Reserve actuators
            aT = opti_validation.variable(DatStore(trial).nDOF,Mesh(trial).N);
            opti_validation.subject_to(-1 < aT <1);
            opti_validation.set_initial(aT,Results.RActivation(trial).MTE/Misc.Topt);
            %   - Time derivative of muscle-tendon forces (states)
            vMtilde = opti_validation.variable(NMuscles(trial),Mesh(trial).N);
            opti_validation.subject_to(vMtilde_min < vMtilde < vMtilde_max);
            opti_validation.set_initial(vMtilde,Results.vMtilde(trial).MTE);

            %   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
            lM_projected = opti_validation.variable(NMuscles(trial),Mesh(trial).N+1);
            opti_validation.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber has strict positive length
            opti_validation.set_initial(lM_projected,Results.lM_projected_opt(trial).MTE);

            % constraint on projected fiber length
            lMo = optimized_params(Misc.idx_allMuscleList{trial},2);
            alphao = optimized_params(Misc.idx_allMuscleList{trial},4);
            lM = lMtilde.*lMo;
            w = lMo.*sin(alphao);
            opti_validation.subject_to(lM.^2 - w.^2 == lM_projected.^2);

            % Time bounds
            N = Mesh(trial).N;
            h = Mesh(trial).step;
            for k=1:N
                % Variables within current mesh interval
                ak = a(:, k); lMtildek = lMtilde(:, k);
                vMtildek = vMtilde(:, k); aTk = aT(:, k); ek = e(:, k);
                lM_projectedk = lM_projected(:, k);

                % Integration   Uk = (X_(k+1) - X_k)/*dt
                Xk = [ak; lMtildek];
                Zk = [a(:,k + 1);lMtilde(:,k + 1)];
                Uk = [ActivationDynamics(ek,ak,Misc.tau_act,Misc.tau_deact,Misc.b); vMtildek];
                opti_validation.subject_to(eulerIntegrator(Xk,Zk,Uk,h) == 0);

                % Impose that auxilary variable lM_projected behaves as defined
                [Hilldiffk,FTk] = ForceEquilibrium_lMtildeState(ak,lMtildek,vMtildek,lM_projectedk,DatStore(trial).LMTinterp(k,:)',optimized_params(Misc.idx_allMuscleList{trial},:),optimized_kT(Misc.idx_allMuscleList{trial},1),optimized_shift(Misc.idx_allMuscleList{trial},1));
                % Hill-equilibrium constraint
                opti_validation.subject_to(Hilldiffk == 0);

                % Add path constraints
                % Moment constraints
                for dof = 1:DatStore(trial).nDOF
                    T_exp = DatStore(trial).IDinterp(k,dof);
                    index_sel = (dof-1)*(NMuscles(trial))+1:dof*NMuscles(trial);
                    T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
                    opti_validation.subject_to(T_exp - T_sim == 0);
                end
            end
            % Cost function
            J = 0.5*(sumsqr(e)/Mesh(trial).N/NMuscles(trial) + sumsqr(a)/Mesh(trial).N/NMuscles(trial)) + ...
                Misc.wTres*sumsqr(aT)/Mesh(trial).N/DatStore(trial).nDOF + ...
                Misc.wVm*sumsqr(vMtilde)/Mesh(trial).N/NMuscles(trial);
            opti_validation.minimize(J); % Define cost function in opti

            % Create an NLP solver
            opti_validation.solver(output.setup.nlp.solver,optionssol);

            % Solve
            diary(fullfile(Misc.OutPath,[Misc.subjectName '_ValidationMRS.txt']));
            tic
            sol = opti_validation.solve();
            dt = toc;
            disp(['Computation time solving OCP: ' num2str(dt) ' s'])
            diary off

            %% Extract results

            % Variables at mesh points
            % Muscle activations and muscle-tendon forces
            a_opt{trial} = sol.value(a);
            lMtilde_opt{trial} = sol.value(lMtilde);
            % Muscle excitations
            e_opt{trial} = sol.value(e);
            % Reserve actuators
            aT_opt{trial} = sol.value(aT);
            % Time derivatives of muscle-tendon forces
            vMtilde_opt{trial} = sol.value(vMtilde);

            lM_projected_opt{trial} = sol.value(lM_projected);

            % Save results
            t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
            N = round((tf-t0)*Misc.Mesh_Frequency);
            % Time grid
            tgrid = linspace(t0,tf,N+1)';
            % Save results
            Results.Time(trial).validationMRS = tgrid;
            Results.MActivation(trial).validationMRS = a_opt{trial};
            Results.lMtildeopt(trial).validationMRS = lMtilde_opt{trial};
            Results.lM(trial).validationMRS = lMtilde_opt{trial}.*repmat(optimized_params(Misc.idx_allMuscleList{trial},2),1,length(tgrid));
            Results.vMtilde(trial).validationMRS = vMtilde_opt{trial};
            Results.MExcitation(trial).validationMRS = e_opt{trial};
            Results.RActivation(trial).validationMRS = aT_opt{trial}*Misc.Topt;
            Results.lM_projected_opt(trial).validationMRS = lM_projected_opt{trial};
            % Tendon forces from lMtilde
            Results.lMTinterp(trial).validationMRS = DatStore(trial).LMTinterp';
            [TForcetilde,TForce] = TendonForce_lMtilde(Results.lMtildeopt(trial).validationMRS',optimized_params(Misc.idx_allMuscleList{trial},:)',Results.lMTinterp(trial).validationMRS',optimized_kT(Misc.idx_allMuscleList{trial},1)',optimized_shift(Misc.idx_allMuscleList{trial},1)');
            Results.TForcetilde(trial).validationMRS = TForcetilde';
            Results.TForce(trial).validationMRS = TForce';
            [Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).MTE',Results.vMtilde(trial).MTE',Misc.params(5,Misc.idx_allMuscleList{trial}));
            FMo = ones(N+1,1)*Misc.params(1,Misc.idx_allMuscleList{trial});
            Results.Fpe(trial).validationMRS = (Fpe.*FMo)';
            Results.FMltilde(trial).validationMRS = FMltilde';
            Results.FMvtilde(trial).validationMRS = FMvtilde';

        
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


