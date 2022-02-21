function [Results,Misc,DatStore,lMo_scaling_param_opt,lTs_scaling_param_opt,kT_scaling_param_opt,EMGscale_opt] = runParameterEstimation(Misc,DatStore,Mesh,output,optionssol,Results,NMuscles)
%%
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
% opti_MTE.set_initial(EMGscale,1);
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

end