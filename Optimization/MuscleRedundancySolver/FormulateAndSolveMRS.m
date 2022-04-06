function [Results] = FormulateAndSolveMRS(Misc,DatStore,Mesh,trial,SolverSetup,...
    Results,NMuscles,IG,MuscProperties,PrefixOutFile)
% Problem bounds
e_min = 0; e_max = 1;                   % bounds on muscle excitation
a_min = 0; a_max = 1;                   % bounds on muscle activation
vMtilde_min = -10; vMtilde_max = 10;    % bounds on normalized muscle fiber velocity
lMtilde_min = 0.1; lMtilde_max = 1.7;   % bounds on normalized muscle fiber length

% CasADi setup
import casadi.*
opti    = casadi.Opti();    % create opti structure

a = opti.variable(NMuscles(trial),Mesh(trial).N+1);      % Variable at mesh points
opti.subject_to(a_min < a < a_max);             % Bounds
opti.set_initial(a,IG.a);                 % Initial guess (static optimization)
%   - Muscle fiber lengths
lMtilde = opti.variable(NMuscles(trial),Mesh(trial).N+1);
opti.subject_to(lMtilde_min < lMtilde < lMtilde_max);
opti.set_initial(lMtilde,IG.lMtilde);
%   - Controls
e = opti.variable(NMuscles(trial),Mesh(trial).N);
opti.subject_to(e_min < e < e_max);
opti.set_initial(e, IG.e);
%   - Reserve actuators
aT = opti.variable(DatStore(trial).nDOF,Mesh(trial).N);
opti.subject_to(-1 < aT <1);
if isfield(IG,'aT') && ~isempty(IG.aT)
    opti.set_initial(aT,IG.aT);
end
%   - Time derivative of muscle-tendon forces (states)
vMtilde = opti.variable(NMuscles(trial),Mesh(trial).N);
opti.subject_to(vMtilde_min < vMtilde < vMtilde_max);
opti.set_initial(vMtilde,IG.vMtilde);    
%   - Projected muscle fiber length - Auxilary variable to avoid muscle buckling & square root expression in muscle dynamics
lM_projected = opti.variable(NMuscles(trial),Mesh(trial).N+1);
opti.subject_to(1e-4 < lM_projected(:)); % We impose that projected muscle fiber length has strict positive length
opti.set_initial(lM_projected,IG.lM_projected);
% Initial guess for this variable is retrieved from lMtilde guess
% and geometric relationship between pennation angle, muscle length
% and width
% constraint on projected fiber length
lMo = MuscProperties.params(Misc.idx_allMuscleList{trial},2);
alphao = MuscProperties.params(Misc.idx_allMuscleList{trial},4);
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
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
        MuscProperties.params(Misc.idx_allMuscleList{trial},:),...
        MuscProperties.kT(Misc.idx_allMuscleList{trial},1),...
        MuscProperties.shift(Misc.idx_allMuscleList{trial},1));
    % Hill-equilibrium constraint
    opti.subject_to(Hilldiffk == 0);

    % Add path constraints
    % Moment constraints
    for dof = 1:DatStore(trial).nDOF
        T_exp = DatStore(trial).IDinterp(k,dof);
        index_sel = (dof-1)*(NMuscles(trial))+1:(dof*NMuscles(trial)); % moment is a vector with the different dofs "below" each other
        T_sim = FTk'*DatStore(trial).MAinterp(k,index_sel)' + Misc.Topt*aTk(dof);
        opti.subject_to(T_exp - T_sim == 0);
    end
end
% Cost function
J = Misc.wAct*0.5*(sumsqr(e)/N/NMuscles(trial) + sumsqr(a)/N/NMuscles(trial)) + ...
    Misc.wTres*sumsqr(aT)/N/DatStore(trial).nDOF + ...
    Misc.wVm*sumsqr(vMtilde)/N/NMuscles(trial);

opti.minimize(J); % Define cost function in opti
    
% Create an NLP solver
opti.solver(SolverSetup.nlp.solver,SolverSetup.optionssol);
    
% Solve
diary(fullfile(Misc.OutPath,[Misc.OutName '_' PrefixOutFile '.txt']));
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
saveName = PrefixOutFile;

% append the results structure
Results.Time(trial).(saveName) = tgrid;
Results.lM_projected_opt(trial).(saveName) = lM_projected_opt;
Results.MActivation(trial).(saveName) = a_opt;
Results.lMtildeopt(trial).(saveName) = lMtilde_opt;
Results.lM(trial).(saveName) = lMtilde_opt.*repmat(MuscProperties.params(Misc.idx_allMuscleList{trial},2),1,length(tgrid));
Results.vMtilde(trial).(saveName) = vMtilde_opt;
Results.MExcitation(trial).(saveName) = e_opt;
Results.RActivation(trial).(saveName) = aT_opt*Misc.Topt;
Results.MuscleNames{trial} = DatStore(trial).MuscleNames;
Results.OptInfo = SolverSetup;
% Tendon force
Results.lMTinterp(trial).(saveName) = DatStore(trial).LMTinterp';
[TForcetilde,TForce] = TendonForce_lMtilde(Results.lMtildeopt(trial).(saveName)',...
    MuscProperties.params(Misc.idx_allMuscleList{trial},:)',Results.lMTinterp(trial).(saveName)',...
    MuscProperties.kT(Misc.idx_allMuscleList{trial},1)',MuscProperties.shift(Misc.idx_allMuscleList{trial},1)');    
Results.TForcetilde(trial).(saveName) = TForcetilde';
Results.TForce(trial).(saveName) = TForce';
% get information F/l and F/v properties
[Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(Results.lMtildeopt(trial).(saveName)',...
    Results.vMtilde(trial).(saveName)',MuscProperties.params(Misc.idx_allMuscleList{trial},5)');
FMo = ones(N+1,1)*Misc.params(1,Misc.idx_allMuscleList{trial});
Results.Fpe(trial).(saveName) = (Fpe.*FMo)';
Results.FMltilde(trial).(saveName) = FMltilde';
Results.FMvtilde(trial).(saveName) = FMvtilde';
end