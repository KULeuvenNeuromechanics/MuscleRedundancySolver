function [Results] = runMRS(Misc,DatStore,Mesh,trial,output,optionssol,Results,NMuscles)
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