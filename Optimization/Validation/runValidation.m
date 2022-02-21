function [Results] = runValidation(Misc,DatStore,Mesh,trial,output,optionssol,Results,NMuscles,lMo_scaling_param_opt,lTs_scaling_param_opt,kT_scaling_param_opt)
%%
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