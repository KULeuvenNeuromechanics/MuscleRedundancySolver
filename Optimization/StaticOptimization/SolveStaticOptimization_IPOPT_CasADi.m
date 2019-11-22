function DatStore = SolveStaticOptimization_IPOPT_CasADi(DatStore)

% STEP 1. Inputs
% --------------

time = DatStore.time;
N = length(time);
M = DatStore.nMuscles;
nDOF = DatStore.nDOF;

load('ActiveFVParameters.mat','ActiveFVParameters');
load('PassiveFLParameters','PassiveFLParameters');
load('Faparam.mat','Faparam');

act = ones(N,M);
FMltilde = ones(N,M);
FMvtilde = ones(N,M);
Fpe = ones(N,M);
cos_alpha = ones(N,M);
for m = 1:M
    pp_y = spline(time,DatStore.LMT(:,m));
    [LMTg,vMTg,~] = SplineEval_ppuval(pp_y,time,1);
    [~, ~, FMltilde(:,m), FMvtilde(:,m), Fpe(:,m), cos_alpha(:,m)] = ...
        HillModel_RigidTendon(act(:,m),LMTg,vMTg,DatStore.params(:,m),...
        ActiveFVParameters,PassiveFLParameters,Faparam);
    clear pp_y 
end
FMo = ones(size(act,1),1)*DatStore.Fiso;
Fpas = FMo.*Fpe.*cos_alpha;
Fact = FMo.*FMltilde.*FMvtilde.*cos_alpha;

% Add optimal force for reserve torques
Topt = 150/sqrt(1000);
Fpas = [Fpas zeros(N,nDOF)];
Fact = [Fact Topt*ones(N,nDOF)];

ID_data = DatStore.T_exp;

I = N*(M+nDOF);
MomentArms = zeros(I,nDOF);
temp = zeros(N,M);
for i = 1:nDOF    
    temp(:,:) = DatStore.dM(:,i,:);
    MomentArms(:,i) = ...
        reshape([temp zeros(N,i-1) ones(N,1) zeros(N,nDOF-i)]',I,1);    
end

% STEP 2. Optimization
% --------------------

import casadi.*
opti = casadi.Opti(); % Create opti instance
x = opti.variable((M+nDOF)*N); % Design variables
lb = repmat([zeros(M,1); -1500000*ones(nDOF,1)],N,1);
ub = repmat([ones(M,1); 1500000*ones(nDOF,1)],N,1);
opti.subject_to(lb < x < ub); % Bounds
auxdata = ...
    {M,N,nDOF,reshape(Fact',I,1),reshape(Fpas',I,1),ID_data,MomentArms};
opti.minimize(0.5*sumsqr(x));
[M, N, nDOF, Fmax, Fpas, ID_data, MomentArm] = deal(auxdata{:});
  F = Fmax .* x + Fpas;  
  c = MX(nDOF*N,1);  
  for k = 1:nDOF
      F_matrix = reshape(F, M+nDOF, N)';
      MomentArm_matrix = reshape(MomentArm(:,k), M+nDOF, N)';
      c(k:nDOF:end) = sum(F_matrix.*MomentArm_matrix, 2) - ID_data(:,k); 
  end
opti.subject_to(c == 0); % Constraints
% Settings
optionssol.ipopt.nlp_scaling_method = 'gradient-based';
optionssol.ipopt.linear_solver = 'mumps';
optionssol.ipopt.tol = 1e-6;
optionssol.ipopt.max_iter = 1000;
% Solver
opti.solver('ipopt',optionssol);
sol = opti.solve();
% Extract results
x_opt = reshape(sol.value(x), M+nDOF, N)';
act = x_opt(:,1:M);
eT = x_opt(:, M+1:M+nDOF)*Topt;
DatStore.SoAct = act;
DatStore.SoRAct = eT;
SoForce = FMo.*(act.*FMltilde.*FMvtilde + Fpe); 
DatStore.SoForce = SoForce;
DatStore.cos_alpha = cos_alpha;
