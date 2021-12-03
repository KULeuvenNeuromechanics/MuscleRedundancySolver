function [functions] = generateCasADiFunctions(Misc)
import casadi.*

% Euler integrator
X_MX = MX.sym('X_MX',Misc.NMuscles*2);
Z_MX = MX.sym('Z_MX',Misc.NMuscles*2);
U_MX = MX.sym('U_MX',Misc.NMuscles*2);
h = 1/Misc.Mesh_Frequency;

f_eulerIntegrator = Function('f_eulerIntegrator',{X_MX,Z_MX,U_MX},{eulerIntegrator(X_MX,Z_MX,U_MX,h)});  
functions.f_eulerIntegrator = f_eulerIntegrator;

% Activation dynamics
e_MX = MX.sym('e_MX',Misc.NMuscles);
a_MX = MX.sym('a_MX',Misc.NMuscles);

f_ActivationDynamics = Function('f_ActivationDynamics',{e_MX,a_MX},{ActivationDynamics(e_MX,a_MX,Misc.tauAct,Misc.tauDeact,Misc.b)});  
functions.f_ActivationDynamics = f_ActivationDynamics;

% ForceEquilibrium_lMtildeState
a_MX = MX.sym('a_MX',Misc.NMuscles);
lMtilde_MX = MX.sym('lMtilde_MX',Misc.NMuscles);     
vMtilde_MX = MX.sym('vMtilde_MX',Misc.NMuscles);           
lM_projected_MX = MX.sym('lM_projected_MX',Misc.NMuscles);
LMTinterp_MX = MX.sym('LMTinterp_MX',Misc.NMuscles);
params_MX = MX.sym('params_MX',Misc.NMuscles,5);
kT_MX = MX.sym('kT_MX',Misc.NMuscles);
shift_MX = MX.sym('shift_MX',Misc.NMuscles);

[Hilldiffk_MX,FT_MX] = ForceEquilibrium_lMtildeState(a_MX,lMtilde_MX,vMtilde_MX,lM_projected_MX,LMTinterp_MX,params_MX,kT_MX,shift_MX);

f_ForceEquilibrium_lMtildeState = Function('f_ForceEquilibrium_lMtildeState',{a_MX,lMtilde_MX,vMtilde_MX,lM_projected_MX,LMTinterp_MX,params_MX,kT_MX,shift_MX},{Hilldiffk_MX,FT_MX});  
functions.f_ForceEquilibrium_lMtildeState = f_ForceEquilibrium_lMtildeState;


% ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree
a_MX = MX.sym('a_MX',Misc.NMuscles);
lMtilde_MX = MX.sym('lMtilde_MX',Misc.NMuscles);     
vMtilde_MX = MX.sym('vMtilde_MX',Misc.NMuscles);           
lM_projected_MX = MX.sym('lM_projected_MX',Misc.NMuscles);
LMTinterp_MX = MX.sym('LMTinterp_MX',Misc.NMuscles);
lMo_scaling_param_MX = MX.sym('lMo_scaling_param_MX',Misc.NMuscles);
lTs_scaling_param_MX = MX.sym('lTs_scaling_param_MX',Misc.NMuscles);
kT_scaling_param_MX = MX.sym('kT_scaling_param_MX',Misc.NMuscles);
[Hilldiffk_MX,FT_MX] = ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree(a_MX,lMtilde_MX,vMtilde_MX,lM_projected_MX,LMTinterp_MX,[lMo_scaling_param_MX lTs_scaling_param_MX kT_scaling_param_MX],Misc.params',Misc.kT');
            
f_ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree = Function('f_ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree',{a_MX,lMtilde_MX,vMtilde_MX,lM_projected_MX,LMTinterp_MX,[lMo_scaling_param_MX lTs_scaling_param_MX kT_scaling_param_MX]},{Hilldiffk_MX,FT_MX});  
functions.f_ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree = f_ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree;


% compute T_sim
FT_MX = SX.sym('FT_MX',Misc.NMuscles);
aT_MX = SX.sym('aT_MX',Misc.nDOF);     
MAinterp_MX = SX.sym('MAinterp_MX',Misc.NMuscles*Misc.nDOF);           
f_compute_T_sim = Function('f_compute_T_sim',{FT_MX,aT_MX,MAinterp_MX},{compute_T_sim(FT_MX,aT_MX,MAinterp_MX,Misc)});  
functions.f_compute_T_sim = f_compute_T_sim;
end

