% This script contains several CasADi-based functions that are
% used when solving the OCPs
import casadi.*

%% Functions for sum of squared values
% # muscles
e_ss = SX.sym('etemp_NMuscles',auxdata.NMuscles);
J_ss = 0;
for i=1:length(e_ss)
    J_ss = J_ss + e_ss(i).^2;
end
f_ssNMuscles = Function('f_ssNMuscles',{e_ss},{J_ss});
% # dofs
e_ss = SX.sym('e_Ndof',auxdata.Ndof);
J_ss = 0;
for i=1:length(e_ss)
    J_ss = J_ss + e_ss(i).^2;
end
f_ssNdof = Function('f_ssNdof',{e_ss},{J_ss});

%% Function for sum of products 
% # muscles
ma_sp = SX.sym('ma_NMuscles',auxdata.NMuscles);
ft_sp = SX.sym('ft_NMuscles',auxdata.NMuscles);
J_sp = 0;
for i=1:length(ma_sp)
    J_sp = J_sp + ma_sp(i,1)*ft_sp(i,1);    
end
f_spNMuscles = Function('f_spNMuscles',{ma_sp,ft_sp},{J_sp});

%% Muscle contraction dynamics: tendon force Ft as a state
% States 
FTtilde_SX = SX.sym('FTtilde',auxdata.NMuscles); % Normalized tendon forces
a_SX = SX.sym('a',auxdata.NMuscles); % Muscle activations
% Controls 
dFTtilde_SX = SX.sym('dFTtilde',auxdata.NMuscles); % Derivative tendon forces
% Results
Hilldiff_SX = SX(auxdata.NMuscles,1); % Hill equilibrium
FT_SX = SX(auxdata.NMuscles,1); % Tendon forces
lMT_SX = SX.sym('lMT',auxdata.NMuscles); % Muscle-tendon lengths
vMT_SX = SX.sym('vMT',auxdata.NMuscles); % Muscle-tendon velocities
for m = 1:auxdata.NMuscles 
    [Hilldiff_SX(m),FT_SX(m)] = ForceEquilibrium_FtildeState(a_SX(m),FTtilde_SX(m),dFTtilde_SX(m),lMT_SX(m),vMT_SX(m),auxdata.params(:,m),auxdata.Fvparam,auxdata.Fpparam,auxdata.Faparam,auxdata.Atendon(m),auxdata.shift(m));        
end
f_forceEquilibrium_FtildeState = Function('f_forceEquilibrium_FtildeState',{a_SX,FTtilde_SX,dFTtilde_SX,lMT_SX,vMT_SX,},{Hilldiff_SX,FT_SX});
% Test function
% atest = rand(auxdata.NMuscles,1);
% FTtildetest = rand(auxdata.NMuscles,1);
% dFTtildetest = rand(auxdata.NMuscles,1);
% lMTtest = rand(auxdata.NMuscles,1);
% vMTtest = rand(auxdata.NMuscles,1);
% for m = 1:auxdata.NMuscles
%     [Hilldifftest(m,1),FTtest(m,1)] = ForceEquilibrium_FtildeState(...
%         atest(m),FTtildetest(m),dFTtildetest(m),lMTtest(m),...
%         vMTtest(m),auxdata.params(:,m),auxdata.Fvparam,auxdata.Fpparam,auxdata.Faparam,auxdata.Atendon(m));
% end
% [Hilldifftest2,FTtest2] = f_forceEquilibrium_FtildeState(atest,FTtildetest,dFTtildetest,lMTtest,vMTtest);
% assertHill = max(abs(Hilldifftest-Hilldifftest2));
% assertF = max(abs(FTtest-FTtest2));

%% Muscle contraction dynamics: normalized muscle fiber length as a state
% States 
lMtilde_SX = SX.sym('lMtilde',auxdata.NMuscles); % Normalized muscle fiber lengths
a_SX = SX.sym('a',auxdata.NMuscles); % Muscle activations
% Controls 
vMtilde_SX = SX.sym('vMtilde',auxdata.NMuscles); % Derivative muscle fiber lengths
% Results
Hilldiff_SX = SX(auxdata.NMuscles,1); % Hill equilibrium
FT_SX = SX(auxdata.NMuscles,1); % Tendon forces
lMT_SX = SX.sym('lMT',auxdata.NMuscles); % Muscle-tendon lengths
for m = 1:auxdata.NMuscles     
    [Hilldiff_SX(m),FT_SX(m)] = ForceEquilibrium_lMtildeState(a_SX(m),lMtilde_SX(m),vMtilde_SX(m),lMT_SX(m),auxdata.params(:,m),auxdata.Fvparam,auxdata.Fpparam,auxdata.Faparam,auxdata.Atendon(m),auxdata.shift(m));
end
f_forceEquilibrium_lMtildeState = Function('f_forceEquilibrium_lMtildeState',{a_SX,lMtilde_SX,vMtilde_SX,lMT_SX},{Hilldiff_SX,FT_SX});
% Test function
% atest = rand(auxdata.NMuscles,1);
% lMtildetest = rand(auxdata.NMuscles,1);
% vMtildetest = rand(auxdata.NMuscles,1);
% lMTtest = rand(auxdata.NMuscles,1);
% for m = 1:auxdata.NMuscles
%     [Hilldifftest(m,1),FTtest(m,1)] = ForceEquilibrium_lMtildeState(...
%         atest(m),lMtildetest(m),vMtildetest(m),lMTtest(m),...
%         auxdata.params(:,m),auxdata.Fvparam,auxdata.Fpparam,auxdata.Faparam,auxdata.Atendon(m));
% end
% [Hilldifftest2,FTtest2] = f_forceEquilibrium_lMtildeState(atest,lMtildetest,vMtildetest,lMTtest);
% assertHill = max(abs(Hilldifftest-Hilldifftest2));
% assertF = max(abs(FTtest-FTtest2));

%% Muscle activation dynamics
e_SX = SX.sym('e',auxdata.NMuscles); % Muscle excitations
a_SX = SX.sym('a',auxdata.NMuscles); % Muscle activations
dadt_SX = SX(auxdata.NMuscles,1); % Hill equilibrium
for m = 1:auxdata.NMuscles
    dadt_SX(m) = ActivationDynamics(e_SX(m),a_SX(m),auxdata.tauAct(m),auxdata.tauDeact(m),auxdata.b);
end
f_ActivationDynamics = Function('f_ActivationDynamics',{e_SX,a_SX},{dadt_SX});
% Test function
% etest = rand(auxdata.NMuscles,1);
% atest = rand(auxdata.NMuscles,1);
% for m = 1:auxdata.NMuscles
%     [datdttest(m,1)] = ActivationDynamics(etest(m),atest(m),auxdata.tauAct(m),auxdata.tauDeact(m),auxdata.b);
% end
% datdttest2 = f_ActivationDynamics(etest,atest);
% assertdadt = max(abs(datdttest-datdttest2));
