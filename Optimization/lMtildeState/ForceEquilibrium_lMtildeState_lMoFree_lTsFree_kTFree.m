% Hill-type muscle model: equilibrium between muscle and tendon forces
% All muscle-tendon characteristics are fully described in the publication
% and its online supplement

function [err, FT] = ForceEquilibrium_lMtildeState_lMoFree_lTsFree_kTFree(a,lMtilde,vMtilde,lMT,lMo_lTs_kT_scaling,params,Atendon)

FMo = params(:,1);
lMo = lMo_lTs_kT_scaling(:,1).*params(:,2);
lTs = lMo_lTs_kT_scaling(:,2).*params(:,3);
alphao = params(:,4);
Atendon = lMo_lTs_kT_scaling(:,3).*Atendon;

shift = (exp(Atendon.*(1 - 0.995)))/5 - (exp(35.*(1 - 0.995)))/5;


% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
lT = lMT - sqrt((lM.^2 - w.^2));
lTtilde = lT./lTs;

% Tendon force-length characteristic
fse = (exp(Atendon.*(lTtilde - 0.995)))/5-0.25+shift;

% get muscle force-length-velocity characteristic
[Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lMtilde,vMtilde);

% Active muscle force
Fce = a.*FMltilde.*FMvtilde;

% Muscle force
FM = Fce+Fpe;
% Tendon force
FT = FMo.*fse;

% Equilibrium between muscle and tendon forces
% Fm*cos(alpha) = Ft
cos_alpha = (lMT-lT)./lM;
err =  FM.*cos_alpha-fse;

end