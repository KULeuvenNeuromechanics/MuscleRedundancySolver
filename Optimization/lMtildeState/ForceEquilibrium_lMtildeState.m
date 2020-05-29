% Hill-type muscle model: equilibrium between muscle and tendon forces
% All muscle-tendon characteristics are fully described in the publication
% and its online supplement

function [err, FT, Fpe, FMltilde, FMvtilde, cos_alpha] = ForceEquilibrium_lMtildeState(a,lMtilde,vMtilde,aux,lMT,params,Atendon,shift)

FMo = params(:,1);
lMo = params(:,2);
lTs = params(:,3);

% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
lT = lMT - aux;
lTtilde = lT./lTs;

% Tendon force-length characteristic
fse = (exp(Atendon.*(lTtilde - 0.995)))/5-0.25+shift;

% get muscle force-length characteristic
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