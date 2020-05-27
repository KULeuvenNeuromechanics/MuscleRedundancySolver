% Hill-type muscle model: equilibrium between muscle and tendon forces
% All muscle-tendon characteristics are fully described in the publication
% and its online supplement

function [err, FT] = ForceEquilibrium_lMtildeState(a,lMtilde,vMtilde,lMT,params,Atendon,shift)

FMo = params(:,1);
lMo = params(:,2);
lTs = params(:,3);
alphao = params(:,4);
Atendon = Atendon;
shift = shift;

% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
lT = lMT - sqrt((lM.^2 - w.^2));
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