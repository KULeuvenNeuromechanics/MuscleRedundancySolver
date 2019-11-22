% This function computes the tendon force from the normalized muscle fiber
% length

function [Ftilde,F] = TendonForce_lMtilde(lMtilde,params,lMT,Atendon,shift)

FMo = ones(size(lMtilde,1),1)*params(1,:);
lMo = ones(size(lMtilde,1),1)*params(2,:);
lTs = ones(size(lMtilde,1),1)*params(3,:);
alphao = ones(size(lMtilde,1),1)*params(4,:);
Atendon = ones(size(lMtilde,1),1)*Atendon;
shift = ones(size(lMtilde,1),1)*shift;

% Hill-type muscle model: geometric relationships
lM = lMtilde.*lMo;
w = lMo.*sin(alphao);
lT = lMT - sqrt((lM.^2 - w.^2));
lTtilde = lT./lTs;

% Tendon force-length characteristic
Ftilde = (exp(Atendon.*(lTtilde - 0.995)))/5-0.25+shift;
F = FMo.*Ftilde;

end
