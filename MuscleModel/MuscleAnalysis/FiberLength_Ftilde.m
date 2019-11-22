% This function computes the muscle fiber length from the normalized tendon
% force

function [lM,lMtilde ] = FiberLength_Ftilde(Ftilde,params,lMT,Atendon,shift)

lMo = ones(size(Ftilde,1),1)*params(2,:);
lTs = ones(size(Ftilde,1),1)*params(3,:);
alphao = ones(size(Ftilde,1),1)*params(4,:);
Atendon = ones(size(Ftilde,1),1)*Atendon;
shift = ones(size(Ftilde,1),1)*shift;

% Non-linear tendon
lTtilde = real(log(5*(Ftilde + 0.25 - shift))./Atendon + 0.995);

% Hill-model relationship
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;
end

