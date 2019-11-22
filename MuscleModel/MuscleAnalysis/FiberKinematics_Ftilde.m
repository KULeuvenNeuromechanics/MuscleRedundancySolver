

function [lM,lMtilde,vM,vMtilde,lTtilde ] = FiberKinematics_Ftilde(Ftilde,dfse,params,lMT,vMT,Atendon,shift)
% This function computes the muscle fiber length and velocity with force as a state (force and rate of force)

% input arguments
lMo = ones(size(Ftilde,1),1)*params(2,:);
lTs = ones(size(Ftilde,1),1)*params(3,:);
alphao = ones(size(Ftilde,1),1)*params(4,:);
vMmax = ones(size(Ftilde,1),1)*params(5,:);
Atendon = ones(size(Ftilde,1),1)*Atendon;
shift   = ones(size(Ftilde,1),1)*shift;

% Non-linear tendon
lTtilde = real(log(5*(Ftilde + 0.25 - shift))./Atendon + 0.995);

% Hill-model relationship
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;

% Hill-model relationship
vT = lTs.*dfse./(0.2*Atendon.*exp(Atendon.*(lTtilde-0.995)));
cos_alpha = (lMT-lTs.*lTtilde)./lM;
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;
end

