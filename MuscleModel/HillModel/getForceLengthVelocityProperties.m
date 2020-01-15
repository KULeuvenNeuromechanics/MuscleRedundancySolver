function [Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties(lMtilde,vMtilde,varargin)
% gets the force-length and force-velocity properties of DeGroote 2016 muscle model 

% Make plot with the F/L and F/v relation ?
boolPlot= 0;
if ~isempty(varargin)
    boolPlot = varargin{1};
end


% Parameters of active muscle force-velocity characteristic
ActiveFVParameters = [-0.215812499592629,-32.5966241739010,-1.12412150864786,0.912644059915004];
Fvparam(1) = 1.475*ActiveFVParameters(1); Fvparam(2) = 0.25*ActiveFVParameters(2);
Fvparam(3) = ActiveFVParameters(3) + 0.75; Fvparam(4) = ActiveFVParameters(4) - 0.027;

% Parameters of active muscle force-length characteristic
Faparam = [0.814483478343008;1.05503342897057;0.162384573599574;0.0633034484654646;0.433004984392647;0.716775413397760;-0.0299471169706956;0.200356847296188];

% Parameters of passive muscle force-length characteristic
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
Fpparam = [pp1;pp2];

% compute the % of the maximal force generating capacity at the current
% state

% get passive force
e0 = 0.6;
kpe = 4;
t5 = exp(kpe * (lMtilde - 0.10e1) / e0);
Fpe = ((t5 - 0.10e1) - Fpparam(1)) / Fpparam(2);

% get active force-length
% Active muscle force-length characteristic
b11 = Faparam(1);
b21 = Faparam(2);
b31 = Faparam(3);
b41 = Faparam(4);
b12 = Faparam(5);
b22 = Faparam(6);
b32 = Faparam(7);
b42 = Faparam(8);

b13 = 0.1;
b23 = 1;
b33 = 0.5*sqrt(0.5);
b43 = 0;
num3 = lMtilde-b23;
den3 = b33+b43*lMtilde;
FMtilde3 = b13*exp(-0.5*num3.^2./den3.^2);

num1 = lMtilde-b21;
den1 = b31+b41*lMtilde;
FMtilde1 = b11*exp(-0.5*num1.^2./den1.^2);

num2 = lMtilde-b22;
den2 = b32+b42*lMtilde;
FMtilde2 = b12*exp(-0.5*num2.^2./den2.^2);

FMltilde = FMtilde1+FMtilde2+FMtilde3;

% get active force-velocity
e1 = Fvparam(1);
e2 = Fvparam(2);
e3 = Fvparam(3);
e4 = Fvparam(4);
FMvtilde = e1*log((e2*vMtilde+e3)+sqrt((e2*vMtilde+e3).^2+1))+e4;




%% Plot default graph
% plot FL and FV properties
if boolPlot
  
   figure();
   subplot(1,2,1)
   plot(lMtilde,FMltilde,'k'); hold on;
   plot(lMtilde,Fpe,'b'); legend('active','passive');
   xlabel('norm fiber length');
   ylabel('maximal force production');
   
   subplot(1,2,2)
   plot(vMtilde,FMvtilde ,'k'); hold on;
   xlabel('norm fiber velocity');
   ylabel('maximal force production');
   
end

end