

function shift = getShift(kT)
% This script returns the value used to shift the tendon force-length curve
% when changing the tendon stiffness. 
% With the standard stiffness (35), the shift is 0. For a different
% stiffness, the curve is shifted so that the normalized tendon force is
% the same as with the standard stiffness when the normalized tendon length
% is 1.
kT35 = 35;
lTtilde = 1;

fse_kt35 = (exp(kT35.*(lTtilde - 0.995)))/5 - 0.25; 
fse_kt = (exp(kT.*(lTtilde - 0.995)))/5 - 0.25; 
shift = fse_kt35-fse_kt;
end

