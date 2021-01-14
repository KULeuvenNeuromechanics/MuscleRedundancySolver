function [] = Warnings_FiberLength_RigidTendon(lmtilde,Mnames)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% get the minmal and maximal norm fiber lengths
lm_min = round(min(lmtilde),3);
lm_max = round(max(lmtilde),3);

% display muscle names with short or long length

iSel = find(lm_min< 0.4 | lm_max > 1.5);
nM = length(iSel);

% display warnings
for i=1:nM
    m = iSel(i);
   disp([' Warning muscle length rigid tendon - ' Mnames{m}, ...
       ' Min-Max: ' num2str(lm_min(m)) ' ' num2str(lm_max(m))]);
    
end



end

