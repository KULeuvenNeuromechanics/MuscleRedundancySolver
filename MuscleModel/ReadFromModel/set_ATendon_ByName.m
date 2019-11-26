function [Misc,DatStore] = set_ATendon_ByName(Misc,DatStore)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nit = length(Misc.Set_ATendon_ByName(:,1));
for i=1:nit
    NameSel = Misc.Set_ATendon_ByName{i,1};
    k = Misc.Set_ATendon_ByName{i,2};
    IndMus = strcmp(NameSel,DatStore.MuscleNames);
    if any(IndMus)
        Misc.Atendon(IndMus)= k;
    else
        disp(['Cannot set stiffness of ' NameSel ' because this muscle in not used']);
    end
end
end

