function [Misc,DatStore] = set_kT_ByName(Misc,DatStore)
% Set stifness of different muscles
nit = length(Misc.Set_kT_ByName(:,1));
for i=1:nit
    NameSel = Misc.Set_kT_ByName{i,1};
    kT = Misc.Set_kT_ByName{i,2};
    IndMus = strcmp(NameSel,DatStore.MuscleNames);
    if any(IndMus)
        Misc.kT(IndMus)= kT;
    else
        disp(['Cannot set stiffness of ' NameSel ' because this muscle in not used']);
    end
end
end

