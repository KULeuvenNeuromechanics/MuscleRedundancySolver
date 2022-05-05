function [Misc] = set_kT_ByName(Misc)
% Set stifness of muscles that was user provided (not default)
nit = length(Misc.Set_kT_ByName(:,1));
for i=1:nit
    NameSel = Misc.Set_kT_ByName{i,1};
    kT = Misc.Set_kT_ByName{i,2};
    IndMus = strcmp(NameSel,Misc.allMuscleList);
    if any(IndMus)
        Misc.kT(IndMus)= kT;
    else
        warning(['Cannot set stiffness of ' NameSel ' because this muscle in not used']);
        disp(' ')
    end
end
end

