function [free_lMo,free_lTs,free_kT] = getFreeIndecies(Misc,DatStore,trials_sel)
mNames = [];
for t=trials_sel
    mNames = [mNames DatStore(t).MuscleNames];
end
mNamesUnique = unique(mNames);

free_lMo = getFreeIndex(mNamesUnique,Misc.Estimate_OptimalFiberLength,Misc.Coupled_fiber_length,Misc);
free_lTs = getFreeIndex(mNamesUnique,Misc.Estimate_TendonSlackLength,Misc.Coupled_slack_length,Misc);
free_kT = getFreeIndex(mNamesUnique,Misc.Estimate_TendonStiffness,Misc.Coupled_TendonStiffness,Misc);
