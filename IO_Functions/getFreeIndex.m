function [free_idx,free_Musc] = getFreeIndex(mNamesUnique,estMusc,coupleMusc,Misc)
estMusc_all = unique([estMusc reshape(coupleMusc,1,numel(coupleMusc))]);
mOpt = intersect(mNamesUnique,estMusc_all);

free_idx = nan(1,length(mOpt));
free_Musc = cell(1,length(mOpt));
ct=0;
for m=1:Misc.nAllMuscList
    if any(ismember(mOpt,Misc.allMuscleList{m}))
        ct=ct+1;
        free_idx(ct) = m;
        free_Musc{ct} = Misc.allMuscleList{m};
    end
end
