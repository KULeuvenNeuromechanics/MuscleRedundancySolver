function [free_idx,free_Musc] = getFreeIndex(mNamesUnique,estMusc,coupleMusc,Misc)

if ~isempty(estMusc)
    estMusc_all = unique([estMusc; reshape(coupleMusc,numel(coupleMusc),1)]);
    mOpt = intersect(mNamesUnique,estMusc_all');

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
else
    free_idx = [];
    free_Musc = [];
end
