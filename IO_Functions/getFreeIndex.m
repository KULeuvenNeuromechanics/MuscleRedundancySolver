function [free_idx,free_Musc] = getFreeIndex(mNamesUnique,estMusc,coupleMusc,Misc)
% --------------------------------------------------------------------------
%getFreeIndex
%     This function finds all muscles for which paramters are to be estimated
%     based on both Misc.Estimate_<parameter> and Misc.Coupled_<parameter>
% 
% INPUT:
%     mNamesUnique
%     All muscles invloved in all trials
% 
%     estMusc
%     Misc.Estimate_<parameter>
%     
%     coupleMusc
%     Misc.Coupled_<parameter>
%     
%     Misc
%     Miscellaneous info used through the code
%     
% OUTPUT:
%     free_idx
%     Index of all muscles for which <paramter> is to be estimated
%     (index with respect to all muscles in the model as listed in Misc.allMuscleList)
% 
%     free_Musc
%     Names of all muscles for which <paramter> is to be estimated
% 
% Original author: Dhruv Gupta
% Original date: May 3, 2022
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

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
