function [free_lMo,free_lTs,free_kT] = getFreeIndecies(Misc,DatStore)
% --------------------------------------------------------------------------
%getFreeIndecies
%     This function finds all muscles for which paramters are to be estimated
%     based on both Misc.Estimate_<parameter> and Misc.Coupled_<parameter>
% 
% INPUT:
%     Misc
%     Miscellaneous info used through the code
% 
%     DatStore
%     Structure of all data
%     
% OUTPUT:
%     free_lMo
%     All muscles for which optimal fiber length is to be estimated
% 
%     free_lTs
%     All muscles for which tendon slack length is to be estimated
% 
%     free_kT
%     All muscles for which tendon stiffness is to be estimated
% 
% Original author: Dhruv Gupta
% Original date: May 3, 2022
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

mNames = [];
for t=Misc.trials_sel
    mNames = [mNames DatStore(t).MuscleNames];
end
mNamesUnique = unique(mNames);

% optimal fiber length
free_lMo = getFreeIndex(mNamesUnique,Misc.Estimate_OptimalFiberLength,Misc.Coupled_fiber_length,Misc);
% tendon slack length
free_lTs = getFreeIndex(mNamesUnique,Misc.Estimate_TendonSlackLength,Misc.Coupled_slack_length,Misc);
% tendon stiffness
free_kT = getFreeIndex(mNamesUnique,Misc.Estimate_TendonStiffness,Misc.Coupled_TendonStiffness,Misc);
