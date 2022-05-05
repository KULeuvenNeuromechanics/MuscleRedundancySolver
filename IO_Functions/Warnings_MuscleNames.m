function [Misc] = Warnings_MuscleNames(DatStore,Misc)
% --------------------------------------------------------------------------
%Warnings_MuscleNames
%     Generate warnings should there be inconsistensies in the user-provided
%     muscle names.
% 
% INPUT:
%     Misc
%     Miscellaneous info used through the code
% 
%     DatStore
%     Structure of all data
%     
% OUTPUT:
%     Misc
%     Miscellaneous info used through the code
% 
% Original author: MaartenAfschrift
% Original date: June 8, 2020
%
% Updated by: Tom Van Wouwe
% 
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

% 1D matrices (list of muscles)
[Misc] = issueWarning1D(Misc.Estimate_OptimalFiberLength,'Estimate_OptimalFiberLength',DatStore,Misc);
[Misc] = issueWarning1D(Misc.Estimate_TendonSlackLength,'Estimate_TendonSlackLength',DatStore,Misc);
[Misc] = issueWarning1D(Misc.Estimate_TendonStiffness,'Estimate_TendonStiffness',DatStore,Misc);
[Misc] = issueWarning1D(Misc.EMGSelection,'EMGSelection',DatStore,Misc);

% 2D matrices (n x 2 matricies, where n is number of couples)
[Misc] = issueWarning2D(Misc.Coupled_fiber_length,'Coupled_fiber_length',DatStore,Misc);
[Misc] = issueWarning2D(Misc.Coupled_TendonStiffness,'Coupled_TendonStiffness',DatStore,Misc);
[Misc] = issueWarning2D(Misc.Coupled_slack_length,'Coupled_slack_length',DatStore,Misc);
[Misc] = issueWarning2D(Misc.EMG_MuscleCopies,'EMG_MuscleCopies',DatStore,Misc);

end

