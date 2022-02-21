function [Misc] = Warnings_MuscleNames(DatStore,Misc)
% Generate warnings should there be inconsistensies in the user-provided
% muscle names.

[Misc] = issueWarning1D(Misc.Estimate_OptimalFiberLength,'Estimate_OptimalFiberLength',DatStore,Misc);
[Misc] = issueWarning1D(Misc.Estimate_TendonSlackLength,'Estimate_TendonSlackLength',DatStore,Misc);
[Misc] = issueWarning1D(Misc.Estimate_TendonStiffness,'Estimate_TendonStiffness',DatStore,Misc);
[Misc] = issueWarning1D(Misc.EMGSelection,'EMGSelection',DatStore,Misc);

[Misc] = issueWarning2D(Misc.Coupled_fiber_length,'Coupled_fiber_length',DatStore,Misc);
[Misc] = issueWarning2D(Misc.Coupled_TendonStiffness,'Coupled_TendonStiffness',DatStore,Misc);
[Misc] = issueWarning2D(Misc.Coupled_slack_length,'Coupled_slack_length',DatStore,Misc);
[Misc] = issueWarning2D(Misc.EMG_MuscleCopies,'EMG_MuscleCopies',DatStore,Misc);

end

