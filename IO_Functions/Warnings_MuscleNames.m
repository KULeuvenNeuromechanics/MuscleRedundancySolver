function [] = Warnings_MuscleNames(DatStore,Misc,i)
% Generate warnings should there be inconsistensies in the user-provided
% muscle names.

% get muscle names
MNames = DatStore(i).MuscleNames;

% optimization parameters muscles
[r,c] = size(Misc.Estimate_OptimalFiberLength);
for i = 1:r
    for j = 1:c
        if ~any(strcmp(MNames,Misc.Estimate_OptimalFiberLength{i,j}))
            warning(['Could not find muscle ' Misc.Estimate_OptimalFiberLength{i,j} ' in the selected dofs of the model',...
                ', please adapt Misc.Estimate_OptimalFiberLength']);
        end
    end
end
[r,c] = size(Misc.Estimate_TendonStiffness);
for i = 1:r
    for j = 1:c
        if ~any(strcmp(MNames,Misc.Estimate_TendonStiffness{i,j}))
            warning(['Could not find muscle ' Misc.Estimate_TendonStiffness{i,j} ' in the selected dofs of the model',...
                ', please adapt Misc.Estimate_TendonStiffness']);
        end
    end
end

[r,c] = size(Misc.Coupled_fiber_length);
for i = 1:r
    for j = 1:c
        if ~any(strcmp(MNames,Misc.Coupled_fiber_length{i,j}))
            warning(['Could not find muscle ' Misc.Coupled_fiber_length{i,j} ' in the selected dofs of the model',...
                ', please adapt Misc.Estimate_TendonStiffness']);
        end
    end
end

[r,c] = size(Misc.Coupled_slack_length);
for i = 1:r
    for j = 1:c
        if ~any(strcmp(MNames,Misc.Coupled_slack_length{i,j}))
            warning(['Could not find muscle ' Misc.Coupled_slack_length{i,j} ' in the selected dofs of the model',...
                ', please adapt Misc.Estimate_TendonStiffness']);
        end
    end
end

[r,c] = size(Misc.Coupled_TendonStiffness);
for i = 1:r
    for j = 1:c
        if ~any(strcmp(MNames,Misc.Coupled_TendonStiffness{i,j}))
            warning(['Could not find muscle ' Misc.Coupled_TendonStiffness{i,j} ' in the selected dofs of the model',...
                ', please adapt Misc.Estimate_TendonStiffness']);
        end
    end
end



% EMG information
[r,c] = size(Misc.EMG_MuscleCopies);
for i = 1:r
    for j = 1:c
        if ~any(strcmp(MNames,Misc.EMG_MuscleCopies{i,j}))
            warning(['Could not find muscle ' Misc.EMG_MuscleCopies{i,j} ' in the selected dofs of the model',...
                ', please adapt Misc.EMG_MuscleCopies']);
        end
    end
end
[r,c] = size(Misc.EMGSelection);
for i = 1:r
    for j = 1:c
        if ~any(strcmp(MNames,Misc.EMGSelection{i,j}))
            warning(['Could not find muscle ' Misc.EMGSelection{i,j} ' in the selected dofs of the model',...
                ', please adapt Misc.EMG_MuscleCopies']);
        end
    end
end



end

