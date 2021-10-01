function [DatStore] = GetIndices_OptMuscleParams(DatStore,Misc)
% Add to DatStore the information on the indices of which muscles parameters will be
% estimated.

% loop over all trials. This is a bit redundant since these indices are the
% same for all trials
if Misc.UStracking == 1 || Misc.EMGconstr == 1
    
   % pre allocate indices
    free_lMo = zeros(length(Misc.Estimate_OptimalFiberLength),1);
    free_kT = zeros(length(Misc.Estimate_TendonStiffness),1);
    coupled_kT = zeros(size(Misc.Coupled_TendonStiffness));
    coupled_lMo = zeros(size(Misc.Coupled_fiber_length));
    coupled_lTs = zeros(size(Misc.Coupled_slack_length));
    
    % muscles with free fiber lengths
    for j = 1:length(free_lMo)
        IndsSel = find(strcmp(DatStore(1).MuscleNames,Misc.Estimate_OptimalFiberLength{j}));
        if ~isempty(IndsSel)
            free_lMo(j) = IndsSel;
        else
            warning(['Cannot optimize optimal fiber length of muscle ' Misc.Estimate_OptimalFiberLength{j},...
                ' because this muscle is not included in the model']);
        end
    end
    % muscles with free tendon length
    for j = 1:length(free_kT)
        IndsSel = find(strcmp(DatStore(1).MuscleNames,Misc.Estimate_TendonStiffness{j}));
        if ~isempty(IndsSel)
            free_kT(j) = IndsSel;
        else
            warning(['Cannot optimize tendon stiffness of muscle ' Misc.Estimate_TendonStiffness{j},...
                ' because this muscle is not included in the model']);
        end
    end
    % couples scale factor for optimized tendon lengths
    for k = 1:size(coupled_kT,1)
        for j = 1:size(coupled_kT,2)
            Inds = find(strcmp(DatStore(1).MuscleNames,Misc.Coupled_TendonStiffness{k,j}));
            if ~isempty(Inds)
                coupled_kT(k,j) = Inds;
            else
                warning(['Could not couple tendon stiffness of muscle ' Misc.Coupled_TendonStiffness{k,j},...
                    'because this muscle in not included in the model']);
            end
        end
    end
    % added coupling of muscle fiber length
    for k = 1:size(coupled_lMo,1)
        for j = 1:size(coupled_lMo,2)
            Inds= find(strcmp(DatStore(1).MuscleNames,Misc.Coupled_fiber_length{k,j}));
            if ~isempty(Inds)
                coupled_lMo(k,j) = Inds;
            else
                warning(['Could not couple fiber length of muscle ' Misc.Coupled_TendonStiffness{k,j},...
                    'because this muscle in not included in the model']);
            end
        end
    end
    % added coupling of tendon slack length
    for k = 1:size(coupled_lTs,1)
        for j = 1:size(coupled_lTs,2)
            Inds= find(strcmp(DatStore(1).MuscleNames,Misc.Coupled_slack_length{k,j}));
            if ~isempty(Inds)
                coupled_lTs(k,j) = Inds;
            else
                warning(['Could not couple tendon slack length of muscle ' Misc.Coupled_TendonStiffness{k,j},...
                    'because this muscle in not included in the model']);
            end
        end
    end
end

% add to DatStore
for i=1:length(DatStore)
    DatStore(i).free_lMo = free_lMo;
    DatStore(i).free_kT = free_kT;
    DatStore(i).coupled_kT = coupled_kT;
    DatStore(i).coupled_lMo = coupled_lMo;
    DatStore(i).coupled_lTs = coupled_lTs;
end
end

