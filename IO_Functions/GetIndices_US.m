function [DatStore] = GetIndices_US(DatStore,Misc,i)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if Misc.UStracking == 1 || Misc.EMGconstr == 1 
    DatStore(i).free_lMo = zeros(length(Misc.Estimate_OptimalFiberLength),1);
    DatStore(i).free_kT = zeros(length(Misc.Estimate_TendonStiffness),1);
    DatStore(i).coupled_kT = zeros(size(Misc.Coupled_TendonStiffness))';
    DatStore(i).coupled_lMo = zeros(size(Misc.Coupled_fiber_length))';
    DatStore(i).coupled_lTs = zeros(size(Misc.Coupled_slack_length))';
    
    for j = 1:length(DatStore(i).free_lMo)
        IndsSel = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_OptimalFiberLength{j}));
        if ~isempty(IndsSel)
            DatStore(i).free_lMo(j) = IndsSel;
        else
            warning(['Cannot optimize optimal fiber length of muscle ' Misc.Estimate_OptimalFiberLength{j},...
                ' because this muscle is not included in the model']);
        end
    end    
    for j = 1:length(DatStore(i).free_kT)
        IndsSel = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_TendonStiffness{j}));
        if ~isempty(IndsSel)
            DatStore(i).free_kT(j) = IndsSel;
        else
            warning(['Cannot optimize tendon stiffness of muscle ' Misc.Estimate_TendonStiffness{j},...
                ' because this muscle is not included in the model']);
        end
    end    
    if Misc.UStracking
        DatStore(i).USsel = zeros(length(Misc.USSelection),1);        
        for j = 1:length(DatStore(i).USsel)
            DatStore(i).USsel(j) = find(strcmp(DatStore(i).MuscleNames,Misc.USSelection{j}));
        end
    end
    for k = 1:size((DatStore(i).coupled_kT),1)
        for j = 1:size((DatStore(i).coupled_kT),2)
            Inds = find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_TendonStiffness{j,k}));
            if ~isempty(Inds)
                DatStore(i).coupled_kT(k,j) = Inds;
            else
                warning(['Could not couple tendon stiffness of muscle ' Misc.Coupled_TendonStiffness{j,k},...
                'because this muscle in not included in the model']);
            end
        end
    end
    % added coupling of muscle fiber length
    for k = 1:size((DatStore(i).coupled_lMo),1)
        for j = 1:size((DatStore(i).coupled_lMo),2)
             Inds= find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_fiber_length{j,k}));
            if ~isempty(Inds)
                DatStore(i).coupled_lMo(k,j) = Inds;
            else
                warning(['Could not couple fiber length of muscle ' Misc.Coupled_TendonStiffness{j,k},...
                'because this muscle in not included in the model']);
            end
        end
    end
    % added coupling of tendon slack length
    for k = 1:size((DatStore(i).coupled_lTs),1)
        for j = 1:size((DatStore(i).coupled_lTs),2)
             Inds= find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_slack_length{j,k}));
            if ~isempty(Inds)
                DatStore(i).coupled_lTs(k,j) = Inds;
            else
                warning(['Could not couple tendon slack length of muscle ' Misc.Coupled_TendonStiffness{j,k},...
                'because this muscle in not included in the model']);
            end
        end
    end
end
end

