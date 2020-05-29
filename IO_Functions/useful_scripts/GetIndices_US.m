function [DatStore] = GetIndices_US(DatStore,Misc,i)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if Misc.UStracking == 1
    DatStore(i).free_lMo = zeros(length(Misc.Estimate_OptFL),1);
    DatStore(i).USsel = zeros(length(Misc.USSelection),1);
    DatStore(i).free_kT = zeros(length(Misc.Estimate_TendonStifness),1);
    DatStore(i).coupled_kT = zeros(size(Misc.Coupled_TendonStifness))';
    DatStore(i).coupled_lMo = zeros(size(Misc.Coupled_fiber_length))';
    DatStore(i).coupled_lTs = zeros(size(Misc.Coupled_slack_length))';
    for j = 1:length(DatStore(i).free_lMo)
        DatStore(i).free_lMo(j) = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_OptFL{j}));
    end
    
    for j = 1:length(DatStore(i).free_kT)
        DatStore(i).free_kT(j) = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_TendonStifness{j}));
    end
    
    for j = 1:length(DatStore(i).USsel)
        DatStore(i).USsel(j) = find(strcmp(DatStore(i).MuscleNames,Misc.USSelection{j}));
    end
    
    for k = 1:size((DatStore(i).coupled_kT),1)
        for j = 1:size((DatStore(i).coupled_kT),2)
            DatStore(i).coupled_kT(k,j) = find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_TendonStifness{j,k}));
        end
    end
    % added coupling of muscle fiber length 24/01/2020
    for k = 1:size((DatStore(i).coupled_lMo),1)
        for j = 1:size((DatStore(i).coupled_lMo),2)
            DatStore(i).coupled_lMo(k,j) = find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_fiber_length{j,k}));
        end
    end
    % added coupling of tendon slack length 24/01/2020
    for k = 1:size((DatStore(i).coupled_lTs),1)
        for j = 1:size((DatStore(i).coupled_lTs),2)
            DatStore(i).coupled_lTs(k,j) = find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_slack_length{j,k}));
        end
    end
end
end

