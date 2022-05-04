function [DatStore,Misc] = GetIndices(DatStore,Misc)
% --------------------------------------------------------------------------
%GetIndecies
%     Add to DatStore the information on the indices of which muscles parameters will be
%     estimated.
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
%     DatStore
%     Structure of all data
%     
% Original author:
% Original date:
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

% if Misc.UStracking == 1 || Misc.EMGconstr == 1
Misc.free_lMo = zeros(length(Misc.Estimate_OptimalFiberLength),1);
Misc.free_lTs = zeros(length(Misc.Estimate_TendonSlackLength),1);
Misc.free_kT = zeros(length(Misc.Estimate_TendonStiffness),1);
Misc.coupled_kT = zeros(size(Misc.Coupled_TendonStiffness));
Misc.coupled_lMo = zeros(size(Misc.Coupled_fiber_length));
Misc.coupled_lTs = zeros(size(Misc.Coupled_slack_length));

for j = 1:length(Misc.free_lMo)
    IndsSel = find(strcmp(Misc.allMuscleList,Misc.Estimate_OptimalFiberLength{j}));
    if ~isempty(IndsSel)
        Misc.free_lMo(j) = IndsSel;
    end
end
for j = 1:length(Misc.free_lTs)
    IndsSel = find(strcmp(Misc.allMuscleList,Misc.Estimate_TendonSlackLength{j}));
    if ~isempty(IndsSel)
        Misc.free_lTs(j) = IndsSel;
    end
end
for j = 1:length(Misc.free_kT)
    IndsSel = find(strcmp(Misc.allMuscleList,Misc.Estimate_TendonStiffness{j}));
    if ~isempty(IndsSel)
        Misc.free_kT(j) = IndsSel;
    end
end
for k = 1:size(Misc.coupled_kT,1)
    for j = 1:size(Misc.coupled_kT,2)
        Inds = find(strcmp(Misc.allMuscleList,Misc.Coupled_TendonStiffness{k,j}));
        if ~isempty(Inds)
            Misc.coupled_kT(k,j) = Inds;
        end
    end
end
% added coupling of muscle fiber length
for k = 1:size(Misc.coupled_lMo,1)
    for j = 1:size(Misc.coupled_lMo,2)
        Inds= find(strcmp(Misc.allMuscleList,Misc.Coupled_fiber_length{k,j}));
        if ~isempty(Inds)
            Misc.coupled_lMo(k,j) = Inds;
        end
    end
end
% added coupling of tendon slack length
for k = 1:size((Misc.coupled_lTs),1)
    for j = 1:size((Misc.coupled_lTs),2)
        Inds= find(strcmp(Misc.allMuscleList,Misc.Coupled_slack_length{k,j}));
        if ~isempty(Inds)
            Misc.coupled_lTs(k,j) = Inds;
        end
    end
end

for i=1:Misc.nTrials
    DatStore(i).free_lMo = zeros(length(Misc.Estimate_OptimalFiberLength),3);
    DatStore(i).free_lTs = zeros(length(Misc.Estimate_TendonSlackLength),3);
    DatStore(i).free_kT = zeros(length(Misc.Estimate_TendonStiffness),3);
    DatStore(i).coupled_kT = zeros([size(Misc.Coupled_TendonStiffness) 3]);
    DatStore(i).coupled_lMo = zeros([size(Misc.Coupled_fiber_length) 3]);
    DatStore(i).coupled_lTs = zeros([size(Misc.Coupled_slack_length) 3]);
    
    for j = 1:size(DatStore(i).free_lMo,1)
        IndsSel = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_OptimalFiberLength{j}));
        if ~isempty(IndsSel)
            DatStore(i).free_lMo(j,3) = IndsSel;
            DatStore(i).free_lMo(j,1) = find(strcmp(Misc.allMuscleList,Misc.Estimate_OptimalFiberLength{j}));
        end
    end
    for j = 1:size(DatStore(i).free_lTs,1)
        IndsSel = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_TendonSlackLength{j}));
        if ~isempty(IndsSel)
            DatStore(i).free_lTs(j,3) = IndsSel;
            DatStore(i).free_lTs(j,1) = find(strcmp(Misc.allMuscleList,Misc.Estimate_TendonSlackLength{j}));
        end
    end
    for j = 1:size(DatStore(i).free_kT,1)
        IndsSel = find(strcmp(DatStore(i).MuscleNames,Misc.Estimate_TendonStiffness{j}));
        if ~isempty(IndsSel)
            DatStore(i).free_kT(j,3) = IndsSel;
            DatStore(i).free_kT(j,1) = find(strcmp(Misc.allMuscleList,Misc.Estimate_TendonStiffness{j}));
        end
    end
    for k = 1:size(DatStore(i).coupled_kT,1)
        for j = 1:size(DatStore(i).coupled_kT,2)
            Inds = find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_TendonStiffness{k,j}));
            if ~isempty(Inds)
                DatStore(i).coupled_kT(k,j,3) = Inds;
                DatStore(i).coupled_kT(k,j,1) = find(strcmp(Misc.allMuscleList,Misc.Coupled_TendonStiffness{k,j}));
            end
        end
    end
    % added coupling of muscle fiber length
    for k = 1:size(DatStore(i).coupled_lMo,1)
        for j = 1:size(DatStore(i).coupled_lMo,2)
            Inds= find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_fiber_length{k,j}));
            if ~isempty(Inds)
                DatStore(i).coupled_lMo(k,j,3) = Inds;
                DatStore(i).coupled_lMo(k,j,1) = find(strcmp(Misc.allMuscleList,Misc.Coupled_fiber_length{k,j}));
            end
        end
    end
    % added coupling of tendon slack length
    for k = 1:size(DatStore(i).coupled_lTs,1)
        for j = 1:size(DatStore(i).coupled_lTs,2)
            Inds= find(strcmp(DatStore(i).MuscleNames,Misc.Coupled_slack_length{k,j}));
            if ~isempty(Inds)
                DatStore(i).coupled_lTs(k,j,3) = Inds;
                DatStore(i).coupled_lTs(k,j,1) = find(strcmp(Misc.allMuscleList,Misc.Coupled_slack_length{k,j}));
            end
        end
    end
end
% end
end

