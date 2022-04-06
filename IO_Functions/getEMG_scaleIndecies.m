function [DatStore,scaledMuscles] = getEMG_scaleIndecies(DatStore,ntrials)
mNamesEMG = [];
for t=1:ntrials
    mNamesEMG = [mNamesEMG DatStore(t).EMG.EMGselection];
end
scaledMuscles = unique(mNamesEMG);

for t=1:ntrials
    for m = 1:DatStore(t).EMG.nEMG
        DatStore(t).EMG.idx_EMGsel(m,4) = find(ismember(scaledMuscles,DatStore(t).EMG.EMGselection{m}));
    end
end
