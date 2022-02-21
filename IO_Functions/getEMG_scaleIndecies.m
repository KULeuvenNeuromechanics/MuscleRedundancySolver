function [DatStore,scaledMuscles] = getEMG_scaleIndecies(DatStore,trials_sel)
mNamesEMG = [];
for t=trials_sel
    mNamesEMG = [mNamesEMG DatStore(t).EMG.EMGselection];
end
scaledMuscles = unique(mNamesEMG);

for t=trials_sel
    for m = 1:DatStore(t).EMG.nEMG
        DatStore(t).EMG.idx_EMGsel(m,4) = find(ismember(scaledMuscles,DatStore(t).EMG.EMGselection{m}));
    end
end
