function [DatStore,scaledMuscles] = getEMG_scaleIndecies(DatStore,ntrials)
% --------------------------------------------------------------------------
%getEMG_scaleIndecies
%     EMGscale of a given muscles is the same in all trials
%     Different trials can have different muscles with EMG
%     This function finds which muscles need to have EMG scaled
%     This fuction also assigns index of the EMG scale factors to the corresponding muscle in each trial
% 
% INPUT:
%     DatStore
%     Structure of all data
% 
%     ntrials
%     number of trials
% 
% OUTPUT:
%     DatStore
%     Structure of all data
% 
%     scaledMuscles
%     All the the muscles that need to have EMG scaled
% 
% Original author: Dhruv Gupta
% Original date: May 3, 2022
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

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
