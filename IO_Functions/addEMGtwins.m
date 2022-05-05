function [EMGselection,EMGsel,Misc] = addEMGtwins(EMGselection,EMGsel,EMGdat,Misc,NameSel,NameCopy,DatStore,EMGheaders,iF)
% --------------------------------------------------------------------------
%addEMGtwins
%     This function adds the EMG twin (muscles for which same EMG information is to be used). 
% 
% INPUT:
%     EMGselection
%     list of all muscles with EMG in that trial
%     
%     EMGsel
%     EMG data of the muscles with EMG in that trial
%     
%     EMGdat
%     Data of EMG channels in that trial
%     
%     Misc
%     Miscellaneous info used through the code
%     
%     NameSel
%     Muscle whose EMG is to be copied
%     
%     NameCopy
%     Muscle that will use NameSel's EMG
%     
%     DatStore
%     Structure of all data
%     
%     EMGheaders
%     Names of all muscles with EMG data
%     
%     iF
%     current trial number  
% 
% OUTPUT:
%     EMGselection
%     list of all muscles with EMG in that trial
%     
%     EMGsel
%     EMG data of the muscles with EMG in that trial
%     
%     Misc
%     Miscellaneous info used through the code
% 
% Original author: Dhruv Gupta
% Original date: May 3, 2022
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

% check if EMG signals we want to copy exists
Ind_ColCopy = strcmp(Misc.EMGsel{iF},NameSel);
% check if twin muscle is in the model
Ind_ColOut = strcmp(DatStore(iF).MuscleNames,NameCopy);
% check if twin muscle has EMG data
BoolTwinHasEMG = strcmp(EMGheaders{iF},NameCopy);
% check if twin muscle has already been added
TwinAlreadyAdded = strcmp(EMGselection,NameCopy);
if any(Ind_ColCopy) && any(Ind_ColOut) && ~any(BoolTwinHasEMG) && ~any(TwinAlreadyAdded)% only if both muscles are selected
    Misc.idx_EMGsel{iF}(end+1,1) = find(ismember(Misc.allMuscleList,NameCopy));
    Misc.idx_EMGsel{iF}(end,2) = find(ismember(EMGheaders{iF},NameSel));
    Misc.idx_EMGsel{iF}(end,3) = find(ismember(DatStore(iF).MuscleNames,NameCopy));
    Misc.EMGsel{iF}{end+1} = NameCopy;
    EMGsel = [EMGsel EMGdat(:,Misc.idx_EMGsel{iF}(end,2))];
    EMGselection = [EMGselection {NameCopy}];
elseif ~any(Ind_ColOut)
    if strcmp(NameCopy(end),Misc.side{iF}) % All muscles end with 'l' or 'r'
        disp([' Cannot copy EMG muscle ' NameSel ' to twin ',...
            NameCopy ' because ' NameCopy ' is not selected in the model',...
            'for file ' Misc.EMGfile{iF}]);
        disp(' ')
    end
elseif any(Ind_ColOut) && any(BoolTwinHasEMG) && ~any(TwinAlreadyAdded)
    disp(['twin muscle ' NameCopy ' has EMG data as input, we therefore did', ... 
        ' not constrain this activity based on ' NameSel,...
        'for file ' Misc.EMGfile{iF},...
        '. We based it on the activity of ' NameCopy]);
    disp(' ')
    Misc.idx_EMGsel{iF}(end+1,1) = find(ismember(Misc.allMuscleList,NameCopy));
    Misc.idx_EMGsel{iF}(end,2) = find(ismember(EMGheaders{iF},NameCopy));
    Misc.idx_EMGsel{iF}(end,3) = find(ismember(DatStore(iF).MuscleNames,NameCopy));
    Misc.EMGsel{iF}{end+1} = NameCopy;
    EMGsel = [EMGsel EMGdat(:,Misc.idx_EMGsel{iF}(end,2))];
    EMGselection = [EMGselection {NameCopy}];
end
end
