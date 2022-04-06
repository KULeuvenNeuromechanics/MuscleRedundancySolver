function [Misc,DatStore] = GetEMGInfo(Misc,DatStore)
%GetEMGInfo Reads the file with EMG information (.sto)., runs activation
%dynamics on the EMG data when asked and puts the EMG data in the right
%format (handles copies and so on).

if Misc.boolEMG
    for trial = 1:Misc.nTrials
        if isempty(Misc.EMGfile{trial}) || isempty(Misc.IKfile{trial}) || isempty(Misc.IDfile{trial})
            warning(['EMG or IK or ID file for trial ' Misc.trialName{trial} 'has not been defined. Please update Misc.EMGfile, Misc.IKfile or Misc.IDfile']);
        end
    end
end

if Misc.boolEMG    
    % file information
    bool_error  = 0;
    IndError=zeros(length(Misc.EMGSelection),Misc.nTrials);
    % Load the data and check for errors
    EMGFile = struct;
    for iF = 1:Misc.nTrials
        % get information for the EMG constraints
        clear emgFile
        EMGFile = ReadMotFile(Misc.EMGfile{iF});         
        % check if we have to update the headers based on user input
        bool_updateheader   = 0;
        if ~isempty(Misc.EMGheaders{iF})        
            bool_updateheader=1;
        end
        % verify if the selected muscles are in the model
        % verify if the muscles in the .mot files are in the model
        % verify if the muscles in Misc.EMGSelection are in the .mot file 
        EMGheaders{iF}  = EMGFile(iF).names;
        if bool_updateheader
            EMGheaders{iF} = Misc.EMGheaders;
        end
        
        ct = 0;
        EMGSel = Misc.EMGSelection{iF};
        for i=1:length(EMGSel)
            if any(ismember(DatStore(iF).MuscleNames,EMGSel{i}))
                if any(ismember(EMGheaders{iF},EMGSel{i}))
                    ct = ct +1;
                    % First column corresponds to index in allMuscleList
                    % Second column corresponds to index in Misc.EMGheaders
                    % Third column corresponds to index in MuscleNames of that trial                    
                    Misc.idx_EMGsel{iF}(ct,1) = find(ismember(Misc.allMuscleList,EMGSel{i}));
                    Misc.idx_EMGsel{iF}(ct,2) = find(ismember(EMGheaders{iF},EMGSel{i}));
                    Misc.idx_EMGsel{iF}(ct,3) = find(ismember(DatStore(iF).MuscleNames,EMGSel{i}));
                    Misc.EMGsel{iF}{ct} = EMGSel{i};
                else
                    if bool_updateheader == 0
                        disp(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Update the headers of file: ' Misc.EMGfile{iF}]);
                    else
                        disp(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Update the headers in:  Misc.EMGheaders{' num2str(iF) '}']);
                    end
                    IndError(i,iF)=1;
                end
            end
        end
    end
    IndErr = zeros(length(Misc.EMGSelection),1);
    for i = 1:length(Misc.EMGSelection)
        if sum(IndError(i,:)) == Misc.nTrials
            bool_error = 1;
            IndErr(i) = 1;
        end
    end
    if bool_error ==1
        warning(['Removed muscles with EMG information from the',...
            ' analysis because these muscles are not in the model, or do not span the selected DOFs for any trial (see above)']);
        Misc.EMGSelection(find(IndErr)) = [];
    end    
   
    %% Process the data 
    for iF = 1:Misc.nTrials
        EMGdat    = EMGFile(iF).data;        
        % get the EMG data
        EMGselection = Misc.EMGSelection{iF};
        EMGsel = EMGdat(:,Misc.idx_EMGsel{iF}(:,2));
        
        % add twins
        if  ~isempty(Misc.EMG_MuscleCopies)            
            nCopy = size(Misc.EMG_MuscleCopies,1);
            for j=1:nCopy
                NameSel = Misc.EMG_MuscleCopies{j,1};
                NameCopy =  Misc.EMG_MuscleCopies{j,2};
                [EMGselection,EMGsel,Misc] = addEMGtwins(EMGselection,EMGsel,EMGdat,...
                    Misc,NameSel,NameCopy,DatStore,EMGheaders,iF);
                NameSel = Misc.EMG_MuscleCopies{j,2};
                NameCopy =  Misc.EMG_MuscleCopies{j,1};
                [EMGselection,EMGsel,Misc] = addEMGtwins(EMGselection,EMGsel,EMGdat,...
                    Misc,NameSel,NameCopy,DatStore,EMGheaders,iF);
            end
        end    
        DatStore(iF).EMG.BoundsScaleEMG = Misc.BoundsScaleEMG; % bounds on scale factors
        DatStore(iF).EMG.EMGbounds      = Misc.EMGbounds; % bounds on EMG tracking
        DatStore(iF).EMG.nEMG           = length(EMGselection); % number of EMG signals
        DatStore(iF).EMG.idx_EMGsel     = Misc.idx_EMGsel{iF}; % indices of EMG in modelled muscles
        DatStore(iF).EMG.EMGsel         = EMGsel; % selected EMG data from .mot file
        %DatStore(iF).EMG.EMGselection   = EMGselection;
        DatStore(iF).EMG.time           = EMGdat(:,1);
        DatStore(iF).EMG.boolEMG        = Misc.boolEMG; % boolean if EMG data is used
        %DatStore(iF).EMG.EMGspline      = spline(DatStore(iF).EMG.time',DatStore(iF).EMG.EMGsel');        
    end   
else
    for iF = 1:Misc.nTrials
        % Boolean in DatStore that EMG info is not used ?
        DatStore(iF).EMG.BoundsScaleEMG = []; % bounds on scale factors
        DatStore(iF).EMG.EMGbounds      = []; % bounds on EMG tracking
        DatStore(iF).EMG.nEMG           = []; % number of EMG signals
        DatStore(iF).EMG.idx_EMGsel     = []; % indices of EMG in modelled muscles
        DatStore(iF).EMG.EMGsel         = []; % selected EMG data from .mot file
        %DatStore(iF).EMG.EMGselection   = [];
        DatStore(iF).EMG.boolEMG        = Misc.boolEMG; % boolean if EMG data is used
    end    
end

