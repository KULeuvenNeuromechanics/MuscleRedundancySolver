function [Misc,DatStore] = GetEMGInfo(Misc,DatStore)
%GetEMGInfo Reads the file with EMG information (.sto)., runs activation
%dynamics on the EMG data when asked and puts the EMG data in the right
%format (handles copies and so on).



% Author: Maarten Afschrift

if Misc.boolEMG
    for trial = Misc.trials_sel
        if isempty(Misc.EMGfile{trial}) || isempty(Misc.IKfile{trial}) || isempty(Misc.IDfile{trial})
            warning(['EMG or IK or ID file for trial ' Misc.trialName{trial} 'has not been defined. Please update Misc.EMGfile, Misc.IKfile or Misc.IDfile']);
        end
    end
end

if Misc.boolEMG    
    % file information
    bool_error  = 0;
    IndError=zeros(length(Misc.EMGSelection),length(Misc.trials_sel));
    % Load the data and check for errors
    EMGFile = struct;
    for t = 1:length(Misc.trials_sel)
        iF = Misc.trials_sel(t);
        % get information for the EMG constraints
        clear emgFile
        emgFile = read_motionFile_v40(Misc.EMGfile{iF}); 
        % prevent errors with the headers
        EMGFile(iF).colheaders_EMGfile = emgFile.labels;
        ct = 0;
        if isfield(Misc,'EMGFileHeaderCorrespondence')
            for c=1:length(EMGFile(iF).colheaders_EMGfile)
                idx_col_c1 = find(ismember(Misc.EMGFileHeaderCorrespondence(:,1),EMGFile(iF).colheaders_EMGfile{c}));
                if isempty(idx_col_c1)
                    warning(['Muscle corresponding to ' EMGFile(iF).colheaders_EMGfile{c} ' in file ' Misc.EMGfile{iF} ' is not defined in Misc.EMGFileHeaderCorrespondence. Update Misc.EMGFileHeaderCorrespondence. Removing EMG column ' EMGFile(iF).colheaders_EMGfile{c} ' from analyses'])
                else
                    ct = ct + 1;
                    EMGFile(iF).colheaders{ct} = Misc.EMGFileHeaderCorrespondence{idx_col_c1,2};
                    EMGFile(iF).data(:,ct) = emgFile.data(:,c);
                end
            end                
        end
        % check if we have to update the headers based on user input
        bool_updateheader   = 0;
        if ~isempty(Misc.EMGheaders{iF})        
            bool_updateheader=1;
        end
        % verify if the selected muscles are in the model
        % verify if the muscles in the .mot files are in the model
        % verify if the muscles in Misc.EMGSelection are in the .mot file 
        EMGheaders{iF}  = EMGFile(iF).colheaders;
        if bool_updateheader
            EMGheaders{iF} = Misc.EMGheaders{iF};
        else
            Misc.EMGheaders{iF} = EMGheaders{iF};
        end
        
        ct = 0;
        for i=1:length(Misc.EMGSelection)
            if any(ismember(DatStore(iF).MuscleNames,Misc.EMGSelection{i}))
                if any(ismember(EMGheaders{iF},Misc.EMGSelection{i}))
                    ct = ct +1;
                    % First column corresponds to index in allMuscleList
                    % Second column corresponds to index in Misc.EMGheaders
                    % Third column corresponds to index in MuscleNames of that trial                    
                    Misc.idx_EMGsel{iF}(ct,1) = find(ismember(Misc.allMuscleList,Misc.EMGSelection{i}));
                    Misc.idx_EMGsel{iF}(ct,2) = find(ismember(EMGheaders{iF},Misc.EMGSelection{i}));
                    Misc.idx_EMGsel{iF}(ct,3) = find(ismember(DatStore(iF).MuscleNames,Misc.EMGSelection{i}));
                    Misc.EMGsel{iF}{ct} = Misc.EMGSelection{i};
                else
                    if bool_updateheader == 0
                        disp(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Update the headers of file: ' Misc.EMGfile{iF}]);
                    else
                        disp(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Update the headers in:  Misc.EMGheaders{' num2str(iF) '}']);
                    end
                    IndError(i,t)=1;
                end
            else
                if strcmp(Misc.EMGSelection{i}(end),Misc.side{iF}) % All muscles end with 'l' or 'r'
                    disp(['Could not find ' Misc.EMGSelection{i} ' in the model for trial ' Misc.EMGfile{iF} ', Update the Misc.EMGSelection structure']);
                    IndError(i,t)=1;
                end
            end
        end
    end
    IndErr = zeros(length(Misc.EMGSelection),1);
    for i = 1:length(Misc.EMGSelection)
        if sum(IndError(i,:)) == length(Misc.trials_sel)
            bool_error = 1;
            IndErr(i) = 1;
        end
    end
    if bool_error ==1
        warning(['Removed muscles with EMG information from the',...
            ' analysis because these muscles are not in the model, or do not span the selected DOFs for any trial (see above)']);
        Misc.EMGSelection(find(IndError)) = [];
    end    
   
    %% Process the data 
    for iF = Misc.trials_sel
        EMGdat    = EMGFile(iF).data;        
        % get the EMG data
        EMGselection = Misc.EMGsel{iF};
        EMGsel = EMGdat(:,Misc.idx_EMGsel{iF}(:,2));
        
        % add twins
        if  ~isempty(Misc.EMG_MuscleCopies)            
            nCopy = size(Misc.EMG_MuscleCopies,1);
            for j=1:nCopy
                NameSel = Misc.EMG_MuscleCopies{j,1};
                NameCopy =  Misc.EMG_MuscleCopies{j,2};
                [EMGselection,EMGsel,Misc] = addEMGtwins(EMGselection,EMGsel,EMGdat,Misc,NameSel,NameCopy,DatStore,EMGheaders,iF);
                NameSel = Misc.EMG_MuscleCopies{j,2};
                NameCopy =  Misc.EMG_MuscleCopies{j,1};
                [EMGselection,EMGsel,Misc] = addEMGtwins(EMGselection,EMGsel,EMGdat,Misc,NameSel,NameCopy,DatStore,EMGheaders,iF);
            end
        end    
        DatStore(iF).EMG.BoundsScaleEMG = Misc.BoundsScaleEMG;
        DatStore(iF).EMG.EMGbounds      = Misc.EMGbounds;
        DatStore(iF).EMG.nEMG           = length(EMGselection);
%         DatStore(iF).EMG.EMGindices     = EMGindices;
        DatStore(iF).EMG.idx_EMGsel     = Misc.idx_EMGsel{iF};
        DatStore(iF).EMG.EMGsel         = EMGsel;
        DatStore(iF).EMG.EMGselection   = EMGselection;
        DatStore(iF).EMG.time           = EMGdat(:,1);
        DatStore(iF).EMG.boolEMG        = Misc.boolEMG;
        DatStore(iF).EMG.EMGspline      = spline(DatStore(iF).EMG.time',DatStore(iF).EMG.EMGsel');        
    end   
else
    for iF = Misc.trials_sel
        % Boolean in DatStore that EMG info is not used ?
        DatStore(iF).EMG.BoundsScaleEMG = [];
        DatStore(iF).EMG.EMGbounds      = [];
        DatStore(iF).EMG.nEMG           = [];
        DatStore(iF).EMG.EMGindices     = [];
        DatStore(iF).EMG.EMGsel         = [];
        DatStore(iF).EMG.EMGselection   = [];
        DatStore(iF).EMG.boolEMG        = Misc.boolEMG;
    end    
end

