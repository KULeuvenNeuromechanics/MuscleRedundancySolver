function [Misc,DatStore] = GetEMGInfo(Misc,DatStore)
% --------------------------------------------------------------------------
%GetEMGInfo
%    GetEMGInfo Reads the file with EMG information and puts the EMG data in the right
%    format (handles copies and so on).
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
% Original author: Maarten Afschrift
% Original date: November 26, 2019
%
% Updated by: Bram Van Den Bosch
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

if Misc.boolEMG
    for trial = 1:Misc.nTrials
        if isempty(Misc.EMGfile{trial}) || isempty(Misc.IKfile{trial}) || isempty(Misc.IDfile{trial})
            warning(['EMG or IK or ID file for trial ' num2str(trial) ' has not been defined. Please update Misc.EMGfile, Misc.IKfile or Misc.IDfile']);
            disp(' ')
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
        emgFile = ReadMotFile(Misc.EMGfile{iF});         
        % prevent errors with the headers
        EMGFile(iF).colheaders_EMGfile = emgFile.names';
        
        if isfield(Misc,'EMGheaders')
            EMGFile(iF).colheaders = Misc.EMGheaders{iF};
            EMGheaders{iF} = Misc.EMGheaders{iF};
            EMGFile(iF).data = emgFile.data;
            bool_updateheader=1;
        else
            if isfield(Misc,'EMGFileHeaderCorrespondence')
                ct = 0;
                for c=1:length(EMGFile(iF).colheaders_EMGfile)
                    idx_col_c1 = find(ismember(Misc.EMGFileHeaderCorrespondence(:,1),EMGFile(iF).colheaders_EMGfile{c}));
                    if isempty(idx_col_c1)
                        warning(['Muscle corresponding to ' EMGFile(iF).colheaders_EMGfile{c} ' in file ' Misc.EMGfile{iF} ' is not defined in Misc.EMGFileHeaderCorrespondence. Update Misc.EMGFileHeaderCorrespondence. Removing EMG column ' EMGFile(iF).colheaders_EMGfile{c} ' from analyses'])
                        disp(' ')
                    else
                        ct = ct + 1;
                        EMGFile(iF).colheaders{ct} = Misc.EMGFileHeaderCorrespondence{idx_col_c1,2};
                        EMGFile(iF).data(:,ct) = emgFile.data(:,c);
                    end
                end
            else
                EMGFile(iF).colheaders = EMGFile(iF).colheaders_EMGfile;
                EMGFile(iF).data = emgFile.data;
            end
            bool_updateheader=0;
            EMGheaders{iF} = EMGFile(iF).colheaders;
        end
        
        % verify if the selected muscles are in the model
        % verify if the muscles in the .mot files are in the model
        % verify if the muscles in Misc.EMGSelection are in the .mot file 
        
        ct = 0;
        for i=1:length(Misc.EMGSelection)
            if any(ismember(DatStore(iF).MuscleNames,Misc.EMGSelection{i}))
                if any(ismember(EMGheaders{iF},Misc.EMGSelection{i}))
                    ct = ct +1;
                    % First column corresponds to index in allMuscleList
                    % Second column corresponds to index in EMGheaders
                    % Third column corresponds to index in MuscleNames of that trial                    
                    Misc.idx_EMGsel{iF}(ct,1) = find(ismember(Misc.allMuscleList,Misc.EMGSelection{i}));
                    Misc.idx_EMGsel{iF}(ct,2) = find(ismember(EMGheaders{iF},Misc.EMGSelection{i}));
                    Misc.idx_EMGsel{iF}(ct,3) = find(ismember(DatStore(iF).MuscleNames,Misc.EMGSelection{i}));
                    Misc.EMGsel{iF}{ct} = Misc.EMGSelection{i};
                else
                    if bool_updateheader == 0
                        warning(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Update the headers of file: ' Misc.EMGfile{iF}]);
                        disp(' ')
                    else
                        warning(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Update the headers in:  Misc.EMGheaders{' num2str(iF) '}']);
                        disp(' ')
                    end
                    IndError(i,iF)=1;
                end
            else
                if strcmp(Misc.EMGSelection{i}(end),Misc.side{iF}) % All muscles end with 'l' or 'r'
                    warning(['Could not find ' Misc.EMGSelection{i} ' in the model for trial ' Misc.EMGfile{iF} ', Update the Misc.EMGSelection structure']);
                    disp(' ')
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
    if bool_error
        warning(['Removed muscles with EMG information from the',...
            ' analysis because these muscles are not in the model, or do not span the selected DOFs for any trial (see above)']);
        disp(' ')
        Misc.EMGSelection(find(IndErr)) = [];
    end    
   
    %% Process the data 
    for iF = 1:Misc.nTrials
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
        DatStore(iF).EMG.EMGselection   = EMGselection;
        DatStore(iF).EMG.time           = EMGdat(:,1);
        DatStore(iF).EMG.boolEMG        = Misc.boolEMG; % boolean if EMG data is used
        DatStore(iF).EMG.EMGspline      = spline(DatStore(iF).EMG.time',DatStore(iF).EMG.EMGsel');        
    end   
else
    for iF = 1:Misc.nTrials
        % Boolean in DatStore that EMG info is not used ?
        DatStore(iF).EMG.BoundsScaleEMG = []; % bounds on scale factors
        DatStore(iF).EMG.EMGbounds      = []; % bounds on EMG tracking
        DatStore(iF).EMG.nEMG           = []; % number of EMG signals
        DatStore(iF).EMG.idx_EMGsel     = []; % indices of EMG in modelled muscles
        DatStore(iF).EMG.EMGsel         = []; % selected EMG data from .mot file
        DatStore(iF).EMG.EMGselection   = [];
        DatStore(iF).EMG.boolEMG        = Misc.boolEMG; % boolean if EMG data is used
    end    
end

