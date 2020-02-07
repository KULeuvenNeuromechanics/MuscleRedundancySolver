function [DatStore] = GetEMGInfo(Misc,DatStore)
%GetEMGInfo Reads the file with EMG information (.sto)., runs activation
%dynamics on the EMG data when asked and puts the EMG data in the right
%format (handles copies and so on).



% Author: Maarten Afschrift

boolEMG = 0;
% check if the input is correct
if isfield(Misc,'EMGconstr') && Misc.EMGconstr == 1
    boolEMG = 1;
    nF = length(Misc.EMGfile);
    if nF ~=length(DatStore)
        disp('Warning: number of EMG files is not equal to the number of IK or ID files.');
    end
end


if boolEMG    
    % file information
    nF = length(Misc.EMGfile);
    % Load the data and check for errors
    for iF = 1:nF        
        % get information for the EMG constraints
        EMGFile(iF)      = importdata(Misc.EMGfile{iF});        
    end    
    % prevent errors with the headers
    for iF = 1:nF
        if ~isfield(EMGFile(iF),'colheaders')
            EMGFile(iF).colheaders = strsplit(EMGFile(1).textdata{end});
        end
    end
    % check if we have to update the headers based on user input
    bool_updateheader   = 0;
    if isfield(Misc,'EMGheaders') && ~isempty(Misc.EMGheaders)        
        bool_updateheader=1;
    end
    % verify if the selected muscles are in the model
    iF       = 1;    
    bool_error  = 0;
    IndError=zeros(length(Misc.EMGSelection),1);
    for i=1:length(Misc.EMGSelection)
        if ~any(strcmp(Misc.EMGSelection{i},DatStore(iF).MuscleNames))
            disp(['Could not find ' Misc.EMGSelection{i} ' in the model, Update the Misc.EMGSelection structure']);
            bool_error=1;
            IndError(i)=1;
        end
    end
    % verify if the muscles in the .mot files are in the model
    EMGheaders  = EMGFile(iF).colheaders;
    if bool_updateheader
       EMGheaders      = Misc.EMGheaders; 
    end
    for i=1:length(Misc.EMGSelection)
        if ~any(strcmp(Misc.EMGSelection{i},EMGheaders))
            if bool_updateheader == 0
                disp(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Updata the headers of file: ' Misc.EMGfile]);
            else
                disp(['Could not find ' Misc.EMGSelection{i} ' in the header of the EMG file, Update the headers in:  Misc.EMGheaders']);
            end
            bool_error=1;
            IndError(i)=1;
        end
    end
    if bool_error ==1
        warning(['Removed several muscles with EMG information from the',...
            ' analysis because these muscles are not in the model, or do not span the selected DOFs (see above)']);
        Misc.EMGSelection(find(IndError)) = [];
    end    
    
    %% Process the data 
    % ERROR IN FOR LOOP PARAMETERS FIXED - 24/01/2020 JPCB
    for iF = 1:nF
        EMGdat              = EMGFile(iF).data;        
        [nfr, nc] = size(EMGdat);  
        % get the EMG data
        nIn = length(Misc.EMGSelection);
        EMGsel = nan(nfr,nIn);   EMGindices = nan(nIn,1);
        EMGselection = Misc.EMGSelection;
        for i=1:length(Misc.EMGSelection)
            ind = strcmp(Misc.EMGSelection{i},EMGheaders);
            EMGsel(:,i) = EMGdat(:,ind);
            EMGindices(i) = find(strcmp(Misc.EMGSelection{i},DatStore(iF).MuscleNames));
        end        
        % add twins
        if  ~isempty(Misc.EMG_MuscleCopies)
            nCopy = length(Misc.EMG_MuscleCopies(:,1));
            EMGsel = [EMGsel zeros(nfr,nCopy)];
            EMGindices = [ EMGindices ; zeros(nCopy,1)];        
            for j=1:length(Misc.EMG_MuscleCopies(:,1))
                NameSel = Misc.EMG_MuscleCopies{j,1};
                NameCopy =  Misc.EMG_MuscleCopies{j,2};

                Ind_ColCopy = strcmp(Misc.EMGSelection,NameSel);
                EMGsel(:,i+j) = EMGsel(:,Ind_ColCopy);
                EMGindices(i+j) = find(strcmp(NameCopy,DatStore(iF).MuscleNames));
                EMGselection = [EMGselection {NameCopy}];
            end
        end    
        DatStore(iF).EMG.BoundsScaleEMG = Misc.BoundsScaleEMG;
        DatStore(iF).EMG.EMGbounds      = Misc.EMGbounds;
        DatStore(iF).EMG.nEMG           = length(EMGindices);
        DatStore(iF).EMG.EMGindices     = EMGindices;
        DatStore(iF).EMG.EMGsel         = EMGsel;
        DatStore(iF).EMG.EMGselection   = EMGselection;
        DatStore(iF).EMG.time           = EMGdat(:,1);
        DatStore(iF).EMG.boolEMG        = boolEMG;
        DatStore(iF).EMG.EMGspline      = spline(DatStore(iF).EMG.time',DatStore(iF).EMG.EMGsel');        
    end   
else
    for iF = 1:length(DatStore)
        % Boolean in DatStore that EMG info is not used ?
        DatStore(iF).EMG.BoundsScaleEMG       = [];
        DatStore(iF).EMG.EMGbounds      = [];
        DatStore(iF).EMG.nEMG           = [];
        DatStore(iF).EMG.EMGindices     = [];
        DatStore(iF).EMG.EMGsel         = [];
        DatStore(iF).EMG.EMGselection   = [];
        DatStore(iF).EMG.boolEMG         = boolEMG;
    end    
end

