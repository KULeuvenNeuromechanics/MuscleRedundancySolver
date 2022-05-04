function [Misc,DatStore] = GetUSInfo(Misc,DatStore)
% --------------------------------------------------------------------------
%GetEMGInfo
%    GetEMGInfo Reads the file with EMG information and puts the EMG data in the right
%    format.
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
% Original date:
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

if Misc.UStracking
    for trial = 1:Misc.nTrials
        if isempty(Misc.USfile{trial}) || isempty(Misc.IKfile{trial}) || isempty(Misc.IDfile{trial})
            warning(['US or IK or ID file for trial ' num2str(trial) ' has not been defined. Please update Misc.USfile, Misc.IKfile or Misc.IDfile']);
        end
    end
end
Misc.boolUS = Misc.UStracking;

if Misc.boolUS    
    % file information
    bool_error  = 0;
    IndError=zeros(length(Misc.USSelection),Misc.nTrials);
    % Load the data and check for errors
    USFile = struct;
    for iF = 1:Misc.nTrials      
        % get information for the US constraints
        clear usFile
        usFile = ReadMotFile(Misc.USfile{iF});         
        % prevent errors with the headers
        USFile(iF).colheaders_USfile = usFile.names';
        
        if isfield(Misc,'USheaders')
            USFile(iF).colheaders = Misc.USheaders{iF};
            USheaders{iF} = Misc.USheaders{iF};
            USFile(iF).data = usFile.data;
            bool_updateheader=1;
        else
            if isfield(Misc,'USFileHeaderCorrespondence')
                ct = 0;
                for c=1:length(USFile(iF).colheaders_USfile)
                    idx_col_c1 = find(ismember(Misc.USFileHeaderCorrespondence(:,1),USFile(iF).colheaders_USfile{c}));
                    if isempty(idx_col_c1)
                        warning(['Muscle corresponding to ' USFile(iF).colheaders_USfile{c} ' in file ' Misc.USfile{iF} ' is not defined in Misc.USFileHeaderCorrespondence. Update Misc.USFileHeaderCorrespondence. Removing US column ' USFile(iF).colheaders_USfile{c} ' from analyses'])
                    else
                        ct = ct + 1;
                        USFile(iF).colheaders{ct} = Misc.USFileHeaderCorrespondence{idx_col_c1,2};
                        USFile(iF).data(:,ct) = usFile.data(:,c);
                    end
                end
            else
                USFile(iF).colheaders = USFile(iF).colheaders_USfile;
                USFile(iF).data = usFile.data;
            end
            bool_updateheader=0;
            USheaders{iF} = USFile(iF).colheaders;
        end
        
        % verify if the selected muscles are in the model
        % verify if the muscles in the .mot files are in the model
        % verify if the muscles in Misc.USSelection are in the .mot file 
        
        ct = 0;
        for i=1:length(Misc.USSelection)
            if any(ismember(DatStore(iF).MuscleNames,Misc.USSelection{i}))
                if any(ismember(USheaders{iF},Misc.USSelection{i}))
                    ct = ct +1;
                    % First column corresponds to index in allMuscleList
                    % Second column corresponds to index in USheaders
                    % Third column corresponds to index in MuscleNames of that trial                    
                    Misc.idx_USsel{iF}(ct,1) = find(ismember(Misc.allMuscleList,Misc.USSelection{i}));
                    Misc.idx_USsel{iF}(ct,2) = find(ismember(USheaders{iF},Misc.USSelection{i}));
                    Misc.idx_USsel{iF}(ct,3) = find(ismember(DatStore(iF).MuscleNames,Misc.USSelection{i}));
                    Misc.USsel{iF}{ct} = Misc.USSelection{i};
                else
                    if bool_updateheader == 0
                        disp(['Could not find ' Misc.USSelection{i} ' in the header of the US file, Update the headers of file: ' Misc.USfile{iF}]);
                    else
                        disp(['Could not find ' Misc.USSelection{i} ' in the header of the US file, Update the headers in:  Misc.USheaders{' num2str(iF) '}']);
                    end
                    IndError(i,iF)=1;
                end
            else
                if strcmp(Misc.USSelection{i}(end),Misc.side{iF}) % All muscles end with 'l' or 'r'
                    disp(['Could not find ' Misc.USSelection{i} ' in the model for trial ' Misc.USfile{iF} ', Update the Misc.USSelection structure']);
                    IndError(i,iF)=1;
                end
            end
        end
    end
    IndErr = zeros(length(Misc.USSelection),1);
    for i = 1:length(Misc.USSelection)
        if sum(IndError(i,:)) == Misc.nTrials
            bool_error = 1;
            IndErr(i) = 1;
        end
    end
    if bool_error
        warning(['Removed muscles with US information from the',...
            ' analysis because these muscles are not in the model, or do not span the selected DOFs for any trial (see above)']);
        Misc.USSelection(find(IndErr)) = [];
    end
    
    %% Process the data
    for iF = 1:Misc.nTrials
        USdat              = USFile(iF).data;        
        % get the US data
        USsel = USdat(:,Misc.idx_USsel{iF}(:,2));
        DatStore(iF).US.nUS           = length(Misc.USsel{iF});
        DatStore(iF).US.idx_USsel     = Misc.idx_USsel{iF};
        DatStore(iF).US.USsel         = USsel;
        DatStore(iF).US.USselection   = Misc.USsel{iF};
        DatStore(iF).US.time          = USdat(:,1);
        DatStore(iF).US.boolUS        = Misc.boolUS;
        DatStore(iF).US.USspline      = spline(DatStore(iF).US.time',DatStore(iF).US.USsel');        
    end   
else
    for iF = 1:Misc.nTrials
        % Boolean in DatStore that US info is not used ?       
        DatStore(iF).US.nUS           = [];
        DatStore(iF).US.USindices     = [];
        DatStore(iF).US.USsel         = [];
        DatStore(iF).US.USselection   = [];
        DatStore(iF).US.boolUS        = Misc.boolUS;
    end    
end

