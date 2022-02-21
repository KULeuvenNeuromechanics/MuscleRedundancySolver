function [Misc,DatStore] = GetUSInfo(Misc,DatStore)
%GetUSInfo Reads the file with US information (.sto)., stores data in the right
%format.



% Author: Maarten Afschrift

if Misc.UStracking
    for trial = Misc.trials_sel
        if isempty(Misc.USfile{trial}) || isempty(Misc.IKfile{trial}) || isempty(Misc.IDfile{trial})
            warning(['US or IK or ID file for trial ' Misc.trialName{trial} 'has not been defined. Please update Misc.USfile, Misc.IKfile or Misc.IDfile']);
        end
    end
end
Misc.boolUS = Misc.UStracking;

if Misc.boolUS    
    % file information
    nFiles = length(Misc.USfile);
    bool_error  = 0;
    IndError=zeros(length(Misc.USSelection),1);
    % Load the data and check for errors
    for iF = Misc.trials_sel      
        % get information for the US constraints
        USfile(iF)      = importdata(Misc.USfile{iF});        
        % prevent errors with the headers
        if ~isfield(USfile(iF),'colheaders')
            USfile(iF).colheaders = strsplit(USfile(iF).textdata{end});
        end
        USheaders{iF}  = USfile(iF).colheaders;
        % verify if the selected muscles are in the model
        % verify if the muscles in the .mot files are in the model
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
                    disp(['Could not find ' Misc.USSelection{i} ' in the header of the US file, Update the headers of file: ' Misc.USfile{iF}]);
                    bool_error=1;
                    IndError(i)=1;
                end
            else
                disp(['Could not find ' Misc.USSelection{i} ' in the model for trial ' Misc.USfile{iF} ', update the Misc.USSelection structure']);
                bool_error=1;
                IndError(i)=1;
            end
        end
    end
    if bool_error ==1
        warning(['Removed several muscles with US information from the',...
            ' analysis because these muscles are not in the model, or do not span the selected DOFs (see above)']);
        Misc.USSelection(find(IndError)) = [];
    end    
    %% Process the data    
    for iF = Misc.trials_sel
        USdat              = USfile(iF).data;        
        % get the US data
        USsel = USdat(:,Misc.idx_USsel{iF}(:,2));
        DatStore(iF).US.nUS           = length(Misc.USsel{iF});
%         DatStore(iF).US.USindices     = USindices;
        DatStore(iF).US.idx_USsel     = Misc.idx_USsel{iF};
        DatStore(iF).US.USsel         = USsel;
        DatStore(iF).US.USselection   = Misc.USsel{iF};
        DatStore(iF).US.time          = USdat(:,1);
        DatStore(iF).US.boolUS        = Misc.boolUS;
        DatStore(iF).US.USspline      = spline(DatStore(iF).US.time',DatStore(iF).US.USsel');        
    end   
else
    for iF = Misc.trials_sel
        % Boolean in DatStore that US info is not used ?       
        DatStore(iF).US.nUS           = [];
        DatStore(iF).US.USindices     = [];
        DatStore(iF).US.USsel         = [];
        DatStore(iF).US.USselection   = [];
        DatStore(iF).US.boolUS        = Misc.boolUS;
    end    
end

