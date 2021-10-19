function [DatStore] = GetUSInfo(Misc,DatStore)
%GetUSInfo Reads the file with US information (.sto)., stores data in the right
%format.



% Author: Maarten Afschrift

boolUS = 0;
% check if the input is correct
if isfield(Misc,'UStracking') && Misc.UStracking == 1
    boolUS = 1;
    nF = length(Misc.USfile);
    if nF ~=length(DatStore)
        disp('Warning: number of US files is not equal to the number of IK or ID files.');
    end
end

if boolUS 
    
    % file information
    nFiles = length(Misc.USfile);
    % Load the data and check for errors
    for iFile = 1:nF        
        % get information for the EMG constraints
        USfile(iFile)      = importdata(Misc.USfile{iFile});        
    end    
    % prevent errors with the headers
    for iFile = 1:nF
%         if ~isfield(USfile(iFile),'colheaders')
            USfile(iFile).colheaders = strsplit(USfile(iFile).textdata{end});
%         end
    end
    % verify if the selected muscles are in the model
    for iFile = 1:nF 
        bool_error  = 0;
        USselection = Misc.USSelection(iFile,:);
        IndError=zeros(length(USselection),1);
        for i=1:length(Misc.USSelection)
            if ~any(strcmp(Misc.USSelection{i},DatStore(iFile).MuscleNames))
                disp(['Could not find ' USselection{i} ' in the model, update the Misc.USSelection structure']);
                bool_error=1;
                IndError(i)=1;
            end
        end
        % verify if the muscles in the .mot files are in the model
        USheaders  = USfile(iFile).colheaders;
        
        for i=1:length(USselection)
            if ~any(strcmp(USselection{i},USheaders))
                if bool_updateheader == 0
                    disp(['Could not find ' USselection{i} ' in the header of the US file, Updata the headers of file: ' Misc.USfile(iFile)]);
                else
                    disp(['Could not find ' USselection{i} ' in the header of the US file, Update the headers in:  Misc.USheaders']);
                end
                bool_error=1;
                IndError(i)=1;
            end
        end
        if bool_error ==1
            warning(['Removed several muscles with US information from the',...
                ' analysis because these muscles are not in the model, or do not span the selected DOFs (see above)']);
            Misc.USSelection(find(IndError)) = [];
        end    
    end
    
    %% Process the data    
    for iF = 1:nF
        USdat              = USfile(iF).data;        
        [nfr, nc] = size(USdat);  
        % get the US data
        nIn = length(Misc.USSelection(iF,:));
        USsel = nan(nfr,nIn);   USindices = nan(nIn,1);
        USselection = Misc.USSelection(iF,:);
        USheaders  = USfile(iF).colheaders;
        for i=1:length(USselection)
            ind = strcmp(USselection{i},USheaders);
            USsel(:,i) = USdat(:,ind);
            USindices(i) = find(strcmp(USselection{i},DatStore(iF).MuscleNames));
        end 
        DatStore(iF).US.nUS           = length(USindices);
        DatStore(iF).US.USindices     = USindices;
        DatStore(iF).US.USsel         = USsel;
        DatStore(iF).US.USselection   = USselection;
        DatStore(iF).US.time           = USdat(:,1);
        DatStore(iF).US.boolUS         = boolUS;
        DatStore(iF).US.USspline      = spline(DatStore(iF).US.time',DatStore(iF).US.USsel');        
    end   
else
    for iF = 1:length(DatStore)
        % Boolean in DatStore that US info is not used ?       
        DatStore(iF).US.nUS           = [];
        DatStore(iF).US.USindices     = [];
        DatStore(iF).US.USsel         = [];
        DatStore(iF).US.USselection   = [];
        DatStore(iF).US.boolUS         = boolUS;
    end    
end

