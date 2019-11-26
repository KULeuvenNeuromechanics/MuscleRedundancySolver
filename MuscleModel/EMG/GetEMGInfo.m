function [DatStore] = GetEMGInfo(Misc,DatStore)
%GetEMGInfo Reads the file with EMG information (.sto)., runs activation
%dynamics on the EMG data when asked and puts the EMG data in the right
%format (handles copies and so on).

% Author: Maarten Afschrift

%   Detailed explanation goes here
% get information for the EMG constraints
EMGFile = importdata(Misc.EMGfile);
EMGdat = EMGFile.data;
EMGheaders = EMGFile.colheaders;
bool_updateheader=0;
if ~isempty(Misc.EMGheaders)
    EMGheaders = Misc.EMGheaders;
    bool_updateheader=1;
end
[nfr, nc]= size(EMGdat);

% process the EMG through activation dynamics

if isfield(Misc,'ActDynEMG') && Misc.ActDynEMG == 1
    disp('Running Activation dynamics on EMG data');
    % integrator settings
    AbsTol = 1e-5;      RelTol = 1e-3;
    options = odeset('AbsTol',AbsTol,'RelTol',RelTol);
    % pre-allcation
    actEMG = nan(nfr,nc);
    time = EMGdat(:,1);
    actEMG(:,1) = time;  % copy the time column
    tsim = [time(1) time(end)];
    if isfield(Misc,'PlotBool') && Misc.PlotBool ==1
        figure();
    end    
    for i=2:nc        
        e = EMGdat(:,i);
        x0 = mean(e(1:3));      % initial state based on excitations in first frames
        [t,x] = ode15s(@ActivationOde,tsim,x0, options, time, e);
        actEMG(:,i) = interp1(t,x,time);        
        if isfield(Misc,'PlotBool') && Misc.PlotBool ==1
            subplot(4,ceil((nc-1)./4),i-1);
            plot(time,e,'r'); hold on;
            plot(t,x,'--b');
            if i == nc
                legend('EMG','actEMG');
                title(EMGheaders{i});
            end
        end        
    end    
    EMGdat = actEMG;    
end
if isfield(Misc,'PlotBool') && Misc.PlotBool ==1
    suptitle('Activation dynamics');
end

% safety check
bool_error = 0;
IndError=zeros(length(Misc.EMGSelection),1);
for i=1:length(Misc.EMGSelection)
    if ~any(strcmp(Misc.EMGSelection{i},DatStore.MuscleNames))
        disp(['Could not find ' Misc.EMGSelection{i} ' in the model, Update the Misc.EMGSelection structure']);
        bool_error=1;
        IndError(i)=1;
    end
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
    warning(['Removed several muscles with EMG informaion from the',...
        ' analysis because these muscles are not in the model, or do not span the selected DOFs (see above)']);
    Misc.EMGSelection(find(IndError)) = [];
end

% get the EMG data
nIn = length(Misc.EMGSelection);
EMGsel = nan(nfr,nIn);   EMGindices = nan(nIn,1);
EMGselection = Misc.EMGSelection;
for i=1:length(Misc.EMGSelection)
    ind = strcmp(Misc.EMGSelection{i},EMGheaders);
    EMGsel(:,i) = EMGdat(:,ind);
    EMGindices(i) = find(strcmp(Misc.EMGSelection{i},DatStore.MuscleNames));
end

% add twins
nCopy = length(Misc.EMG_MuscleCopies(:,1));
EMGsel = [EMGsel zeros(nfr,nCopy)];
EMGindices = [ EMGindices ; zeros(nCopy,1)];

for j=1:length(Misc.EMG_MuscleCopies(:,1))
    NameSel = Misc.EMG_MuscleCopies{j,1};
    NameCopy =  Misc.EMG_MuscleCopies{j,2};
    
    Ind_ColCopy = strcmp(Misc.EMGSelection,NameSel);
    EMGsel(:,i+j) = EMGsel(:,Ind_ColCopy);
    EMGindices(i+j) = find(strcmp(NameCopy,DatStore.MuscleNames));
    EMGselection = [EMGselection {NameCopy}];
end

EMGint = interp1(EMGdat(:,1),EMGsel,DatStore.time);


DatStore.MaxScale       = Misc.MaxScale;
DatStore.EMGbounds      = Misc.EMGbounds;
DatStore.nEMG           = length(EMGindices);
DatStore.EMGindices     = EMGindices;
DatStore.actEMGintSO    = EMGint;
DatStore.EMGselection   = EMGselection;
end

