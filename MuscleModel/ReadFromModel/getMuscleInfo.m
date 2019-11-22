function [DatStore] = getMuscleInfo(IK_path,ID_path,Misc)
%   Get_dof_MuscleInfo selects the DOF that are acuated by muscles specified by the user and selects for those dof the moment arms of the muscles
%   author: Maarten Afschrift,
%   Last Update: 29 Januari 2019

%% Read muscle analysis
% Pre-allocate loop variables
DOF_inds=nan(length(Misc.DofNames_Input),1);
ct=1;

% Loop over each DOF in the model
for i=1:length(Misc.DofNames_Input)
    
    % read the Muscle Analysis Result
    MA_FileName=fullfile(Misc.MuscleAnalysisPath,[Misc.trialName '_MuscleAnalysis_MomentArm_' Misc.DofNames_Input{i} '.sto']);
    if exist(MA_FileName,'file')
        dm_Data_temp=importdata(fullfile(Misc.MuscleAnalysisPath,[Misc.trialName '_MuscleAnalysis_MomentArm_' Misc.DofNames_Input{i} '.sto']));
    else
        error(['Cannot open muscle analysis results for: ' Misc.DofNames_Input{i}])
    end
    
    % get the indexes for the selected MuscleNames (only needed in first iteration)
    if i==1    
        nfr = length(dm_Data_temp.data(:,1));
        headers=dm_Data_temp.colheaders;
        Inds_muscles=nan(length(Misc.MuscleNames_Input),1);        
        ctm=1;
        for j=1:length(Misc.MuscleNames_Input)
            ind_sel=find(strcmp(Misc.MuscleNames_Input{j},headers));
            if ~isempty(ind_sel)
                Inds_muscles(ctm)=ind_sel; IndsNames_sel(ctm)=j;
                ctm=ctm+1;
            else
                disp(['Warning: The selected muscle ' Misc.MuscleNames_Input{j} ' does not exist in the selected model. This muscles is removed from the program']);                
            end
        end
        Misc.MuscleNames=Misc.MuscleNames_Input(IndsNames_sel);
        Inds_muscles(isnan(Inds_muscles))=[];                               % Delete the muscles names that are not selected by the user
        dM_temp=nan(nfr,length(Misc.DofNames_Input),length(Misc.MuscleNames));    % pre-allocate moment arms
        
        % read indexes in time frame for muscle analysis
        t_Mus=dm_Data_temp.data(:,1);           t_Mus=round(t_Mus*10000)/10000;
        ind0=find(t_Mus>=Misc.time(1),1,'first'); ind_end=find(t_Mus<=Misc.time(2),1,'last');
        Mus_inds=ind0:ind_end;
    end
    
    % Evaluate if one of the muscles spans this DOF (when moment arms > 0.001)
    dM=dm_Data_temp.data(Mus_inds,Inds_muscles);    
    if any(any(abs(dM)>0.001))
        Misc.DofNames_muscles{ct}=Misc.DofNames_Input{i};
        dM_temp(:,i,:)=dM;
        DOF_inds(ct)=i;
        ct=ct+1;   
    end     
end

% Combine DOFs_actuated by muscles and the DOFS selected by the user
ct=1;
Inds_deleteDOFS=zeros(length(Misc.DofNames_muscles),1);
for i=1:length(Misc.DofNames_muscles)
    if ~any(strcmp(Misc.DofNames_Input,Misc.DofNames_muscles{i}))
         Inds_deleteDOFS(ct)=i;ct=ct+1;
    end
end
Inds_deleteDOFS(ct:end)=[];
DOF_inds(Inds_deleteDOFS)=[];
Misc.DofNames=Misc.DofNames_muscles; Misc.DofNames(Inds_deleteDOFS)=[];

% warnings when not all the input DOFS are actuated by muscles
for i=1:length(Misc.DofNames_Input)
    if ~any(strcmp(Misc.DofNames_Input{i},Misc.DofNames))
        disp(['Warning DOF: The input dof: ' Misc.DofNames_Input{i} ' is not actuated by the selected muscles and therefore removed from the analysis']);
    end
end

% Filter the moment arms information and store them in DatStore.dM
dM_raw=dM_temp(:,DOF_inds,:);
t_dM = dm_Data_temp.data(:,1);
t_dM=round(t_dM*10000)/10000;
fs=1/mean(diff(t_dM));
[B,A] = butter(Misc.f_order_dM, Misc.f_cutoff_dM/(fs/2));
DatStore.dM = filtfilt(B,A,dM_raw);

% filter Muscle-tendon lengths and store them in DatStore.LMT
LMT_dat=importdata(fullfile(Misc.MuscleAnalysisPath,[Misc.trialName '_MuscleAnalysis_Length.sto']));
LMT_raw=LMT_dat.data(Mus_inds,Inds_muscles);
t_lMT = LMT_dat.data(:,1);
t_lMT=round(t_lMT*10000)/10000;
fs=1/mean(diff(t_lMT));             % sampling frequency
[B,A] = butter(Misc.f_order_lMT,Misc.f_cutoff_lMT/(fs/2));
DatStore.LMT = filtfilt(B,A,LMT_raw);

% store information in the DatStore structure
DatStore.MuscleNames = Misc.MuscleNames;
DatStore.DOFNames    = Misc.DofNames;
DatStore.nMuscles    = length(Misc.MuscleNames);
DatStore.nDOF        = length(Misc.DofNames);

%% Filter IK

% load joint kinematics
[~,Misc.trialName,~]=fileparts(IK_path);
IK_data=importdata(IK_path);
if ~isfield(IK_data,'colheaders')
    IK_data.colheaders=strsplit(IK_data.textdata{end});
end

% select the IK information between the selected time frames
t_IK=IK_data.data(:,1);     t_IK=round(t_IK*10000)/10000; 
ind0=find(t_IK>=Misc.time(1),1,'first'); ind_end=find(t_IK<=Misc.time(2),1,'last');
IK_inds=ind0:ind_end;

% filter the kinematics and kinetics
fs=1/mean(diff(t_IK));
[B, A] = butter(Misc.f_order_IK, Misc.f_cutoff_IK/(fs/2));
[~,DOFS] = size(IK_data.data);
for i = 2:DOFS         % filtering time vector not needed (start from 2)
    IK_data.data(:,i) = filtfilt(B,A,IK_data.data(:,i));
end

%% Filter ID
% Get the ID data
ID_data=importdata(ID_path);
t_ID=ID_data.data(:,1); t_ID=round(t_ID*10000)/10000;

% get the ID index
if ~isfield(ID_data,'colheaders')
    ID_data.colheaders=strsplit(ID_data.textdata{end});
end
ID_header=ID_data.colheaders;     IK_header = IK_data.colheaders;
ID_Header_inds=zeros(size(DOF_inds));  IK_Header_inds = zeros(size(DOF_inds));
for i=1:length(Misc.DofNames)
   ID_Header_inds(i)=find(strcmp([Misc.DofNames{i} '_moment'],ID_header));
   IK_Header_inds(i)=find(strcmp(Misc.DofNames{i},IK_header)); 
end

% filter the ID data and store in Datstore.T_exp
fs=1/mean(diff(t_ID));
[B, A] = butter(Misc.f_order_ID, Misc.f_cutoff_ID/(fs/2));
[~,M] = size(ID_data.data);
for i = 2:M         % filtering time vector not needed (start from 2)
    ID_data.data(:,i) = filtfilt(B,A,ID_data.data(:,i));
end

% select ID data between start and end
ID_data_int=interp1(ID_data.data(:,1),ID_data.data,IK_data.data(:,1));       % interpolate data for IK sampling frequency
t_ID=ID_data_int(:,1); t_ID=round(t_ID*10000)/10000;
ind0=find(t_ID>=Misc.time(1),1,'first'); ind_end=find(t_ID<=Misc.time(2),1,'last');
ID_inds=ind0:ind_end;

%% store the data
DatStore.T_exp = ID_data_int(ID_inds,ID_Header_inds);
DatStore.q_exp = IK_data.data(IK_inds,IK_Header_inds); 
DatStore.time = t_IK(IK_inds);


% check if size of IK and ID matrices are equal
if length(ID_inds) ~= length(IK_inds)    
    disp(['Time range IK in the solution file: first time frame ' num2str(IK_data.data(1,1)) '  last time frame:' num2str(IK_data.data(end,1))]);
    disp(['Time range ID in the solution file: first time frame ' num2str(ID_data.data(1,1)) '  last time frame:' num2str(ID_data.data(end,1))]);
    disp(['Selected time range is: ' num2str(Misc.time)]);
    error('There is something wrong with the time frames in your IK or ID file');
end
