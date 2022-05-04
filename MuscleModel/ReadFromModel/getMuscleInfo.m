function [Misc,DatStore] = getMuscleInfo(Misc,DatStore)
%   Get_dof_MuscleInfo selects the DOF that are acuated by muscles specified by the user and selects for those dof the moment arms of the muscles
%   author: Maarten Afschrift,
%   Last Update: 29 Januari 2019

%% Read muscle analysis
for t = 1:Misc.nTrials
    trial = t;
    IK_path = Misc.IKfile{t};
    ID_path = Misc.IDfile{t};
    
    % Pre-allocate loop variables
    DOF_inds=nan(length(Misc.DofNames_Input{t}),1);
    ct=1;

    % Loop over each DOF in the model
    for i=1:length(Misc.DofNames_Input{t})
        % read the Muscle Analysis Result
        MA_FileName=fullfile(Misc.MuscleAnalysisPath,[Misc.MAtrialName{t} '_MuscleAnalysis_MomentArm_' Misc.DofNames_Input{t}{i} '.sto']);
        if exist(MA_FileName,'file')
            dm_Data_temp=ReadMotFile(fullfile(Misc.MuscleAnalysisPath,...
                [Misc.MAtrialName{t} '_MuscleAnalysis_MomentArm_' Misc.DofNames_Input{t}{i} '.sto']));
        else
            error(['Cannot open muscle analysis results for: ' Misc.DofNames_Input{t}{i}])
        end

        % get the indexes for the selected MuscleNames (only needed in first iteration)
        if i==1
            nfr = length(dm_Data_temp.data(:,1));
            headers=dm_Data_temp.names;
            Inds_muscles=nan(length(Misc.MuscleNames_Input{t}),1);
            IndsNames_sel=nan(length(Misc.MuscleNames_Input{t}),1);
            ctm=1;
            for j=1:length(Misc.MuscleNames_Input{t})
                ind_sel=find(strcmp(Misc.MuscleNames_Input{t}{j},headers));
                if ~isempty(ind_sel)
                    Inds_muscles(ctm)=ind_sel; IndsNames_sel(ctm)=j;
                    ctm=ctm+1;
                else
                    warning(['The selected muscle ' Misc.MuscleNames_Input{t}{j} ' does not exist in the selected model of trial ' num2str(t) '. This muscles is removed from the program']);
                    disp(' ')
                end
            end
            Misc.MuscleNames{t}=Misc.MuscleNames_Input{t}(IndsNames_sel);
            Misc.idx_allMuscleList{t}=Misc.idx_allMuscleList_Input{t}(IndsNames_sel);
            Inds_muscles(isnan(Inds_muscles))=[];                               % Delete the muscles names that are not selected by the user
            dM_temp=nan(nfr,length(Misc.DofNames_Input{t}),length(Misc.MuscleNames{t}));    % pre-allocate moment arms
        end

        % Evaluate if one of the muscles spans this DOF (when moment arms > 0.001)
        dM=dm_Data_temp.data(:,Inds_muscles);    
        if any(any(abs(dM)>0.0001))
            Misc.DofNames_muscles{t}{ct}=Misc.DofNames_Input{t}{i};
            dM_temp(:,i,:)=dM;
            DOF_inds(ct)=i;
            ct=ct+1;   
        end     
    end

    % Combine DOFs_actuated by muscles and the DOFS selected by the user
    ct=1;
    Inds_deleteDOFS=zeros(length(Misc.DofNames_muscles{t}),1);
    for i=1:length(Misc.DofNames_muscles{t})
        if ~any(strcmp(Misc.DofNames_Input{t},Misc.DofNames_muscles{t}{i}))
             Inds_deleteDOFS(ct)=i;ct=ct+1;
        end
    end
    Inds_deleteDOFS(ct:end)=[];
    DOF_inds(Inds_deleteDOFS)=[];
    Misc.DofNames{t}=Misc.DofNames_muscles{t}; Misc.DofNames{t}(Inds_deleteDOFS)=[];

    % warnings when not all the input DOFS are actuated by muscles
    for i=1:length(Misc.DofNames_Input{t})
        if ~any(strcmp(Misc.DofNames_Input{t}{i},Misc.DofNames{t}))
            warning(['The input dof: ' Misc.DofNames_Input{t}{i} ' is not actuated by the selected muscles and therefore removed from the analysis']);
            disp(' ')
        end
    end

    % Filter the moment arms information and store them in DatStore.dM
    dM_raw=dM_temp(:,DOF_inds,:);
    t_dM = dm_Data_temp.data(:,1);
    fs=1/mean(diff(t_dM));
    [B,A] = butter(Misc.f_order_dM, Misc.f_cutoff_dM/(fs/2));
    DatStore(trial).dM = filtfilt(B,A,dM_raw);

    % filter Muscle-tendon lengths and store them in DatStore.LMT
    LMT_dat=ReadMotFile(fullfile(Misc.MuscleAnalysisPath,[Misc.MAtrialName{t} '_MuscleAnalysis_Length.sto']));
    LMT_raw=LMT_dat.data(:,Inds_muscles);
    t_lMT = LMT_dat.data(:,1);
    fs=1/mean(diff(t_lMT));             % sampling frequency
    [B,A] = butter(Misc.f_order_lMT,Misc.f_cutoff_lMT/(fs/2));
    DatStore(trial).LMT = filtfilt(B,A,LMT_raw);

    % store information in the DatStore structure
    DatStore(trial).MuscleNames = Misc.MuscleNames{t};
    DatStore(trial).DOFNames    = Misc.DofNames{t};
    DatStore(trial).nMuscles    = length(Misc.MuscleNames{t});
    DatStore(trial).nDOF        = length(Misc.DofNames{t});

    %% Filter IK

    % load joint kinematics
    IK_data=ReadMotFile(IK_path);

    % select the IK information between the selected time frames
    t_IK=IK_data.data(:,1);
    ind0=find(t_IK>=Misc.time(trial,1),1,'first'); ind_end=find(t_IK<=Misc.time(trial,2),1,'last');
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
    ID_data=ReadMotFile(ID_path);
    t_ID=ID_data.data(:,1);

    % get the ID index
    ID_header=ID_data.names;     IK_header = strtrim(IK_data.names);
    ID_Header_inds=zeros(size(DOF_inds));  IK_Header_inds = zeros(size(DOF_inds));
    for i=1:length(Misc.DofNames{t})
       ID_Header_inds(i)=find(strcmp([Misc.DofNames{t}{i} '_moment'],ID_header));
       IK_Header_inds(i)=find(strcmp(Misc.DofNames{t}{i},IK_header)); 
    end

    % filter the ID data and store in Datstore.T_exp
    fs=1/mean(diff(t_ID));
    [B, A] = butter(Misc.f_order_ID, Misc.f_cutoff_ID/(fs/2));
    [~,M] = size(ID_data.data);
    for i = 2:M         % filtering time vector not needed (start from 2)
        ID_data.data(:,i) = filtfilt(B,A,ID_data.data(:,i));
    end

    % select ID data between start and end
    ID_data_int = interp1(ID_data.data(:,1),ID_data.data,IK_data.data(:,1));       % interpolate data for IK sampling frequency
    t_ID = IK_data.data(:,1);
    ind0 = find(t_ID>=Misc.time(trial,1),1,'first'); ind_end=find(t_ID<=Misc.time(trial,2),1,'last');
    ID_inds=ind0:ind_end;

    %% store the data
    DatStore(trial).T_exp = ID_data_int(ID_inds,ID_Header_inds);
    DatStore(trial).q_exp = IK_data.data(IK_inds,IK_Header_inds); 
    DatStore(trial).time = t_IK(IK_inds);


    % check if size of IK and ID matrices are equal
    if length(ID_inds) ~= length(IK_inds)    
        disp(['Time range IK in the solution file: first time frame ' num2str(IK_data.data(1,1)) '  last time frame:' num2str(IK_data.data(end,1))]);
        disp(['Time range ID in the solution file: first time frame ' num2str(ID_data.data(1,1)) '  last time frame:' num2str(ID_data.data(end,1))]);
        disp(['Selected time range is: ' num2str(Misc.time(trial,:))]);
        error('There is something wrong with the time frames in your IK or ID file');
    end
end