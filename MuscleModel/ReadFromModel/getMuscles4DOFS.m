function Misc=getMuscles4DOFS(Misc)
%   getMuscles4DOFS Selects all the muscles that actuate the DOFS specified by the user

for t = 1:Misc.nTrials
    if isempty(Misc.MuscleNames_Input{t})
        % Loop over each DOF in the model
        for i=1:length(Misc.DofNames_Input{t})
            % Get the moment arms from the muscle analysis results
            dm_Data_temp=ReadMotFile(fullfile(Misc.MuscleAnalysisPath,...
                [Misc.MAtrialName{t} '_MuscleAnalysis_MomentArm_' Misc.DofNames_Input{t}{i} '.sto']));    
            dM = dm_Data_temp.data(1,2:end);    
            dM_store(i,:)=dM;
        end

        % get the muscles that actuate the selected DOFS (moment arms > 0.001)
        MuscleNames=dm_Data_temp.names(2:end);
        [~, nm]=size(dM_store);
        ct=1;
        for i=1:nm
            if any(abs(dM_store(:,i))>0.0001)
                Misc.MuscleNames_Input{t}{ct}=MuscleNames{i};
                Misc.idx_allMuscleList_Input{t}(ct) = find(ismember(Misc.allMuscleList,MuscleNames{i}));
                ct=ct+1;
            end    
        end
    
        % print to screen
        disp(['MusclesNames Selected automatically for ' Misc.IKfile{t} ':']);
        disp(Misc.MuscleNames_Input{t}');
        disp(' ')
        clear dM_store
    end
    for i=1:length(Misc.MuscleNames_Input{t})
        Misc.idx_allMuscleList_Input{t}(i) = find(ismember(Misc.allMuscleList,Misc.MuscleNames_Input{t}{i}));
    end
end

end