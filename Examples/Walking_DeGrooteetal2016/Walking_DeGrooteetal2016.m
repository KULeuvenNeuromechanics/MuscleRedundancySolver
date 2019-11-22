clear all;close all;clc

%% Choose optimization framework
% framework = 'GPOPS';
framework = 'CasADi';

%% Choose contraction dynamics formulation
% formulation_contdyn = 'lMtildeState';
formulation_contdyn = 'FtildeState';

%% Choose activation dynamics formulation
% formulation_actdyn = 'DeGroote2016';
formulation_actdyn = 'DeGroote2009';

%% Example
% Add main folder and subfolder to matlab path (installation)
filepath=which('Walking_DeGrooteetal2016.m');
[DirExample_Walking,~,~]=fileparts(filepath); 
[DirExample,~]=fileparts(DirExample_Walking);
[MainDir,~]=fileparts(DirExample);
addpath(genpath(MainDir));

% Needed Input Arguments
IK_path=fullfile(MainDir,'Examples','Walking_DeGrooteetal2016','WalkingData','Walking_IK.mot');
ID_path=fullfile(MainDir,'Examples','Walking_DeGrooteetal2016','WalkingData','Walking_ID.sto');
model_path=fullfile(MainDir,'Examples','Walking_DeGrooteetal2016','WalkingData','subject1.osim');
time=[0.516 1.95]; % Right stance phase (+50ms beginning and end of time interval, more details see manual and publication)
Out_path=fullfile(MainDir,'Examples','Walking_DeGrooteetal2016','Results');

Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_rotation_r','hip_adduction_r'};
Misc.MuscleNames_Input={}; % Selects all muscles for input DOFs when empty

% Optional Input Arguments
Misc.f_cutoff_ID = 6;     % cutoff frequency filtering ID
Misc.f_order_ID = 4;      % order frequency filtering ID
Misc.f_cutoff_lMT = 6;    % cutoff frequency filtering lMT
Misc.f_order_lMT = 4;     % order frequency filtering lMT
Misc.f_cutoff_dM= 6;      % cutoff frequency filtering MA
Misc.f_order_dM = 4;      % order frequency filtering MA
Misc.f_cutoff_IK= 6;      % cutoff frequency filtering IK
Misc.f_order_IK = 4;      % order frequency filtering IK

%% Solve the problem
switch framework
    case 'GPOPS'
        switch formulation_actdyn
            case 'DeGroote2016'     
                switch formulation_contdyn
                    case 'lMtildeState'
                        [Time,MExcitation,MActivation,RActivation,TForcetilde,...
                            TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=...
                            SolveMuscleRedundancy_lMtildeState_GPOPS(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                    case 'FtildeState'   
                        [Time,MExcitation,MActivation,RActivation,TForcetilde,...
                            TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=...
                            SolveMuscleRedundancy_FtildeState_GPOPS(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                end

            case 'DeGroote2009' 
                switch formulation_contdyn
                    case 'lMtildeState'
                        [Time_actdyn,MExcitation_actdyn,MActivation_actdyn,...
                            RActivation_actdyn,TForcetilde_actdyn,TForce_actdyn,...
                            lMtilde_actdyn,lM_actdyn,MuscleNames_actdyn,...
                            OptInfo_actdyn,DatStore_actdyn]=...
                            SolveMuscleRedundancy_lMtildeState_actdyn_GPOPS(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                    case 'FtildeState'   
                        [Time_actdyn,MExcitation_actdyn,MActivation_actdyn,...
                            RActivation_actdyn,TForcetilde_actdyn,TForce_actdyn,...
                            lMtilde_actdyn,lM_actdyn,MuscleNames_actdyn,...
                            OptInfo_actdyn,DatStore_actdyn]=...
                            SolveMuscleRedundancy_FtildeState_actdyn_GPOPS(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                end
        end
    case 'CasADi'
        switch formulation_actdyn
            case 'DeGroote2016'      
                switch formulation_contdyn
                    case 'lMtildeState'
                        [Time,MExcitation,MActivation,RActivation,TForcetilde,...
                            TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=...
                            SolveMuscleRedundancy_lMtildeState_CasADi(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                    case 'FtildeState'   
                        [Time,MExcitation,MActivation,RActivation,TForcetilde,...
                            TForce,lMtilde,lM,MuscleNames,OptInfo,DatStore]=...
                            SolveMuscleRedundancy_FtildeState_CasADi(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                end

            case 'DeGroote2009'
                switch formulation_contdyn
                    case 'lMtildeState'
                        [Time_actdyn,MExcitation_actdyn,MActivation_actdyn,...
                            RActivation_actdyn,TForcetilde_actdyn,TForce_actdyn,...
                            lMtilde_actdyn,lM_actdyn,MuscleNames_actdyn,...
                            OptInfo_actdyn,DatStore_actdyn]=...
                            SolveMuscleRedundancy_lMtildeState_actdyn_CasADi(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                    case 'FtildeState'   
                        [Time_actdyn,MExcitation_actdyn,MActivation_actdyn,...
                            RActivation_actdyn,TForcetilde_actdyn,TForce_actdyn,...
                            lMtilde_actdyn,lM_actdyn,MuscleNames_actdyn,...
                            OptInfo_actdyn,DatStore_actdyn]=...
                            SolveMuscleRedundancy_FtildeState_actdyn_CasADi(...
                            model_path,IK_path,ID_path,time,Out_path,Misc);
                end
        end
end
