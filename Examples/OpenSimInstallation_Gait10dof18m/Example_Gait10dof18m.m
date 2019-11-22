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
filepath=which('Example_Gait10dof18m.m'); 
[DirExample,~,~]=fileparts(filepath); 
[DirExample2,~,~]=fileparts(DirExample); 
[MainDir,~]=fileparts(DirExample2);
addpath(genpath(MainDir));

% Needed Input Arguments
Datapath='C:\OpenSim 3.3\Models\Gait10dof18musc\OutputReference';
IK_path=fullfile(Datapath,'IK','subject01_walk_IK.mot');
ID_path=[]; % Compute ID from the external loads
model_path=fullfile(Datapath,'subject01.osim');
time=[0.7 1.4]; % Part of the right stance phase
Out_path=fullfile(MainDir,'Examples','OpenSimInstallation_Gait10dof18m','Results');

Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r'};
Misc.Loads_path=fullfile(Datapath,'ExperimentalData','subject01_walk_grf.xml');
Misc.ID_ResultsPath=fullfile(Datapath,'ID','inversedynamics.sto');

% Optional Input Arguments
% Here is an example of how to adjust the Achilles tendon stiffness.
% We first add the input argument MuscleNames_Input with ALL muscles 
% that actuate the degrees of freedom listed in DofNames_Input.
Misc.MuscleNames_Input={'hamstrings_r','bifemsh_r','glut_max_r',...
    'iliopsoas_r','rect_fem_r','vasti_r','gastroc_r','soleus_r',...
    'tib_ant_r'}; 
% We then change the compliance of the Achilles tendon by changing the 
% parameter Atendon of the gastrocnemius and the soleus. The default 
% value is 35 and a lower value will result in a more compliant tendon.
Misc.Atendon=[35,35,35,35,35,35,15,15,35];

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
