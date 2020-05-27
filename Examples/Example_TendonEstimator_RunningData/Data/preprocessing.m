% clear variables and commond window
clear all; clc;

%% Input information
% Install instructions:
%   Add the main path ...\solvemuscleredundancy_dev to your matlab folder,
%   using the following command (in my case)
ExamplePath = pwd;


tagList = {'PSF', '+8', '+15', '8', '15'}

for i = 1:length(tagList)
    
    tag = char(tagList(i));
    DataPath = [pwd '\KS'];
    IKdata = importdata([DataPath '\KS_' tag '_GM_03.mot']);
    jointPositionNames = IKdata.colheaders;
    
    DataPath = [pwd '\ID'];
    IDdata = importdata([DataPath '\subject12_' tag '_GM_03.sto']);
    jointMomentNames = IDdata.colheaders;
    
    
    
    %% We generate motion files for each stride for the IK, ID and US
    
    
    % PSF condition
    load([DataPath '\' tag '_GM_03.mat']);
    
    load([tag '_SOL_03.mat']);
    resultsSOL = results; clear results;
    load([tag '_GM_03.mat']);
    resultsGM = results; clear results;
    
    % Cycle 1
    time_cycle = resultsSOL.SOL.time.stride.stride2;
    t0 = round(time_cycle(1)+0.00499999,2); tend = round(time_cycle(end)-0.00499999,2);
    time_cycleInterp = t0:0.005:tend; time_cycleInterp = time_cycleInterp';
    
    ankleAngle_cycle = interp1(time_cycle,resultsSOL.IK.ankle.angle.stride.stride2,time_cycleInterp,'pchip');
    kneeAngle_cycle = interp1(time_cycle,resultsSOL.IK.knee.angle.stride.stride2,time_cycleInterp,'pchip');
    hipAngle_cycle = interp1(time_cycle,resultsSOL.IK.hip.angle.stride.stride2,time_cycleInterp,'pchip');
    
    
    IDtorque_cycle = interp1(time_cycle,ID.moments.stride.stride2,time_cycleInterp,'pchip');
    
    % Inverse kinematics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointPositionNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,19) = ankleAngle_cycle; %left ankle
    dataMatrix(:,18) = kneeAngle_cycle; %left knee
    dataMatrix(:,15) = hipAngle_cycle; %left hip
    generateMotFile(dataMatrix, jointPositionNames, [tag '_IK_stride2.mot']);
    
    % Inverse dynamics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointMomentNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2:end) = IDtorque_cycle;
    generateMotFile(dataMatrix, jointMomentNames, [tag '_ID_stride2.mot']);
    
    % Ultrasound
    dataMatrix = zeros(size(time_cycleInterp,1),3);
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2) = 1*interp1(time_cycle,resultsSOL.SOL.length.stride.stride2,time_cycleInterp,'pchip');
    dataMatrix(:,3) = 1*interp1(time_cycle,resultsGM.GM.length.stride.stride2,time_cycleInterp,'pchip');
    USnames = {'time','soleus_l','med_gas_l'};
    generateMotFile(dataMatrix, USnames, [tag '_US_stride2.mot']);
    
    
    
    
    % Cycle 2
    time_cycle = resultsSOL.SOL.time.stride.stride3;
    t0 = round(time_cycle(1)+0.00499999,2); tend = round(time_cycle(end)-0.00499999,2);
    time_cycleInterp = t0:0.005:tend; time_cycleInterp = time_cycleInterp';
    
    ankleAngle_cycle = interp1(time_cycle,resultsSOL.IK.ankle.angle.stride.stride3,time_cycleInterp,'pchip');
    kneeAngle_cycle = interp1(time_cycle,resultsSOL.IK.knee.angle.stride.stride3,time_cycleInterp,'pchip');
    hipAngle_cycle = interp1(time_cycle,resultsSOL.IK.hip.angle.stride.stride3,time_cycleInterp,'pchip');
    
    
    IDtorque_cycle = interp1(time_cycle,ID.moments.stride.stride3,time_cycleInterp,'pchip');
    
    % Inverse kinematics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointPositionNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,19) = ankleAngle_cycle; %left ankle
    dataMatrix(:,18) = kneeAngle_cycle; %left knee
    dataMatrix(:,15) = hipAngle_cycle; %left hip
    generateMotFile(dataMatrix, jointPositionNames, [tag '_IK_stride3.mot']);
    
    % Inverse dynamics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointMomentNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2:end) = IDtorque_cycle;
    generateMotFile(dataMatrix, jointMomentNames, [tag '_ID_stride3.mot']);
    
    % Ultrasound
    dataMatrix = zeros(size(time_cycleInterp,1),3);
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2) = 1*interp1(time_cycle,resultsSOL.SOL.length.stride.stride3,time_cycleInterp,'pchip');
    dataMatrix(:,3) = 1*interp1(time_cycle,resultsGM.GM.length.stride.stride3,time_cycleInterp,'pchip');
    USnames = {'time','soleus_l','med_gas_l'};
    generateMotFile(dataMatrix, USnames, [tag '_US_stride3.mot']);
    
    
    
    % Cycle 3
    time_cycle = resultsSOL.SOL.time.stride.stride4;
    t0 = round(time_cycle(1)+0.00499999,2); tend = round(time_cycle(end)-0.00499999,2);
    time_cycleInterp = t0:0.005:tend; time_cycleInterp = time_cycleInterp';
    
    ankleAngle_cycle = interp1(time_cycle,resultsSOL.IK.ankle.angle.stride.stride4,time_cycleInterp,'pchip');
    kneeAngle_cycle = interp1(time_cycle,resultsSOL.IK.knee.angle.stride.stride4,time_cycleInterp,'pchip');
    hipAngle_cycle = interp1(time_cycle,resultsSOL.IK.hip.angle.stride.stride4,time_cycleInterp,'pchip');
    
    
    IDtorque_cycle = interp1(time_cycle,ID.moments.stride.stride4,time_cycleInterp,'pchip');
    
    % Inverse kinematics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointPositionNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,19) = ankleAngle_cycle; %left ankle
    dataMatrix(:,18) = kneeAngle_cycle; %left knee
    dataMatrix(:,15) = hipAngle_cycle; %left hip
    generateMotFile(dataMatrix, jointPositionNames, [tag '_IK_stride4.mot']);
    
    % Inverse dynamics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointMomentNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2:end) = IDtorque_cycle;
    generateMotFile(dataMatrix, jointMomentNames, [tag '_ID_stride4.mot']);
    
    % Ultrasound
    dataMatrix = zeros(size(time_cycleInterp,1),3);
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2) = 1*interp1(time_cycle,resultsSOL.SOL.length.stride.stride4,time_cycleInterp,'pchip');
    dataMatrix(:,3) = 1*interp1(time_cycle,resultsGM.GM.length.stride.stride4,time_cycleInterp,'pchip');
    USnames = {'time','soleus_l','med_gas_l'};
    generateMotFile(dataMatrix, USnames, [tag '_US_stride4.mot']);
    
    
    
    % Cycle 4
    time_cycle = resultsSOL.SOL.time.stride.stride5;
    t0 = round(time_cycle(1)+0.00499999,2); tend = round(time_cycle(end)-0.00499999,2);
    time_cycleInterp = t0:0.005:tend; time_cycleInterp = time_cycleInterp';
    
    ankleAngle_cycle = interp1(time_cycle,resultsSOL.IK.ankle.angle.stride.stride5,time_cycleInterp,'pchip');
    kneeAngle_cycle = interp1(time_cycle,resultsSOL.IK.knee.angle.stride.stride5,time_cycleInterp,'pchip');
    hipAngle_cycle = interp1(time_cycle,resultsSOL.IK.hip.angle.stride.stride5,time_cycleInterp,'pchip');
    
    
    IDtorque_cycle = interp1(time_cycle,ID.moments.stride.stride5,time_cycleInterp,'pchip');
    
    % Inverse kinematics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointPositionNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,19) = ankleAngle_cycle; %left ankle
    dataMatrix(:,18) = kneeAngle_cycle; %left knee
    dataMatrix(:,15) = hipAngle_cycle; %left hip
    generateMotFile(dataMatrix, jointPositionNames, [tag '_IK_stride5.mot']);
    
    % Inverse dynamics
    dataMatrix = zeros(size(time_cycleInterp,1),size(jointMomentNames,2));
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2:end) = IDtorque_cycle;
    generateMotFile(dataMatrix, jointMomentNames, [tag '_ID_stride5.mot']);
    
    % Ultrasound
    dataMatrix = zeros(size(time_cycleInterp,1),3);
    dataMatrix(:,1) = time_cycleInterp;
    dataMatrix(:,2) = 1*interp1(time_cycle,resultsSOL.SOL.length.stride.stride5,time_cycleInterp,'pchip');
    dataMatrix(:,3) = 1*interp1(time_cycle,resultsGM.GM.length.stride.stride5,time_cycleInterp,'pchip');
    USnames = {'time','soleus_l','med_gas_l'};
    generateMotFile(dataMatrix, USnames, [tag '_US_stride5.mot']);
    
    
    if i < 5
        % Cycle 5
        time_cycle = resultsSOL.SOL.time.stride.stride6;
        t0 = round(time_cycle(1)+0.00499999,2); tend = round(time_cycle(end)-0.00499999,2);
        time_cycleInterp = t0:0.005:tend; time_cycleInterp = time_cycleInterp';
        
        ankleAngle_cycle = interp1(time_cycle,resultsSOL.IK.ankle.angle.stride.stride6,time_cycleInterp,'pchip');
        kneeAngle_cycle = interp1(time_cycle,resultsSOL.IK.knee.angle.stride.stride6,time_cycleInterp,'pchip');
        hipAngle_cycle = interp1(time_cycle,resultsSOL.IK.hip.angle.stride.stride6,time_cycleInterp,'pchip');
        
        
        IDtorque_cycle = interp1(time_cycle,ID.moments.stride.stride6,time_cycleInterp,'pchip');
        
        % Inverse kinematics
        dataMatrix = zeros(size(time_cycleInterp,1),size(jointPositionNames,2));
        dataMatrix(:,1) = time_cycleInterp;
        dataMatrix(:,19) = ankleAngle_cycle; %left ankle
        dataMatrix(:,18) = kneeAngle_cycle; %left knee
        dataMatrix(:,15) = hipAngle_cycle; %left hip
        generateMotFile(dataMatrix, jointPositionNames, [tag '_IK_stride6.mot']);
        
        % Inverse dynamics
        dataMatrix = zeros(size(time_cycleInterp,1),size(jointMomentNames,2));
        dataMatrix(:,1) = time_cycleInterp;
        dataMatrix(:,2:end) = IDtorque_cycle;
        generateMotFile(dataMatrix, jointMomentNames, [tag '_ID_stride6.mot']);
        
        % Ultrasound
        dataMatrix = zeros(size(time_cycleInterp,1),3);
        dataMatrix(:,1) = time_cycleInterp;
        dataMatrix(:,2) = 1*interp1(time_cycle,resultsSOL.SOL.length.stride.stride6,time_cycleInterp,'pchip');
        dataMatrix(:,3) = 1*interp1(time_cycle,resultsGM.GM.length.stride.stride6,time_cycleInterp,'pchip');
        USnames = {'time','soleus_l','med_gas_l'};
        generateMotFile(dataMatrix, USnames, [tag '_US_stride6.mot']);
    end
    
end
