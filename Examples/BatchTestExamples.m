
%% Batch run all examples 

ExPath = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\MuscleTendonEstimator\Examples';
cd(fullfile(ExPath,'Walking_DeGrooteetal2016'));
run(fullfile(ExPath,'Walking_DeGrooteetal2016','Walking_MRSexample')); 

ExPath = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\MuscleTendonEstimator\Examples';
cd(fullfile(ExPath,'JPCB_ORLAU_Example'));
run(fullfile(ExPath,'JPCB_ORLAU_Example','ng10203_optimised'));

ExPath = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\MuscleTendonEstimator\Examples';
cd(fullfile(ExPath,'Example_trackUS'));
run(fullfile(ExPath,'Example_trackUS','Example_MRS'));
run(fullfile(ExPath,'Example_trackUS','Example_MRS_minimal'));
run(fullfile(ExPath,'Example_trackUS','Example_US_Tracking'));

ExPath = 'C:\Users\u0088756\Documents\FWO\Software\GitProjects\MuscleTendonEstimator\Examples';
cd(fullfile(ExPath,'Example_TendonEstimator_RunningData'));
run(fullfile(ExPath,'Example_TendonEstimator_RunningData','Example_UStracking'));
cd(ExPath);
