%% Batch Test examples

% you can test all examples in batch using this script. This is mainly used
% during development to evaluate if all examples (still) work.
RunBatch = {'Example_EMGWalking','EMGconstrained' ;
    'Example_EMGWalking','EMGdriven_SimpleAnkle';
    'Example_trackUS', 'Example_MRS';
    'Example_trackUS', 'Example_US_Tracking';
    'Walking_DeGrooteetal2016', 'Walking_MRSexample'};

StartPath = pwd;

for i=1:length(RunBatch(:,1))
    cd(fullfile(StartPath,RunBatch{i,1}));
    run(RunBatch{i,2});    
end
cd(StartPath);