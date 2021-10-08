function ParamsToOsim(muscleParams,muscleNames,modelPath,outPath,newModelFile)
%---------------------------------------------------------------------------
%ParamsToOsim
%   This function writes muscle parameters to a new osim model, based on a
%   given model.
%
% INPUT:
%   -muscleParams-
%   * 4 x nMuscles structure with
%       * FMo         : maximal isometric force
%       * lMo         : optimal fiber length
%       * lTs         : tendon slack length 
%       * alphao      : penation angle at lMo
%  
%   -muscleNames-  
%   * cell structure with the muscles from which muscle parameters are 
%   available
%
%   -modelPath-
%   * path to the original model file
%
%   -outPath-
%   * output file path
% 
%   -newModelFile- 
%   * name of the new osim model
% 
% OUTPUT:
%   * new osim model based on loaded model, saved in the specified folder
%
% Original author: Bram Van Den Bosch
% Original date: 12/05/2021
%
% Last edit by: Bram Van Den Bosch
% Last edit date: 18/05/2021
%---------------------------------------------------------------------------

import org.opensim.modeling.* % import OpenSim libraries

FMo    = muscleParams.FMo;
lMo    = muscleParams.lMo;
lTs    = muscleParams.lTs;
alphao = muscleParams.alphao;

listing       = dir(modelPath);
modelFilePath = listing.folder;
modelFile     = listing.name;

% load the original model and initialize
model = Model(fullfile(modelFilePath, modelFile));
model.initSystem();

muscles  = model.getMuscles(); % get the muscles of the original model
nMuscles = length(muscleNames); % count the muscles with new parameters

% loop through muscles with new parameters
for i = 1:nMuscles 
    % get the muscle in the generic model
    currentMuscle = muscles.get(muscleNames{i});

    % set new parameters
    currentMuscle.setMaxIsometricForce(FMo(i));
    currentMuscle.setOptimalFiberLength(lMo(i));
    currentMuscle.setTendonSlackLength(lTs(i));
    currentMuscle.setPennationAngleAtOptimalFiberLength(alphao(i));    
end
 
% save the new model to an osim file
model.print(fullfile(outPath, newModelFile));

end


