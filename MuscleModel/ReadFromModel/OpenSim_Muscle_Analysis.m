function [] = OpenSim_Muscle_Analysis(motion_file,model_sel,output_path,event)
%OpenSim_Muscle_Analysis Executes a muscle analsysis from the command line
%   Detailed explanation goes here

import org.opensim.modeling.*

% get the analysis tool and change variables
[FunctionPath,~]=fileparts(mfilename('fullpath'));
path_generic_file=fullfile(FunctionPath,'settings_Muscle_analysis.xml');
tool=AnalyzeTool(path_generic_file,false);
tool.setLoadModelAndInput(true)
osimModel=Model(model_sel);
tool.setModel(osimModel);
tool.setResultsDir(output_path);
tool.setInitialTime(event(1));
tool.setFinalTime(event(2));
[~, name, ~]=fileparts(motion_file);
tool.setName(name);

% run the analysis
tool.setModelFilename(model_sel);
tool.setCoordinatesFileName(motion_file);
out_path_xml=fullfile(output_path,['muscle_analysis_' name '.xml']);
tool.print(out_path_xml);
tool.run;


end

