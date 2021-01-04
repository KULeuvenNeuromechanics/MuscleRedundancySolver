
%% Example of how to compute the metabolic energy from a solution of the MRS
%---------------------------------------------------------------------------

% Load solution (run Walking_MRSexample.m first)
R = load(fullfile(pwd,'Results','Walking3_Results.mat'));

% get the mass of the subject using the function GetModelMass
modelmass = getModelMass(R.Misc.model_path);

% use the post processing function to compute the metabolic energy consumption
E = GetMetabFromMRS(R.Results,R.Misc,modelmass);

% plot with metabolic energy in the three models
figure();
t = R.Results.Time.genericMRS(1:end-1); % get time indexes on mesh points
plot(t,sum(E.genericMRS.Bargh2004.Edot)); hold on; % sum of all muscle metab powers
plot(t,sum(E.genericMRS.Umb2003.Edot)); hold on; % sum of all muscle metab powers
plot(t,sum(E.genericMRS.Umb2010.Edot)); hold on; % sum of all muscle metab powers
plot(t,sum(E.genericMRS.Uch2016.Edot)); hold on; % sum of all muscle metab powers
legend('Bargh','Umb2003','Umb2010','Uch2016');
xlabel('Time [s]');
ylabel('Metabolic power [W]');



