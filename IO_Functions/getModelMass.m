function [mass,segments,segment_mass] = getModelMass(model_path)
%getModelMass Reads the mass of the model from an .osim file
%   Detailed explanation goes here

% import opensim libraries
import org.opensim.modeling.*

% import the model
m = Model(model_path);

% get number of segments
bodies = m.getBodySet();
nbody = bodies.getSize();

% loop over all bodies
segment_mass = nan(nbody,1);
segments = cell(nbody,1);
for m = 1:nbody
    % get the mass of the segment
    segment_mass(m) = bodies.get(m-1).getMass();
    % get the name of the segment
    segments{m} = char(bodies.get(m-1).getName());    
end

% get total body mass
mass = sum(segment_mass);


end

