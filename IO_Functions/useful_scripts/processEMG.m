function output = processEMG(inputs, AnalogRate)
% Purpose:  From raw EMG data as input, the method performs the
%           filtering, rectifying, and normalizing operations needed to
%           obtain an EMG envelope.
% ADDED 23/01/2020 JPCB: 
% The script will now process multiple EMG signles 
% 
% USAGE:  [proEMG, rectEMG] = processEMG(rawEMG, AnalogRate)
%
% Input:   inputs is a structure with multiple raw EMG signals with columns of EMG channels
%          AnalogRate is the sampling rate of the raw EMG data
%
% Output:  output is a strucutre with:
%           proEMG is column matrix of normalized EMG envelopes
%           rectEMG (optional) of just norm rectified EMG
%
% NOTE:    rawEMG signals must contain at least 12 data points
%          in order to apply a fourth order lowpass filter
%
% ASeth Nov-07, streamlined version of ASA, 9-05
% Stanford University
total_nAnalogFrames = 0;

number_of_trials = length(inputs);
% number_of_channels = size_of_signals(2);
for i = 1:number_of_trials
    rawEMG = inputs(i).selected_column_signals;

    [processing(i).nAnalogFrames, processing(i).nc] = size(rawEMG);
    total_nAnalogFrames = total_nAnalogFrames + processing(i).nAnalogFrames;
    
    % Apply 4th order, 0-lag, Butterworth band-pass filter to raw EMG signal.
    order = 4;
    fs = AnalogRate;
    if fs >= 1000
        cutoff = [20 400];                % default
        % cutoff = [80 400];              % use when there is 60 Hz noise
    elseif fs == 600                      % For Delaware EMG data
        cutoff = [11 222];                % default
        % cutoff = [44 222];              % use when there is 60 Hz noise
    end
    [b, a] = butter(order/2, cutoff/(0.5*fs));
    bandEMG = filter(b, a, rawEMG, [], 1);
    
    % Rectify the filtered EMG signal.
    processing(i).rectEMG = abs(bandEMG);
    
    % Apply 4th order, 0-lag, Butterworth low-pass filter to rectified signal.
    order = 4;
    cutoff = 10;
    [b, a] = butter(order, cutoff/(0.5*fs));
    processing(i).lowEMG = filter(b, a, processing(i).rectEMG, [], 1);
    
    processing(i).onesMat = ones(processing(i).nAnalogFrames,1);
    
    if i ~= 1 
       all_rectEMG = cat(1, all_rectEMG, processing(i).rectEMG);
       all_lowEMG = cat(1, all_lowEMG, processing(i).lowEMG);
    else 
        all_rectEMG = processing(i).rectEMG;
        all_lowEMG = processing(i).lowEMG;
    end
end



% Normalize rectified and low-pass filtered EMG signals.
nMinSamples = round(0.01 * total_nAnalogFrames);  
                                % average 1% of samples to get "min".
sortRect = sort(all_rectEMG);           
minRect = mean(sortRect(1:nMinSamples,:));
maxRect = max(all_rectEMG);
sortLow = sort(all_lowEMG);
minLow = mean(sortLow(1:nMinSamples,:));
maxLow = max(all_lowEMG);


for i = 1:number_of_trials
    processing(i).normRect = (processing(i).rectEMG - processing(i).onesMat*minRect)./(processing(i).onesMat*(maxRect - minRect));
    processing(i).normLow = (processing(i).lowEMG - processing(i).onesMat*minLow)./ (processing(i).onesMat*(maxLow - minLow));
    
    % Return processed EMG data.        
    output(i).rectEMG = processing(i).normRect;
    output(i).proEMG = processing(i).normLow;   
end         
             

