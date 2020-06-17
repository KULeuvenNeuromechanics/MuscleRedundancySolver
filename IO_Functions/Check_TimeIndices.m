function [time] = Check_TimeIndices(Misc,time)
% Adapt user-provided time window to match closest time points in the IK
% solution
for i = 1: Misc.nTrials
    IK = importdata(Misc.IKfile{i});
    tIK = IK.data(:,1);
    t0 = tIK(tIK>= time(i,1)); t0 = t0(1);
    tend = tIK(tIK<= time(i,2)); tend = tend(end);
    time_new = [t0 tend];
    if sum(time(i,:)-time_new) ~= 0
        disp(['Adapted time window to framerate IK solution: ' num2str(time_new)]);
        time(i,:) = time_new;
    end
end

end

