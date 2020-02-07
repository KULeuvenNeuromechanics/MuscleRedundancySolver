function [h] = PlotEMGTracking(Results,DatStore)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

lw = 2; % linewidth

Misc = Results.Misc;
h = figure('Name','Tracking EMG');
nPhases = length(Misc.IKfile);
if nPhases > 1
    hTabGroup = uitabgroup;
end
for i=1:nPhases
    % set the name of the tab
    if nPhases>1
        [path,file,ext]=fileparts(Misc.IKfile{i});
        tab = uitab(hTabGroup, 'Title', file);
        axes('parent',tab);
    end
    % plot the EMG signals for this file
    nEMG = DatStore(i).EMG.nEMG;
    p    = numSubplots(nEMG);
    EMGinds = DatStore(i).EMG.EMGindices;
    % interpolate EMG
    tSim = Results.Time(i).MTE;
    EMG = ppval(DatStore(i).EMG.EMGspline,tSim)';
    for j=1:nEMG
        subplot(p(1),p(2),j);
        % simulated muscle activity - parameter estimation
        Cs = [89, 135, 189]./255;
        eSim = Results.MExcitation(i).MTE(EMGinds(j),:)';
        tSim = Results.Time(i).MTE;
        plot(tSim(1:end-1),eSim,'Color',Cs,'LineWidth',lw); hold on;
        % measured EMG
        EMGscaled = EMG(:,j).*Results.Param.EMGscale(j);
        plot(tSim,EMGscaled,'--k','LineWidth',lw);
        % simulated muscle activity - parameter estimation
        if Misc.ValidationBool
            Cs = [161, 116, 64]./255;
            eSim = Results.MExcitation(i).validationMRS(EMGinds(j),:);
            tSim = Results.Time(i).validationMRS;
            plot(tSim(1:end-1),eSim,'Color',Cs,'LineWidth',lw);
        end
        if j==nEMG
            if Misc.ValidationBool
                legend('Parameter estimation','EMG','Validation');
            else
                legend('Parameter estimation','EMG');
            end
        end
        title(DatStore(i).EMG.EMGselection{j});
    end
end



end

