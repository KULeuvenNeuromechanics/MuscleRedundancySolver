function [h] = PlotEMGTracking(Results,DatStore,Misc)
% Plot results regarding tracking of provided EMG signals.

lw = 2; % linewidth

h = figure('Name','Tracking EMG');
nPhases = length(Misc.trials_sel);
if nPhases > 1
    hTabGroup = uitabgroup;
end
Cs = linspecer(3);
for i=Misc.trials_sel
    % set the name of the tab
    if nPhases>1
        [path,file,ext]=fileparts(Misc.IKfile{i});
        tab = uitab(hTabGroup, 'Title', file);
        axes('parent',tab);
    end
    % plot the EMG signals for this file
    nEMG = DatStore(i).EMG.nEMG;
    p    = numSubplots(nEMG);
    EMGinds = DatStore(i).EMG.idx_EMGsel(:,3);
    % interpolate EMG
    tSim = Results.Time(i).MTE;
    EMG = ppval(DatStore(i).EMG.EMGspline,tSim)';
    legend_title = {};
    for j=1:nEMG
        subplot(p(1),p(2),j);
        % simulated muscle activity - parameter estimation
        eSim = Results.MExcitation(i).MTE(EMGinds(j),:)';
        tSim = Results.Time(i).MTE;
        plot(tSim(1:end-1),eSim,'Color',Cs(1,:),'LineWidth',lw); hold on;
        if j == 1
            legend_title = [legend_title,'Parameter Estimation'];
        end
        % measured EMG
        EMGscaled = EMG(:,j).*Results.Param.EMGscale(DatStore(i).EMG.idx_EMGsel(j,4));
        plot(tSim,EMGscaled,'--k','LineWidth',lw);
        if j == 1
            legend_title = [legend_title,'Measured'];
        end
        % simulated muscle activity - parameter estimation
        if Misc.ValidationBool
            eSim = Results.MExcitation(i).validationMRS(EMGinds(j),:);
            tSim = Results.Time(i).validationMRS;
            plot(tSim(1:end-1),eSim,'Color',Cs(2,:),'LineWidth',lw);
            if j == 1
                legend_title = [legend_title,'Validation MRS'];
            end
        end
        
        if Misc.MRSBool
            eSim = Results.MExcitation(i).genericMRS(EMGinds(j),:);
            tSim = Results.Time(i).genericMRS;
            plot(tSim(1:end-1),eSim,'Color',Cs(3,:),'LineWidth',lw);
            if j == 1
                legend_title = [legend_title,'Generic MRS'];
            end
        end
        
        if j==nEMG
            legend(legend_title);
        end
        title(DatStore(i).EMG.EMGselection{j});
    end
end

end

