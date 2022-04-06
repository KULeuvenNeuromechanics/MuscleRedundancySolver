function [ h ] = PlotStates(Results,DatStore,Misc)

% PlotStates: plots the states of all the hill-type muscles in the model.
% Allows for comparison between the solution of the different solved
% redundancy problems.

BoolParamOpt = 0;
if Misc.UStracking == 1 || Misc.EMGconstr == 1
    BoolParamOpt = 1;
end

h = figure('Name','Optimal states and controls');
hTabGroup = uitabgroup;
% colors for different simulations
Cs = linspecer(3);
for trial = 1:Misc.nTrials
    % plot optimal activations
    tab = uitab(hTabGroup, 'Title', [Misc.trialName{trial} ': activation']);
    axes('parent',tab);
    nMus = length(DatStore(trial).MuscleNames);
    lw   = 4;
    p    = numSubplots(nMus);
    legend_title = {};
    for i=1:nMus
        subplot(p(1),p(2),i)
        if Misc.UStracking==1 || Misc.EMGconstr==1
            plot(Results.Time(trial).MTE,Results.MActivation(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',lw); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1 && BoolParamOpt
            plot(Results.Time(trial).validationMRS,Results.MActivation(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',lw./2);  hold on;
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS,Results.MActivation(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',lw./2);  hold on;
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end        
        title(DatStore(trial).MuscleNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title)
    
    % plot optimal lMtilde
    tab = uitab(hTabGroup, 'Title', [Misc.trialName{trial} ': normalized FL']);
    axes('parent',tab);
    legend_title = {};
    for i=1:nMus
        subplot(p(1),p(2),i)
        if Misc.UStracking==1 || Misc.EMGconstr
            plot(Results.Time(trial).MTE,Results.lMtildeopt(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',lw); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1 && BoolParamOpt
            plot(Results.Time(trial).validationMRS,Results.lMtildeopt(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',lw/2); hold on;
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS,Results.lMtildeopt(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',lw/2);  hold on;
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end
        
        title(DatStore(trial).MuscleNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title)
    
    % plot optimal lM
    tab = uitab(hTabGroup, 'Title', [Misc.trialName{trial} ': FL']);
    axes('parent',tab);
    lw = 2;
    legend_title = {};
    for i=1:nMus
        subplot(p(1),p(2),i)
        if Misc.UStracking==1 || Misc.EMGconstr
            plot(Results.Time(trial).MTE,Results.lM(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',lw); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1 && BoolParamOpt
            plot(Results.Time(trial).validationMRS,Results.lM(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',lw/2);  hold on;
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS,Results.lM(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',lw/2);  hold on;
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end
        
        title(DatStore(trial).MuscleNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title)
    
    
    % plot residual actuators
    tab = uitab(hTabGroup, 'Title', [Misc.trialName{trial} ': Residual actuators']);
    axes('parent',tab);
    nDOF = DatStore(trial).nDOF;
    p    = numSubplots(nDOF);    
    legend_title = {};
    for i=1:nDOF
        subplot(p(1),p(2),i)
        if Misc.UStracking==1 || Misc.EMGconstr
            plot(Results.Time(trial).MTE(1:end-1),Results.RActivation(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',lw); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1 && BoolParamOpt
            plot(Results.Time(trial).validationMRS(1:end-1),Results.RActivation(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',lw/2);  hold on;
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS(1:end-1),Results.RActivation(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',lw/2);  hold on;
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end        
        title(DatStore(trial).DOFNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title);    
end
end

