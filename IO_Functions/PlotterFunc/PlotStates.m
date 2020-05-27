function [ h ] = PlotStates( Results,DatStore,Misc )
h = figure('Name','Optimal states and controls');
hTabGroup = uitabgroup;
    Cs = linspecer(3);


for trial = 1:length(Results.MActivation(:))
    % plot optimal activations
    tab = uitab(hTabGroup, 'Title', ['Trial ' num2str(trial) ': activation']);
    axes('parent',tab);
    lw = 2;
    nMus = DatStore.NMuscles;
    p    = numSubplots(nMus);
    legend_title = {};
    for i=1:nMus
        subplot(p(1),p(2),i)
        if Misc.UStracking==1
            plot(Results.Time(trial).MTE,Results.MActivation(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',4); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1
            plot(Results.Time(trial).validationMRS,Results.MActivation(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',2);
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS,Results.MActivation(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',1);
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end
        
        title(DatStore(1).MuscleNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title)
    
    % plot optimal lMtilde
    tab = uitab(hTabGroup, 'Title', ['Trial ' num2str(trial) ': normalized FL']);
    axes('parent',tab);
    lw = 2;
    nMus = DatStore.NMuscles;
    p    = numSubplots(nMus);
    legend_title = {};
    for i=1:nMus
        subplot(p(1),p(2),i)
        if Misc.UStracking==1
            plot(Results.Time(trial).MTE,Results.lMtildeopt(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',4); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1
            plot(Results.Time(trial).validationMRS,Results.lMtildeopt(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',2);
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS,Results.lMtildeopt(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',1);
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end
        
        title(DatStore(1).MuscleNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title)
    
    % plot optimal lM
    tab = uitab(hTabGroup, 'Title', ['Trial ' num2str(trial) ': FL']);
    axes('parent',tab);
    lw = 2;
    nMus = DatStore.NMuscles;
    p    = numSubplots(nMus);
    legend_title = {};
    for i=1:nMus
        subplot(p(1),p(2),i)
        if Misc.UStracking==1
            plot(Results.Time(trial).MTE,Results.lM(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',4); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1
            plot(Results.Time(trial).validationMRS,Results.lM(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',2);
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS,Results.lM(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',1);
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end
        
        title(DatStore(1).MuscleNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title)
    
    
    % plot residual actuators
    tab = uitab(hTabGroup, 'Title', ['Trial ' num2str(trial) ': Residual actuators']);
    axes('parent',tab);
    lw = 2;
    nDOF = DatStore.nDOF;
    p    = numSubplots(nDOF);
    
    legend_title = {};
    for i=1:nDOF
        subplot(p(1),p(2),i)
        if Misc.UStracking==1
            plot(Results.Time(trial).MTE(1:end-1),Results.RActivation(trial).MTE(i,:),'Color',Cs(1,:),'LineWidth',4); hold on;
            if i == 1
                legend_title = [legend_title,'Parameter Estimation'];
            end
        end
        if Misc.ValidationBool==1
            plot(Results.Time(trial).validationMRS(1:end-1),Results.RActivation(trial).validationMRS(i,:),'Color',Cs(2,:),'LineWidth',2);
            if i == 1
                legend_title = [legend_title;'Validation MRS'];
            end
        end
        if Misc.MRSBool==1
            plot(Results.Time(trial).genericMRS(1:end-1),Results.RActivation(trial).genericMRS(i,:),'Color',Cs(3,:),'LineWidth',1);
            if i == 1
                legend_title = [legend_title;'Generic MRS'];
            end
        end
        
        title(DatStore(1).MuscleNames{i},'interpreter','none');
        xlim([DatStore(trial).time(1) DatStore(trial).time(end)])
    end
    legend(legend_title)
    
end
end

