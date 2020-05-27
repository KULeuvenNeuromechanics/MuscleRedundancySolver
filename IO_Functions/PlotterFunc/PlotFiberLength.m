function [h] = PlotFiberLength(Results,DatStore)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


lw = 2;

Misc = Results.Misc;
h = figure('Name','Muscle fiber length');
nPhases = length(Misc.IKfile);
if nPhases > 1
    hTabGroup = uitabgroup;
end

for trial = 1:nPhases
    
    % set the name of the tab
    if nPhases>1
        [path,file,ext]=fileparts(Misc.IKfile{trial});
        tab = uitab(hTabGroup, 'Title', file);
        axes('parent',tab);
    end
    
    % get the number muscles with fiber length tracking
    Minds = DatStore(trial).free_lMo(:);
    nMus  = length(Minds);
    p     = numSubplots(nMus);
    
    for j=1:nMus
       subplot(p(1),p(2),j)
       % simulated fiber length -- parameter estimation
       Cs = [89, 135, 189]./255;
       lMSim = Results.lM(trial).MTE(Minds(j),:)';
       tSim = Results.Time(trial).MTE;
       plot(tSim,lMSim,'Color',Cs,'LineWidth',lw); hold on;
       
       % simulated fiber length -- validation simulation
       if Misc.ValidationBool
            Cs = [161, 116, 64]./255;
            lMSim = Results.lM(trial).validationMRS(Minds(j),:);
            tSim = Results.Time(trial).validationMRS;
            plot(tSim,lMSim,'Color',Cs,'LineWidth',lw);
       end 
       % measured fiber length
       if Misc.UStracking
           plot(Time(trial).MTE,USTracking(trial).data(:,Minds(j))/1000,'--k','LineWidth',lw);
       end
       title(DatStore(trial).MuscleNames{Minds(j)});
    end
    if Misc.ValidationBool==1 && Misc.UStracking==1        
        legend('Parameter estimation','Validation','EMG');
    elseif Misc.ValidationBool==0 && Misc.UStracking==1  
        legend('Parameter estimation','EMG');
    elseif Misc.ValidationBool==0 && Misc.UStracking==0
        legend('Parameter estimation');
    end    
end



end

