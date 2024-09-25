function [h] = PlotEstimatedParameters_1side(Results,Misc)
% Plot optimized muscle-tendon parameters relative to their generic value
% and provided bounds.

h = figure('Name','Opt Param');
hTabGroup = uitabgroup;

% plot estimated lmOpt
tab = uitab(hTabGroup, 'Title', 'Opt Fiber length');
axes('parent',tab);
nMus = Misc.nAllMuscList;
p    = numSubplots(nMus);
Cs = [89, 135, 189]./255;
for i=1:nMus
    subplot(p(1),p(2),i)
    or = Results.Param.Original.lMo(i);
    est = Results.Param.Estimated.lMo(i);
    lb = Results.Param.Bound.lMo.lb;    ub = Results.Param.Bound.lMo.ub;
    b = bar(1,or);  b.FaceColor = [0 0 0]; hold on;
    b = bar(2,est); b.FaceColor = Cs;
    l = line([0 3],repmat(or*lb,1,2)); l.Color = [0 0 0];
    l = line([0 3],repmat(or*ub,1,2)); l.Color = [0 0 0];
    set(gca,'XTick',[]);
    title(Misc.allMuscleList{i},'interpreter','none');
end
legend('Original','Estimated','Bound');

% plot estimated tendon stiffness
tab = uitab(hTabGroup, 'Title', 'Tendon Stiffness');
axes('parent',tab);
Cs = [89, 135, 189]./255;
for i=1:nMus
    subplot(p(1),p(2),i)
    or = Results.Param.Original.kT(i);
    est = Results.Param.Estimated.kT(i);
    lb = Results.Param.Bound.kT.lb;    ub = Results.Param.Bound.kT.ub;
    b = bar(1,or);  b.FaceColor = [0 0 0]; hold on;
    b = bar(2,est); b.FaceColor = Cs;
    l = line([0 3],repmat(or*lb,1,2)); l.Color = [0 0 0];
    l = line([0 3],repmat(or*ub,1,2)); l.Color = [0 0 0];
    set(gca,'XTick',[]);
    title(Misc.allMuscleList{i},'interpreter','none');
end
legend('Original','Estimated','Bound');

% plot estimated tendon slack length
tab = uitab(hTabGroup, 'Title', 'Tendon slack length');
axes('parent',tab);
Cs = [89, 135, 189]./255;
for i=1:nMus
    subplot(p(1),p(2),i)
    or = Results.Param.Original.lTs(i);
    est = Results.Param.Estimated.lTs(i);
    lb = Results.Param.Bound.lTs.lb;    ub = Results.Param.Bound.lTs.ub;
    b = bar(1,or);  b.FaceColor = [0 0 0]; hold on;
    b = bar(2,est); b.FaceColor = Cs;
    l = line([0 3],repmat(or*lb,1,2)); l.Color = [0 0 0];
    l = line([0 3],repmat(or*ub,1,2)); l.Color = [0 0 0];
    set(gca,'XTick',[]);
    title(Misc.allMuscleList{i},'interpreter','none');
end
legend('Original','Estimated','Bound');

% Plot the bounds on EMG
if Misc.EMGconstr 
    tab = uitab(hTabGroup, 'Title', 'EMG scale Factor');
    axes('parent',tab);
    nEMG = length(Misc.scaledEMGmusc);
    bar(Results.Param.EMGscale); hold on;  % scale factor
    plot(1:nEMG,repmat(Results.Param.Bound.EMG.lb,1,nEMG),'--k');  % lower bound
    plot(1:nEMG,repmat(Results.Param.Bound.EMG.ub,1,nEMG),'--k');  % upper bound
    set(gca,'XTick',1:nEMG);
    set(gca,'XTickLabel',Misc.scaledEMGmusc);
    set(gca,'XTickLabelRotation',60);
end

end
