function [h] = PlotEstimatedParameters(Results,DatStore,Misc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


h = figure('Name','Opt Param');
hTabGroup = uitabgroup;

% plot estimated lmOpt
tab = uitab(hTabGroup, 'Title', 'Opt Fiber length');
axes('parent',tab);
iM = DatStore(1).free_lMo(:); % Muscle indexes
nMus = length(iM);
p    = numSubplots(nMus);
Cs = [89, 135, 189]./255;
for i=1:nMus
    subplot(p(1),p(2),i)
    or = Results.Param.Original.lOpt(iM(i));
    est = Results.Param.Estimated.lOpt(iM(i));
    lb = Results.Param.Bound.lOp.lb;    ub = Results.Param.Bound.lOp.ub;
    b = bar(1,or);  b.FaceColor = [0 0 0]; hold on;
    b = bar(2,est); b.FaceColor = Cs;
    l = line([0 3],repmat(or*lb,1,2)); l.Color = [0 0 0];
    l = line([0 3],repmat(or*ub,1,2)); l.Color = [0 0 0];
    set(gca,'XTick',[]);
    title(DatStore(1).MuscleNames{iM(i)},'interpreter','none');
end
legend('Original','Estimated','Bound');

% plot estimated tendon stiffness
tab = uitab(hTabGroup, 'Title', 'Tendon Stiffness');
axes('parent',tab);
iM = DatStore(1).free_kT(:); % Muscle indexes
nMus = length(iM);
p    = numSubplots(nMus);
Cs = [89, 135, 189]./255;
for i=1:nMus
    subplot(p(1),p(2),i)
    or = Results.Param.Original.ATendon(iM(i));
    est = Results.Param.Estimated.ATendon(iM(i));
    lb = Results.Param.Bound.kT.lb;    ub = Results.Param.Bound.kT.ub;
    b = bar(1,or);  b.FaceColor = [0 0 0]; hold on;
    b = bar(2,est); b.FaceColor = Cs;
    l = line([0 3],repmat(or*lb,1,2)); l.Color = [0 0 0];
    l = line([0 3],repmat(or*ub,1,2)); l.Color = [0 0 0];
    set(gca,'XTick',[]);
    title(DatStore(1).MuscleNames{iM(i)},'interpreter','none');
end
legend('Original','Estimated','Bound');

% plot estimated tendon slack length
tab = uitab(hTabGroup, 'Title', 'Tendon slack length');
axes('parent',tab);
iM = DatStore(1).free_lMo(:); % Muscle indexes
nMus = length(iM);
p    = numSubplots(nMus);
Cs = [89, 135, 189]./255;
for i=1:nMus
    subplot(p(1),p(2),i)
    or = Results.Param.Original.L_Slack(iM(i));
    est = Results.Param.Estimated.L_Slack(iM(i));
    lb = Results.Param.Bound.lTs.lb;    ub = Results.Param.Bound.lTs.ub;
    b = bar(1,or);  b.FaceColor = [0 0 0]; hold on;
    b = bar(2,est); b.FaceColor = Cs;
    l = line([0 3],repmat(or*lb,1,2)); l.Color = [0 0 0];
    l = line([0 3],repmat(or*ub,1,2)); l.Color = [0 0 0];
    set(gca,'XTick',[]);
    title(DatStore(1).MuscleNames{iM(i)},'interpreter','none');
end
legend('Original','Estimated','Bound');

% Plot the bounds on EMG
if Misc.EMGconstr 
    tab = uitab(hTabGroup, 'Title', 'EMG scale Factor');
    axes('parent',tab);
    nEMG = DatStore(1).EMG.nEMG;
    bar(Results.Param.EMGscale); hold on;  % scale factor
    plot(1:nEMG,repmat(Results.Param.Bound.EMG.lb,1,nEMG),'--k');  % lower bound
    plot(1:nEMG,repmat(Results.Param.Bound.EMG.ub,1,nEMG),'--k');  % upper bound
    set(gca,'XTick',1:DatStore(1).EMG.nEMG);
    set(gca,'XTickLabel',DatStore(1).EMG.EMGselection);
    set(gca,'XTickLabelRotation',60);
end







end

