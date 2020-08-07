function [Energy] = GetMetabFromMRS(Results,Misc)
% GetMetabFromMR Computes metabolic energy consumption from the muscle
% redundancy solver results.
% input arguments:
%   Results: Results structure (output argument of the MRS)
%   Misc:M isc structure (output argument of the MRS)
% output arguments:
%   Energy 	... Edot: metabolic power without basal rate for each muscle [W]
%           ... Edot_model: sum of metabolic power of all selected muscles with basal rate [W]
%           ... W: total metabolic work (without basal rate) in selected time window [J]

Names = fieldnames(Results.MActivation);
for n = 1:length(Names)
    
    % muscle parameters
    Fiso    = Misc.params(1,:)';
    lMo     = Misc.params(2,:)';
    
    % indexes for analysis (all but the last index)
    iSel = 1:length(Results.MActivation.genericMRS(1,:))-1;
    
    % get states and controls
    exc     = Results.MExcitation.(Names{n})(:,iSel);       % muscle excitation
    act     = Results.MActivation.(Names{n})(:,iSel);
    lMtilde = Results.lMtildeopt.(Names{n})(:,iSel);
    
    % contraction velocity of the muscle fibers
    % vMtilde is in the redudnancy solver the time derivative of lMtilde (i.e expressed in fiberlengths[m]/s)
    vMtilde = Results.vMtilde.(Names{n})(:,iSel);
    vM 		= Results.vMtilde.(Names{n})(:,iSel).*lMo; 	% filling in the fiber length [m]
    
    
    % force contractile element (F/l * F/v * a * Fiso)
    Fce     = Results.FMltilde.(Names{n})(:,iSel).*Results.FMvtilde.(Names{n})(:,iSel).*act.*Fiso;       % act * F/L * F/v* Fiso
    
    % passive muscle forces
    Fpass   = Results.Fpe.(Names{n})(:,iSel);
    
    % muscle mass
    volM = Fiso.*lMo;
    tension = getSpecificTensions_lr(Misc.MuscleNames_Input);
    musclemass = volM.*(1059.7)./(tension*1e6);
    
    % percentage slow and fast twitch fibers
    pctst = getSlowTwitchRatios_lr(Misc.MuscleNames_Input);
    
    % multiplier F/l curve
    FL_mult = Results.FMltilde.(Names{n})(:,iSel);
    
    % model mass
    modelmass = 72;
    
    % scale factor for tanh smoothing (can be very high here)
    b = 1000;
    
    % Barghava 2004
    E_Bargh = nan(size(act));
    E_Barh_model = nan(1,size(act,2));
    for i = iSel
        [E_Bargh(:,i),~,~,~,~,E_Barh_model(i)] = ...
            getMetabolicEnergySmooth2004all(exc(:,i),act(:,i),lMtilde(:,i),vM(:,i),Fce(:,i),Fpass(:,i),...
            musclemass,pctst,FL_mult(:,i),Fiso,modelmass,b);
    end
    
    % Umberger 2003
    E_Umb2003 = nan(size(act));
    E_Umb2003_model = nan(1,size(act,2));
    for i = iSel
        [E_Umb2003(:,i),~,~,~,E_Umb2003_model(i)] = ...
            getMetabolicEnergySmooth2003all(exc(:,i),act(:,i),lMtilde(:,i),vMtilde(:,i),vM(:,i),Fce(:,i),...
            musclemass,pctst,10,FL_mult(:,i),modelmass,b);
    end
    
    % Umberger 2010
    E_Umb2010 = nan(size(act));
    E_Umb2010_model = nan(1,size(act,2));
    for i = iSel
        [E_Umb2010(:,i),~,~,~,E_Umb2010_model(i)] = ...
            getMetabolicEnergySmooth2010all(exc(:,i),act(:,i),lMtilde(:,i),vMtilde(:,i),vM(:,i),Fce(:,i),...
            musclemass,pctst,10,FL_mult(:,i),modelmass,b);
    end
    
    % Uchida 2016
    E_Uch2016 = nan(size(act));
    E_Uch2016_model = nan(1,size(act,2));
    for i = iSel
        [E_Uch2016(:,i),~,~,~,E_Uch2016_model(i)] = ...
            getMetabolicEnergySmooth2016all(exc(:,i),act(:,i),lMtilde(:,i),vMtilde(:,i),vM(:,i),Fce(:,i),...
            musclemass,pctst,10,FL_mult(:,i),modelmass,b);
    end
    
    
    % store metabolic power in a structure
    Energy.(Names{n}).Bargh2004.Edot = E_Bargh;
    Energy.(Names{n}).Bargh2004.Edot_model= E_Barh_model;
    
    % store metabolic power in a structure
    Energy.(Names{n}).Umb2003.Edot = E_Umb2003;
    Energy.(Names{n}).Umb2003.Edot_model= E_Umb2003_model;
    
    % store metabolic power in a structure
    Energy.(Names{n}).Umb2010.Edot = E_Umb2010;
    Energy.(Names{n}).Umb2010.Edot_model= E_Umb2010_model;
    
    % store metabolic power in a structure
    Energy.(Names{n}).Uch2016.Edot = E_Uch2016;
    Energy.(Names{n}).Uch2016.Edot_model= E_Uch2016_model;
    
    % compute the metabolic work
    t = Results.Time.(Names{n})(1:end-1);
    Energy.(Names{n}).Bargh2004.W = trapz(t,E_Barh_model);
    Energy.(Names{n}).Umb2003.W = trapz(t,E_Umb2003_model);
    Energy.(Names{n}).Umb2010.W = trapz(t,E_Umb2010_model);
    Energy.(Names{n}).Uch2016.W = trapz(t,E_Uch2016_model);
end

end

