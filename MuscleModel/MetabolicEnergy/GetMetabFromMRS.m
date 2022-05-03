function [Energy] = GetMetabFromMRS(Results,Misc,DatStore,modelmass)
% GetMetabFromMR Computes metabolic energy consumption from the muscle
% redundancy solver results.
% input arguments:
%   Results: matlab structure "Results" (output argument of the MRS)
%   Misc: matlab structure "Misc" (output argument of the MRS)
%   modelmass: mass of the musculoskeltal model (in kg)
% output arguments:
%   Energy 	... Edot: metabolic power without basal rate for each muscle [W]
%           ... Edot_model: sum of metabolic power of all selected muscles with basal rate [W]
%%
Names = fieldnames(Results.MActivation);
for n = 1:length(Names)
    for tr = 1:Misc.nTrials
        
        % muscle parameters
        Fiso    = Misc.params(1,Misc.idx_allMuscleList{tr})';
        lMo     = Misc.params(2,Misc.idx_allMuscleList{tr})';
        
        % indexes for analysis (all but the last index)
        iSel = 1:length(Results.MActivation(tr).(Names{n})(1,:))-1;
        
        % get states and controls
        exc     = Results.MExcitation(tr).(Names{n})(:,iSel);       % muscle excitation
        act     = Results.MActivation(tr).(Names{n})(:,iSel);
        lMtilde = Results.lMtildeopt(tr).(Names{n})(:,iSel);
        
        % contraction velocity of the muscle fibers
        % vMtilde is in the redudnancy solver the time derivative of
        % lMtilde. (note that this might be confusing !)
        vMtilde = Results.vMtilde(tr).(Names{n})(:,iSel);
        vM 		= Results.vMtilde(tr).(Names{n})(:,iSel).*lMo;        
        
        % force contractile element (F/l * F/v * a * Fiso)
        Fce     = Results.FMltilde(tr).(Names{n})(:,iSel).*Results.FMvtilde.(Names{n})(:,iSel).*act.*Fiso;
        
        % passive muscle forces
        Fpass   = Results.Fpe(tr).(Names{n})(:,iSel);
        
        % muscle mass
        volM = Fiso.*lMo;
        tension = getSpecificTensions_lr(DatStore(tr).MuscleNames); % you can add specific tensions for muscles in this file
        musclemass = volM.*(1059.7)./(tension*1e6);
        
        % percentage slow and fast twitch fibers
        pctst = getSlowTwitchRatios_lr(DatStore(tr).MuscleNames); % you can add fiber type distributions here
        
        % multiplier F/l curve
        FL_mult = Results.FMltilde(tr).(Names{n})(:,iSel);
            
        % scale factor for tanh smoothing (can be very high here)
        % remove tanh implementation here? Is not needed because the metabolic
        % energy equations are not used in the objective function here.
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
        Energy(tr).(Names{n}).Bargh2004.Edot = E_Bargh;
        Energy(tr).(Names{n}).Bargh2004.Edot_model= E_Barh_model;
        
        % store metabolic power in a structure
        Energy(tr).(Names{n}).Umb2003.Edot = E_Umb2003;
        Energy(tr).(Names{n}).Umb2003.Edot_model= E_Umb2003_model;
        
        % store metabolic power in a structure
        Energy(tr).(Names{n}).Umb2010.Edot = E_Umb2010;
        Energy(tr).(Names{n}).Umb2010.Edot_model= E_Umb2010_model;
        
        % store metabolic power in a structure
        Energy(tr).(Names{n}).Uch2016.Edot = E_Uch2016;
        Energy(tr).(Names{n}).Uch2016.Edot_model= E_Uch2016_model;
    end
end
end

