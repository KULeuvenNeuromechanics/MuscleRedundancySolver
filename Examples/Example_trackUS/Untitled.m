    a_opt = sol.value(a);
    lMtilde_opt = sol.value(lMtilde);
    % Muscle excitations
    e_opt = sol.value(e);
    % Reserve actuators
    aT_opt = sol.value(aT);
    % Time derivatives of muscle-tendon forces
    vMtilde_opt = sol.value(vMtilde);
    
    % Save results
    Ntot = 0;
    for trial = 1:nTrials
        t0 = DatStore(trial).time(1); tf = DatStore(trial).time(end);
        N = round((tf-t0)*Misc.Mesh_Frequency);
        Ntot = Ntot + N;
        % Time grid
        tgrid = linspace(t0,tf,N+1)';
        % Save results
        Time(trial).genericMRS = tgrid;
        MActivation(trial).genericMRS = a_opt(:,(trial-1)*(Ntot + trial - 1) + 1:(trial-1)*(Ntot + trial - 1) + N + 1);
        lMtildeopt(trial).genericMRS = lMtilde_opt(:,(trial-1)*(Ntot + trial - 1) + 1:(trial-1)*(Ntot + trial - 1) + N + 1);
        lM(trial).genericMRS = lMtilde_opt(:,(trial-1)*(Ntot + trial - 1) + 1:(trial-1)*(Ntot + trial - 1) + N + 1).*repmat(Misc.lOpt',1,length(tgrid));
        MvMtilde(trial).genericMRS = vMtilde_opt(:,Ntot*(trial-1) + 1:Ntot*(trial-1) + N);
        MExcitation(trial).genericMRS = e_opt(:,Ntot*(trial-1) + 1:Ntot*(trial-1) + N);
        RActivation(trial).genericMRS = aT_opt(:,Ntot*(trial-1) + 1:Ntot*(trial-1) + N)*Misc.Topt;
        MuscleNames = DatStore.MuscleNames;
        OptInfo = output;
        % Tendon forces from lMtilde
        lMTinterp(trial).genericMRS = DatStore(trial).LMTinterp;
        [TForcetilde_,TForce_] = TendonForce_lMtilde(lMtildeopt(trial).genericMRS',Misc.params,lMTinterp(trial).genericMRS,Misc.Atendon,Misc.shift);
        TForcetilde(trial).genericMRS = TForcetilde_';
        TForce(trial).genericMRS = TForce_';
        Ntot = Ntot + N;
    end