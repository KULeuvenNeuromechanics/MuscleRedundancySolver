

%% passive force length properties


%% vary kSF

lMtilde = 0.4:0.01:1.6;
figure();
iV = 0.8:0.01:1.2;
Cols = copper(length(iV));
ct = 1;
for i = iV
    kpe = 4;
    kSF = i;

    % Parameters of passive muscle force-length characteristic
    e0 = 0.6;
    t50 = exp(kpe .* (0.2 - 0.10e1) ./ e0);
    pp1 = (t50 - 0.10e1);
    t7 = exp(kpe);
    pp2 = (t7 - 0.10e1);
    Fpparam = [pp1;pp2];

    % compute the % of the maximal force generating capacity at the current
    % state

    % get passive force
    e0 = 0.6;
    t5 = exp(kpe .* (lMtilde - kSF) ./ e0);
    Fpe = ((t5 - 0.10e1) - Fpparam(1)) ./ Fpparam(2);

    plot(lMtilde,Fpe,'Color',Cols(ct,:)); hold on;
    ct = ct+1;
end




%% vary kpe

lMtilde = 0.4:0.01:1.6;
figure();
iV = 2:0.1:6;
Cols = copper(length(iV));
ct = 1;
for i = iV
    kpe = i;
    kSF = 1.1;

    % Parameters of passive muscle force-length characteristic
    e0 = 0.3;
    t50 = exp(kpe .* (0.2 - 0.10e1) ./ e0);
    pp1 = (t50 - 0.10e1);
    t7 = exp(kpe);
    pp2 = (t7 - 0.10e1);
    Fpparam = [pp1;pp2];

    % compute the % of the maximal force generating capacity at the current
    % state

    % get passive force
    e0 = 0.6;
    t5 = exp(kpe .* (lMtilde - kSF) ./ e0);
    Fpe = ((t5 - 0.10e1) - Fpparam(1)) ./ Fpparam(2);

    plot(lMtilde,Fpe,'Color',Cols(ct,:)); hold on;
    ct = ct+1;
end

%% vary e0

lMtilde = 0.4:0.01:1.6;
figure();
e0V = 0.2:0.01:1;
Cols = copper(length(e0V));
ct = 1;
for i = e0V
    kpe = 4;
    kSF = 1;
    e0 = i;

    % Parameters of passive muscle force-length characteristic
    t50 = exp(kpe .* (0.2 - 0.10e1) ./ e0);
    pp1 = (t50 - 0.10e1);
    t7 = exp(kpe);
    pp2 = (t7 - 0.10e1);
    Fpparam = [pp1;pp2];

    % compute the % of the maximal force generating capacity at the current
    % state

    % get passive force
    t5 = exp(kpe .* (lMtilde - kSF) ./ e0);
    Fpe = ((t5 - 0.10e1) - Fpparam(1)) ./ Fpparam(2);

    plot(lMtilde,Fpe,'Color',Cols(ct,:)); hold on;
    ct = ct+1;
end


%% Test function implementation

nMuscles = 43;
lMtilde = ones(nMuscles,1).*1.4;
vMtilde = 0;
vMtildemax = 10;
kpe = ones(nMuscles,1)*4 + randn(nMuscles,1);
kSF = ones(nMuscles,1);
[Fpe,FMltilde,FMvtilde] = getForceLengthVelocityProperties_setPassiveParam(lMtilde,vMtilde,vMtildemax,kpe,kSF);

figure();
plot(kpe,Fpe,'o')
