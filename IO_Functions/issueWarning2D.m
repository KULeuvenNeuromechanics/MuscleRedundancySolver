function [Misc] = issueWarning2D(struct2check,structName,DatStore,Misc)

[r,c] = size(struct2check);
if c>0
    if c~=2
        error(['struct2check does not have exactly 2 columns',...
            ', please adapt Misc.' structName])
    else
        err_cfl = ones(size(struct2check));
        for t = 1:Misc.nTrials
            % get muscle names
            MNames = DatStore(t).MuscleNames;

            for i = 1:r
                e1 = any(ismember(MNames,struct2check{i,1}));
                e2 = any(ismember(MNames,struct2check{i,2}));
                if e1
                    err_cfl(i,1) = 0;
                end
                if e2
                    err_cfl(i,2) = 0;
                end
                if e1 && ~e2
                    warning(['Could not find muscle' struct2check{i,2} ' coupled with ' struct2check{i,1} ' in model of trial' Misc.trialName{t} ', please adapt Misc.' structName])
                elseif ~e1 && e2
                    warning(['Could not find muscle' struct2check{i,1} ' coupled with ' struct2check{i,2} ' in model of trial' Misc.trialName{t} ', please adapt Misc.' structName])
                end
            end
        end
    end
    if any(any(err_cfl))
        warning(['Could not find muscles ' struct2check{find(err_cfl)} ' in the selected dofs of any trial',...
            ', please adapt Misc.' structName '. Removing this muscle and its couple from Misc.' structName]);
        Misc.(structName) = struct2check(~any(err_cfl,2),:);
    end
end

end
