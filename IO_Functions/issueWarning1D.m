function [Misc] = issueWarning1D(struct2check,structName,DatStore,Misc)
c = length(struct2check);
err_check = ones(1,c);
for t = 1:Misc.nTrials
    % get muscle names
    MNames = DatStore(t).MuscleNames;

    % optimization parameters muscles
    for j = 1:c
        if any(ismember(MNames,struct2check{j}))
            err_check(j) = 0;
        end
    end
end
if any(err_check)
    warning(['Could not find muscles ' struct2check{find(err_check)} ' in the selected dofs of any trial',...
        ', please adapt Misc.' structName '. Removing these muscles from Misc.' structName]);
    Misc.(structName) = struct2check(find(~err_check));
end    

