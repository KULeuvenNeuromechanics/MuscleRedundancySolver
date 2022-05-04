function [Misc] = issueWarning2D(struct2check,structName,DatStore,Misc)
% --------------------------------------------------------------------------
%issueWarning2D
%     Generate warnings should there be inconsistensies in the user-provided
%     muscles that are coupled in some way
% 
% INPUT:
%     struct2check
%     The field to be checked
% 
%     structName
%     Name of the field
%     
%     DatStore
%     Structure of all data
%     
%     Misc
%     Miscellaneous info used through the code
%     
% OUTPUT:
%     Misc
%     Miscellaneous info used through the code
%     
% Original author: Dhruv Gupta
% Original date: May 3, 2022
%
% Last edit by: Dhruv Gupta
% Last edit date: May 3, 2022
% --------------------------------------------------------------------------

[r,c] = size(struct2check);
% need input in as n x 2 matrix, with n representing the numer of pairs
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
                % check if both muscles are in each trial
                if e1 && ~e2
                    warning(['Could not find muscle' struct2check{i,2} ' coupled with ' struct2check{i,1} ' in model of trial_' num2str(t) ', please adapt Misc.' structName])
                elseif ~e1 && e2
                    warning(['Could not find muscle' struct2check{i,1} ' coupled with ' struct2check{i,2} ' in model of trial_' num2str(t) ', please adapt Misc.' structName])
                end
            end
        end
    end
    % check if a muscle is not there in any trial. If it is not there, both
    % that muscle and its couple are removed.
    if any(any(err_cfl))
        warning(['Could not find muscles ' struct2check{find(err_cfl)} ' in the selected dofs of any trial',...
            ', please adapt Misc.' structName '. Removing this muscle and its couple from Misc.' structName]);
        Misc.(structName) = struct2check(~any(err_cfl,2),:);
    end
end

end
