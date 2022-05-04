function [Misc] = issueWarning1D(struct2check,structName,DatStore,Misc)
% --------------------------------------------------------------------------
%issueWarning1D
%     Generate warnings should there be inconsistensies in the user-provided
%     muscle names for fields with 1D list of muscles.
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

% check if muscle is in each trial
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

% check if muscle is not in any trial
if any(err_check)
    warning(['Could not find muscles ' struct2check{find(err_check)} ' in the selected dofs of any trial',...
        ', please adapt Misc.' structName '. Removing these muscles from Misc.' structName]);
    disp(' ')
    Misc.(structName) = struct2check(find(~err_check));
end    

