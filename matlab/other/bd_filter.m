% Utility function for ThreeColorDat2States.m
% Jan 8, 2023 X. Feng
%
% Input: y_state - 1D list with values 0, 1, 2, or 3
% Output: new_list - 1D list with values 0, 1, 2, or 3
% This function finds the longest sub-list not interrupted by 3s
% and sets other elements to 0

function new_list = bd_filter(y_state)

longest_list_start = 1;
longest_list_end = -1;
longest_list_len = 0;

curr_list_start = 1;
curr_list_end = -1;
curr_list_len = 0;

for i = 1: length(y_state)
    
    s = y_state(i);
    
    if s == 3
        
        curr_list_start = i + 1;
        curr_list_end = i;
        curr_list_len = 0;
        continue
        
    else
        
        curr_list_len = curr_list_len + 1;
        curr_list_end = curr_list_end + 1;
        
        if curr_list_len > longest_list_len
            
            % Update longest list
            
            longest_list_start = curr_list_start;
            longest_list_end = curr_list_end;
            longest_list_len = curr_list_len;
        end
    end
end

% set list elements outside of the longest sub-list to 3

new_list = y_state;

for j = 1 : length(y_state)
    
    if (j < longest_list_start) || (j > longest_list_end)
        new_list(j) = 3;
    end
    
end

end