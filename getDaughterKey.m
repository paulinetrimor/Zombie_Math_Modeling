% this function takes in the current nuclide's record :: {1x5} cell array
function [daughterKey] = getDaughterKey(record)
    decaymodes = length(record{2}); % number of possible decay modes/daughter nuclides
    % determine daughter nuclide
    if (decaymodes > 1) % nuclide has multiple daughter nuclides/decay modes
        ratioDecider = rand(1);
        modeA = record{3}(1); 
        % modeB = elemArray{3}(2); % no need to look up 2nd mode
        if (str2double(modeA) <= ratioDecider)
            daughterKey = char(record{5}{1});
        else
            daughterKey = char(record{5}(2));
        end
    else % nuclide has single decay mode
        daughterKey = char(record{5}); %decay to new element
    end
end