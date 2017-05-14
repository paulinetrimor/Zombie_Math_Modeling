%Code Created by: Pauline Trimor and Kendra Sands
%Radioactive Uraninum Project

%decay function takes in nuclide and searches decayTable for values
function[new_elem, newParentAtoms, daughterAtoms, newParentLambda, new_act]=decay(decayTable, elem, currentAtoms, parentAtoms, parentLambda, activity, time)
    if strcmp(char(elem), 'stable')
        return;
    end
    elemArray = decayTable(char(elem));
    decaymodes = length(elemArray{2});
    
    % print current element
    fprintf(strcat('Element: ', char(elem)));
    
    % print mode(s) % corresponding branching ratio(s) of current element
    num_modes=length(elemArray{2});
    if num_modes == 1
        fprintf(strcat('\nDecay: ', elemArray{5}, ' \tMode: ', elemArray{2}, '\tRatio: ', string(elemArray{3}), '%%\tHalf-life: ', string(elemArray{4}), ' years'));
        %fprintf('Decay: ' + elemArray{5} + ' ' + elemArray{2} + ' (' + elemArray{3} + '%)');
    elseif num_modes > 1
        modes = strcat('\nDecay:', elemArray{5}(1), ' \tMode: ', elemArray{2}(1), '\tRatio: ',string(elemArray{3}(1)), '%%\tHalf-Life: ', string(elemArray{4}), ' years');
        fprintf(modes);
        for i=2:num_modes
            modes = strcat('\nDecay: ', elemArray{5}(i) , ' \tMode: ' , elemArray{2}(i), '\tRatio: ', string(elemArray{3}(i)), '%%\tHalf-Life: ', string(elemArray{4}), ' years');
            fprintf(modes);
        end
    end
    
    halflife = elemArray{4};
    %determines new decay rate for current isotope
    lambda = (log(2))/halflife;
    fprintf('\nDecay Constant (lambda): ');
    fprintf('%.4e', lambda');
    
    %time is changing in years
    %update current activity
    new_act = activity*exp(-lambda*time); %units of decays per year
    
    % print new activity
    fprintf('\nNew Activity: ');
    fprintf('%.4e',new_act);
    
    %need to account for time here such that
    rp = rand(1);
    prob_decay = (1-exp(-log(2)/halflife));
    %if (rp >= prob_decay) % nuclide decays * TODO: resolve decay probability
    if (rp >= 0)
        if (decaymodes > 1) % has more than 1 decay mode / ratio
            ratioDecider = rand(1);
            modeA = elemArray{3}(1);
            modeB = elemArray{3}(2);
            if (str2double(modeA) <= ratioDecider)
                new_elem = elemArray{5}{1};
            else
                new_elem = elemArray{5}(2);
            end
        else % single decay mode
            new_elem = elemArray{5}; %decay to new element
        end
        daughterAtoms = (parentLambda / (lambda-parentLambda))*parentAtoms*(exp(-parentLambda*time) - exp(-lambda*time));
        %newParentAtoms = parentAtoms * exp(-parentLambda * time);
        newParentAtoms = parentAtoms - daughterAtoms;
        newParentLambda = lambda;
        %update # of atoms left
     
        % do we reduce the source's atoms?
    %elseif (rp < prob_decay) 
    else % nuclide does not become new nuclide    
        new_elem = elem; % remain the same
        daughterAtoms = parentAtoms*exp(-lambda*time);
        newParentAtoms = parentAtoms * exp(-parentLambda * time);
        newParentLambda = lambda;
    end
    
    fprintf('\nNumber of parent element atoms: ');
    fprintf('%.4e',currentAtoms);
    
    fprintf(strcat('\nNumber of new element atoms (', char(new_elem), '): '));
    fprintf('%.4e',daughterAtoms);
    fprintf('\n\n');
    parentLambda = lambda; 
end

% % compute the order of magnitude of a number
% function [output] = magnitude(num)    
%     i = 0;
%     num = abs(num);
%     if (num > 1) % larger number
%         while num / 10 >= 1
%             num = num / 10;
%             i = i + 1;
%         end
%     else % [0,1]
%         while (num * 10 <= 1)
%             num = num * 10;
%             i = i - 1;
%         end
%     end
%     output = i;
% end

























