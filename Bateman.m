% this function takes in the current nuclide's record :: {1x5} cell array
function [daughterNuclideNum] = Bateman(parentRecord, daughterRecord, parentNuclideNum, t)
    parentLambda = log(2) / parentRecord{4};
    daughterLambda = log(2) / daughterRecord{4};
    N_A = parentNuclideNum;
    lA = parentLambda;
    lB = daughterLambda;
    
    %N_B=N_A0*(lambdaA/(lambdaB-lambdaA))*[e^(-lambdaA*t)-e^(-lambdaB*t)] + N_B0*e^(-lambdaB*t)
    % where lambdaB == daughterLambda && lambdaA == parentLambda
%     exp1 = abs(-lA * t);
%     if (exp1 > 4.3e3)
%         fprintf('%.e4', exp1);
%         % e^-3*t == 10^(-3*2*log10(e))
%         exp1_logd = (-lA*t)*log10(exp(1))
%         fprintf('%.e4', exp1_logd);
%     end
    N_B = N_A * (lA/(lB - lA)) * (exp(-lA * t) - exp(-lB * t));
    daughterNuclideNum = N_B;
end