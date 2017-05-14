%Code created by: Pauline Trimor and Kendra Sands
%decaytable with all possible daughter nuclides and their properties


function [ decayTable ] = constructDecayTable()
% ConstructDecayTable 
    %   Constructs the decay table
    % Implementing decay rate, half-life, branching ratios, and parent/daughter
    % nuclides in a (Hash)Map data structure of type <char,any>. 
    % decayTable('char') = ['isotopeName', [decayModes], [branchingRatios], half-life in years, [daughterNuclides]]; 
    decayTable = containers.Map('KeyType', 'char', 'ValueType', 'any');
    decayTable('U238') = {1, {'alpha'}, {100.0}, 4.468e9, {'Th234'}};
    decayTable('Th234') = {2, {'beta'}, {100.0}, 0.06603, {'Pa234m'}}; 
    decayTable('Pa234m') = {3, {'IT','beta'}, {0.16, 99.84}, 5.2922e-5 , {'Pa234', 'U234'}}; 
    decayTable('Pa234') = {4, {'beta'}, {100.0}, 7.7546e-4, {'U234'}}; 
    decayTable('U234') = {5, {'alpha'}, {100.0}, 2.4550e5, {'Th230'}}; 
    decayTable('Th230') = {6, {'alpha'}, {100.0}, 7.5400e4, {'Ra226'}}; 
    decayTable('Ra226') = {7, {'alpha'}, {100.0}, 1600, {'Rn222'}}; 
    decayTable('Rn222') = {8, {'alpha'}, {100.0}, 0.0105, {'Po218'}}; 
    decayTable('Po218') = {9, {'beta','alpha'}, {0.02, 99.98}, 5.8942e-6, {'At218', 'Pb214'}}; 
    decayTable('At218') = {10, {'beta', 'alpha'}, {0.1, 99.9}, 4.7565e-8, {'Rn218', 'Bi214'}}; 
    decayTable('Rn218') = {11, {'alpha'}, {100.0}, 1.1098e-9, {'Po214'}}; 
    decayTable('Pb214') = {12, {'beta'}, {100.0}, 5.0989e-5, {'Bi214'}}; 
    decayTable('Bi214') = {13, {'alpha', 'beta'}, {0.021, 99.979}, 3.7861e-5, {'Tl210', 'Po214'}}; % TODO: Po214 or Po241? 
    decayTable('Po214') = {14, {'alpha'}, {100.0}, 5.2099e-12, {'Pb210'}}; 
    decayTable('Tl210') = {15, {'beta'}, {100.0}, 2.4734e-6, {'Pb210'}}; 
    decayTable('Pb210') = {16, {'alpha', 'beta'}, {0.0000019, 99.999999}, 22.2, {'Hg206', 'Bi210'}}; 
    decayTable('Bi210') = {17, {'alpha', 'beta'}, {0.0000132, 99.9999868}, 5.7215e-4, {'Tl206', 'Po210'}}; 
    decayTable('Po210') = {18, {'alpha'}, {100.0}, 0.3791, {'Pb206'}}; 
    decayTable('Hg206') = {19, {'beta'}, {100.0}, 1.5829e-5, {'Tl206'}}; 
    decayTable('Tl206') = {20, {'beta'}, {100.0}, 7.9947e-6, {'Pb206'}}; 
    decayTable('Pb206') = {21, {'stable'}, {0}, 0, {'stable'}}; 
end

