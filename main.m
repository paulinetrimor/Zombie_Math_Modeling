%Code created by: Pauline Trimor and Kendra Sands
%Main function, control function

decayTable = constructDecayTable();
N0 = input('Enter number of uranium atoms > '); %units in atoms
U238_lambda = 1.55e-10; %urainum decay rate found from ln(2)/halflife of uranium
U238_halflife = 4.468e9;
init_activity = N0*U238_lambda; %A= kN atoms/year (1.55e-7)
numSpecies = 21;

elem = 'U238';
new_elem = elem;
parentLambda = U238_lambda;
parentAtoms =N0;
parentMass =N0; %units in atoms
currentAtoms =N0; 
activity = init_activity; %units in disintergration/year
control = 1;

massLookup = containers.Map('KeyType', 'int32', 'ValueType', 'any');
massLookup(1) = 'U238';
massLookup(2) = 'Th234'; 
massLookup(3) = 'Pa234m'; 
massLookup(4) = 'Pa234'; 
massLookup(5) = 'U234'; 
massLookup(6) = 'Th230'; 
massLookup(7) = 'Ra226'; 
massLookup(8) = 'Rn222'; 
massLookup(9) = 'Po218'; 
massLookup(10) = 'At218'; 
massLookup(11) = 'Rn218'; 
massLookup(12) = 'Pb214'; 
massLookup(13) = 'Bi214'; 
massLookup(14) = 'Po214'; 
massLookup(15) = 'Tl210'; 
massLookup(16) = 'Pb210'; 
massLookup(17) = 'Bi210'; 
massLookup(18) = 'Po210'; 
massLookup(19) = 'Hg206'; 
massLookup(20) = 'Tl206'; 
massLookup(21) = 'Pb206'; 



nparts = ones(1,N0);
oparts = ones(1,N0);

tstep=1e6;
tvec=tstep:tstep:1e10;

% matrix containing number of nuclides per type per time
discreteNMatrix = zeros(numSpecies,length(tvec)+1);
discreteNMatrix(1,1) = N0;

%matrix for continuous analytic solution
continuousNMatrix = zeros(numSpecies, length(tvec) +1);
continuousNMatrix(1,1) = N0;

% matrix containing all activity values for each nuclide group
discreteAMatrix = zeros(numSpecies,length(tvec)+1);
discreteAMatrix(1,1) = init_activity;

% TODO: Continuous activity matrix

contNucType = 1; % nuclide 'type' for looking up nuclide key for continuous solution
tprev = 0; % previous time

i = 2;
for it=tvec % time to track
    
    discreteNMatrix(:,i) = discreteNMatrix(:,i-1);
    discreteAMatrix(:,i) = discreteAMatrix(:,i-1);
    continuousNMatrix(:,i) = continuousNMatrix(:,i-1);
    
    for nucType=1:(numSpecies-1)
        %for discrete calculations
        particles = find(nparts == nucType); % find all particles of type nucType
        if (length(particles) > 0) % particles of nucType exist
            nucKey = char(massLookup(nucType)); % get the nuclide's name
            record = decayTable(nucKey);
            hlf = record{4};
            
            % generate probability of decay
            prob_decay = (1-exp(-log(2)*tstep/hlf));
            
            % roll the dice for random numbers == to nuclides found
            numDecayed = rand(1,length(particles));
            
            lambda = 0; % scoped outside of Monte Carlo simulation
            % if anything decays...
            decayed = find(numDecayed < prob_decay);
            if (length(decayed) > 0)
                activity = discreteAMatrix(nucType, i);
                [nextNuclide, daughterActivity] = decay2(record, nucKey, activity, it);
                daughterRecord = decayTable(char(nextNuclide)); % get record of next nuclide
                daughterType = daughterRecord{1}; % get numerical type for next nuclide in chain
                activity = daughterActivity;
                
                oparts(particles(decayed)) = daughterType; % modify decayed particle's species
                % update table of # of nuclides
                discreteNMatrix(nucType, i) = discreteNMatrix(nucType, i) - length(decayed); % reduce # of this type by amount decayed
                discreteNMatrix(daughterType, i) = discreteNMatrix(daughterType,i) + length(decayed); % increase # of new type by amount decayed
                             
                % update table of activity values
                discreteAMatrix(daughterType, i) = daughterActivity;
            end
        end    
    end
    % for discrete solution
    nparts = oparts;
    
    % for continuous solution
    % compute number of nuclides from [tprev, it]
    contNucKey = char(massLookup(contNucType));
    record = decayTable(contNucKey); % current nuclide record
    daughterKey = getDaughterKey(record); % determine daughter(B) key  
    
    % continue if still able to decay
    if ~strcmp(char(daughterKey), 'stable')
        daughterRecord = decayTable(daughterKey); % get record for daughter
    
        N_A0 = continuousNMatrix(contNucType, i);
        N_B = Bateman(record, daughterRecord, N_A0, it);

        contNucType = contNucType + 1; % increment type / parent daughter    
        continuousNMatrix(contNucType, i) = N_B; % store num of daugter nuclide
    end

    tprev = it; % update previous time
    i = i + 1; % increment timestep iterator
end

%plot here
figure(1);
loglog(tvec,discreteNMatrix(1,1:end-1),'r','LineWidth',1); hold on;
plot(tvec,discreteNMatrix(2,1:end-1),'m','LineWidth',1);
plot(tvec,discreteNMatrix(3,1:end-1),'c','LineWidth',1);
plot(tvec,discreteNMatrix(4,1:end-1),'k','LineWidth',1);
plot(tvec,discreteNMatrix(5,1:end-1),'b','LineWidth',1);
plot(tvec,discreteNMatrix(6,1:end-1),'y','LineWidth',1); 
plot(tvec,discreteNMatrix(7,1:end-1), 'r', 'LineWidth', 1);
plot(tvec,discreteNMatrix(8,1:end-1), 'b', 'LineWidth', 1);
plot(tvec,discreteNMatrix(9,1:end-1), 'g', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(10,1:end-1), 'r', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(11,1:end-1), 'm', 'LineWidth', 1);
plot(tvec,discreteNMatrix(12,1:end-1), 'c', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(13,1:end-1), 'b', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(14,1:end-1), 'y', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(15,1:end-1), 'r', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(16,1:end-1), 'b', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(17,1:end-1), 'g', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(18,1:end-1), 'r', 'LineWidth', 1);
plot(tvec,discreteNMatrix(19,1:end-1), 'm', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(20,1:end-1), 'c', 'LineWidth', 1); 
plot(tvec,discreteNMatrix(21,1:end-1), 'g', 'LineWidth', 1); hold off;
legend('U238','Th234','Pa234m','Pa234','U234','Th230','Ra226','Rn222',...
    'Po218','At218','Rn218','Pb214','Bi214','Po214','Tl210','Pb210','Bi210',...
    'Po210','Hg206','Tl206','Pb206');
% h=get(g);
% set(h.Parent,'FontSize',16);
title('Monte Carlo - Number of Atoms');
xlabel('Time (years)','FontSize',16);
ylabel('Number of atoms','FontSize',16);

figure(2) 
loglog(tvec,discreteNMatrix(1,1:end-1),'r','LineWidth',1); hold on;
plot(tvec,discreteNMatrix(21,1:end-1),'b','LineWidth',1); hold off;
legend('U238','Pb206');
% h=get(g);
% set(h.Parent,'FontSize',16);
title('Monte Carlo- Number of Atoms');
xlabel('Time (years)','FontSize',16);
ylabel('Number of atoms','FontSize',16);
% 
figure(3)
loglog(tvec,discreteAMatrix(1,1:end-1),'r', 'LineWidth', 1); hold on;
plot(tvec,discreteAMatrix(2,1:end-1),'c','LineWidth', 1);
plot(tvec,discreteAMatrix(3,1:end-1),'k','LineWidth', 1);
plot(tvec,discreteAMatrix(4,1:end-1),'b','LineWidth', 1);
plot(tvec,discreteAMatrix(5,1:end-1),'y','LineWidth', 1);
plot(tvec,discreteAMatrix(6,1:end-1),'r','LineWidth',1);
plot(tvec,discreteAMatrix(7,1:end-1),'b','LineWidth', 1);
plot(tvec,discreteAMatrix(21,1:end-1),'g','LineWidth', 1);hold off;
title('Monte Carlo - Activity');
legend('U238', 'Th234', 'Pa234m', 'Pa234', 'U234', 'Th230', 'Ra226','Pb206');
xlabel('Time (years)');
ylabel('Activity (disintegrations/year)'); 

% continuous N plot
figure(4)
loglog(tvec, continuousNMatrix(1,1:end-1),'r','LineWidth', 1); hold on;
loglog(tvec,discreteNMatrix(1,1:end-1),'c','LineWidth', 1);
loglog(tvec,continuousNMatrix(21,1:end-1),'r','LineWidth', 1);
loglog(tvec,discreteNMatrix(21,1:end-1),'c','LineWidth', 1); hold off;
title('Continuous Solution - Nuclide Decay');
legend('continuous U238', 'discrete U238', 'continuous Pb206', 'discrete Pb206');
xlabel('Time (years)');
ylabel('Number of atoms'); 



































