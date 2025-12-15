%% Maximum number of Knockouts exhaustive code
%
% Each round should return a number of knockout round columns vector
% with each row being a list of knocked out genes that did NOT kill the
% E. Coli
% Uses Gurobi to optimize



clear all; clc;

% run C:\gurobi1300\win64\matlab\gurobi_setup;
load Ec_iJO1366.mat

emodel.A=model.S; 
emodel.obj= model.c; 
emodel.sense = char('='); 
eparams.OutputFlag=0; 
emodel.modelsense = 'max';

survivedKnockouts = [];
num_survivedKnockout = [];
essential_genes = [];


% Sets up timing 
fprintf('Starting for loops...\n');

outerLoopTotal = 1; 
innerLoopTotal = length(model.genes)-1; 
totalOverallIterations = outerLoopTotal * innerLoopTotal;

% these can be changed if you want more/less updates or ant to change the
% size of the progress bar
updateInterval = 200;   
progressBarLength = 40;  

startTime = tic; % Start the timer for nested for loops

for numGene=1:length(model.genes)-1

    geneVector = [numGene];
    rxnList = [];
    for i=1:length(geneVector)
        rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
    end
    rxnList = sort(rxnList);
    if (~isempty(rxnList))
        x = true(size(model.genes));
        x(geneVector) = false;
        removeList{numGene} = [];
        for i = 1:length(rxnList)
            if (~eval(model.rules{rxnList(i)}))
                removeList{numGene} = union(removeList{numGene},rxnList(i));
            end
        end
    end

    emodel.lb=model.lb;
    emodel.ub=model.ub;

    for j = 1:length(rxnList)
        emodel.lb(rxnList(j)) = 0;
        emodel.ub(rxnList(j)) = 0;
    end

    results = gurobi(emodel,eparams);

    if strcmp(results.status, 'INFEASIBLE')
        essential_genes = cat(1,essential_genes,model.genes(numGene));
        continue
    else
        vg=results.objval;
        if abs(vg)>.0001
            survivedKnockouts = cat(1,survivedKnockouts,model.genes(numGene));
            num_survivedKnockout = cat(1,num_survivedKnockout,numGene);
        else
            essential_genes = cat(1,essential_genes,model.genes(numGene));
        end
    end


    % Reporting progress and completion time estimates
    currentOverallIteration = ((1 - 1) * innerLoopTotal) + numGene;
    if mod(currentOverallIteration, updateInterval) == 1 || currentOverallIteration == totalOverallIterations
        updateProgressDisplay(currentOverallIteration, totalOverallIterations, startTime, progressBarLength);
    end

end

fprintf("\n Complete!\n")

r1_knockouts = survivedKnockouts;
r1_num_knockouts = num_survivedKnockout;

save 'r1_knockouts.mat' r1_knockouts
save 'r1_num_knockouts.mat' r1_num_knockouts

%% Round 2 Knockouts

clear all;

load Ec_iJO1366.mat; load r1_knockouts.mat; load r1_num_knockouts.mat;

emodel.A=model.S; emodel.obj= model.c; emodel.sense = char('='); emodel.modelsense = 'max';
eparams.OutputFlag=0; 

survivedKnockouts = [];
num_survivedKnockout = [];
essential_genes2 = [];

% Sets up timing 
fprintf('Starting for loops...\n');

outerLoopTotal = length(r1_knockouts); 
innerLoopTotal = length(r1_knockouts); 
totalOverallIterations = outerLoopTotal*(outerLoopTotal-1)/2;

% these can be changed if you want more/less updates or ant to change the
% size of the progress bar
updateInterval = 1000;   
progressBarLength = 40;  

startTime = tic; % Start the timer for nested for loops

for k=1:outerLoopTotal
    for numGene=k+1:innerLoopTotal

        emodel.lb=model.lb;
        emodel.ub=model.ub;

        geneVector = [r1_num_knockouts(k,:) ,r1_num_knockouts(numGene)];
        rxnList = [];
        for i=1:length(geneVector)
            rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
        end
        rxnList = sort(rxnList);
        if (~isempty(rxnList))
            x = true(size(model.genes));
            x(geneVector) = false;
            removeList{numGene} = [];
            for i = 1:length(rxnList)
                if (~eval(model.rules{rxnList(i)}))
                    removeList{numGene} = union(removeList{numGene},rxnList(i));
                end
            end
        end

        for j = 1:length(rxnList)
            emodel.lb(rxnList(j)) = 0;
            emodel.ub(rxnList(j)) = 0;
        end

        results = gurobi(emodel,eparams);
        if strcmp(results.status, 'INFEASIBLE')
            essential_genes2 = cat(1,essential_genes2,model.genes(geneVector)');
            continue
        else
            vg=results.objval;
            if abs(vg)>.0001
                survivedKnockouts = cat(1,survivedKnockouts,model.genes(geneVector)');
                num_survivedKnockout = cat(1,num_survivedKnockout,geneVector);
            else
                essential_genes2 = cat(1,essential_genes2,model.genes(geneVector)');
            end
        end

    % Reporting progress and completion time estimates
    currentOverallIteration = ((k - 1) * innerLoopTotal) + numGene;
    if mod(currentOverallIteration, updateInterval) == 1 || currentOverallIteration == totalOverallIterations
        updateProgressDisplay(currentOverallIteration, totalOverallIterations, startTime, progressBarLength);
    end

    end
end

fprintf("\n Complete!\n")


temp_num_survivedKnockout = num_survivedKnockout;
temp_survivedKnockouts = survivedKnockouts;

r2_knockoutsSave = temp_survivedKnockouts;
r2_num_knockoutsSave = temp_num_survivedKnockout;
essential_genes2Save=essential_genes2;

save 'r2_knockouts_save.mat' r2_knockoutsSave
save 'r2_num_knockouts_save.mat' r2_num_knockoutsSave
save 'essential_genes2_save.mat' essential_genes2Save

% Round 3 Knockouts


fprintf("Round 3 Start")


load Ec_iJO1366.mat; 
load r1_knockouts.mat; load r1_num_knockouts.mat;
load r2_knockouts_save.mat; load r2_num_knockouts_save.mat;

emodel.A=model.S; emodel.obj= model.c; emodel.sense = char('='); emodel.modelsense = 'max';
eparams.OutputFlag=0; 

survivedKnockouts = [];
num_survivedKnockout = [];
essential_genes3 = [];

% Sets up timing 
fprintf('Starting for loops...\n');

outerLoopTotal = length(r2_knockoutsSave); 
innerLoopTotal = length(r1_knockouts); 
totalOverallIterations = outerLoopTotal * innerLoopTotal;

% these can be changed if you want more/less updates or ant to change the
% size of the progress bar
updateInterval = 1000;   
progressBarLength = 40;  

startTime = tic; % Start the timer for nested for loops

for k=1:outerLoopTotal
    for numGene=1:innerLoopTotal

        if any(r2_num_knockoutsSave(k, :) == numGene)
            % If numGene is found, skip to the next iteration
            continue;
        end

        emodel.lb=model.lb;
        emodel.ub=model.ub;

        geneVector = [r2_num_knockoutsSave(k,:) ,r1_num_knockouts(numGene)];
        rxnList = [];
        for i=1:length(geneVector)
            rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
        end
        rxnList = sort(rxnList);
        if (~isempty(rxnList))
            x = true(size(model.genes));
            x(geneVector) = false;
            removeList{numGene} = [];
            for i = 1:length(rxnList)
                if (~eval(model.rules{rxnList(i)}))
                    removeList{numGene} = union(removeList{numGene},rxnList(i));
                end
            end
        end

        for j = 1:length(rxnList)
            emodel.lb(rxnList(j)) = 0;
            emodel.ub(rxnList(j)) = 0;
        end

        results = gurobi(emodel,eparams);
        if strcmp(results.status, 'INFEASIBLE')
            essential_genes3 = cat(1,essential_genes3,model.genes(geneVector)');
            continue
        else
            vg=results.objval;
            if abs(vg)>.0001
                survivedKnockouts = cat(1,survivedKnockouts,model.genes(geneVector)');
                num_survivedKnockout = cat(1,num_survivedKnockout,geneVector);
            else
                essential_genes3 = cat(1,essential_genes3,model.genes(geneVector)');

            end
        end

        % Reporting progress and completion time estimates
        currentOverallIteration = ((k - 1) * innerLoopTotal) + numGene;
        if mod(currentOverallIteration, updateInterval) == 1 || currentOverallIteration == totalOverallIterations
            updateProgressDisplay(currentOverallIteration, totalOverallIterations, startTime, progressBarLength);
        end

    end
end

fprintf("\n Complete!\n")

r3_knockouts = survivedKnockouts;
r3_num_knockouts = num_survivedKnockout;

save 'r3_knockouts.mat' r3_knockouts
save 'r3_num_knockouts.mat' r3_num_knockouts
save 'essential_genes3.mat' essential_genes3



















