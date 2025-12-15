%% Heuristically finds the maximum number of non-lethal knockouts
% essentially, I will make a 5000 vecotrs of random permutations of
% model.genes that has 10 genes in it. Then I will test to see is all of
% those are feasible. I will filter for the ones that are, then from that I
% will compare for simularities and union the ones that are most simular


load Ec_iJO1366.mat; load r1_num_knockouts.mat;

emodel.A=model.S; emodel.obj= model.c; emodel.sense = char('='); emodel.modelsense = 'max';
eparams.OutputFlag=0;

load lengths.mat; load maxlength.mat; load longestKO.mat;

lengths = lengths;
maxlength = maxlength;
vg = 1;
longestKO = longestKO;
geneVector = [];

k=1;
for k=1:10000

    geneVector = [];
    vg=1;

    while abs(vg)>0.0001
        newGene = r1_num_knockouts(randi(length(r1_num_knockouts)));

        if any(geneVector == newGene)
            % If newgene is found, skip to the next iteration and not add to
            continue;
        else
            geneVector = [geneVector;newGene];
        end

        emodel.lb=model.lb;
        emodel.ub=model.ub;

        rxnList = [];
        for i=1:length(geneVector)
            rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
        end
        rxnList = sort(rxnList);
        if (~isempty(rxnList))
            x = true(size(model.genes));
            x(geneVector) = false;
            removeList{newGene} = [];
            for i = 1:length(rxnList)
                if (~eval(model.rules{rxnList(i)}))
                    removeList{newGene} = union(removeList{newGene},rxnList(i));
                end
            end
        end

        for j = 1:length(rxnList)
            emodel.lb(rxnList(j)) = 0;
            emodel.ub(rxnList(j)) = 0;
        end

        results = gurobi(emodel,eparams);
        if strcmp(results.status, 'INFEASIBLE')
            vg=0;
            continue
        else
            vg=results.objval;
        end




    end

    currentLength = length(geneVector)-1;
    lengths = [lengths;currentLength];
    if currentLength>maxlength
        maxlength = currentLength;
        longestKO = geneVector(1:currentLength,:);
        save 'maxlength' maxlength
        save 'longestKO.mat' longestKO
        maxlength

    end
    save 'lengths.mat' lengths
end

maxlength
% Ensuring correctness

emodel.lb=model.lb;
emodel.ub=model.ub;

rxnList = [];
geneVector = longestKO;
newGene=geneVector(end);
for i=1:length(geneVector)
    rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
end
rxnList = sort(rxnList);
if (~isempty(rxnList))
    x = true(size(model.genes));
    x(geneVector) = false;
    removeList{newGene} = [];
    for i = 1:length(rxnList)
        if (~eval(model.rules{rxnList(i)}))
            removeList{newGene} = union(removeList{newGene},rxnList(i));
        end
    end
end

for j = 1:length(rxnList)
    emodel.lb(rxnList(j)) = 0;
    emodel.ub(rxnList(j)) = 0;
end

results = gurobi(emodel,eparams);
if strcmp(results.status, 'INFEASIBLE')
    vg=0;
else
    vg=results.objval;
end

vg


