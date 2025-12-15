%% Maximum number of Knockouts base code
%
% Each round should return a i (numberof knockout round) columns vector
% with each row being a list of knocked out genes that did NOT kill the
% E. Coli
% uses MATLAB's linprog to optimize













%% Round 1 Knockouts
clear all;
options = optimoptions('linprog','Display','none');

load Ec_iJO1366.mat

survivedKnockouts = [];
num_survivedKnockout = [];

for numGene=1:length(model.genes)-1

    load Ec_iJO1366.mat

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

    for j = 1:length(rxnList)
        model.lb(rxnList(j)) = 0;
        model.ub(rxnList(j)) = 0;
    end

    [v,vg,exitFlag] = linprog(-model.c,[],[],model.S,zeros([1,1805]),model.lb,model.ub,options);

    if exitFlag==-2
        continue
    end

    if abs(vg)>.0001
        survivedKnockouts = cat(1,survivedKnockouts,model.genes(numGene));
        num_survivedKnockout = cat(1,num_survivedKnockout,numGene);
    end
    if mod(numGene,100)==0
        progress = numGene/(length(model.genes));
        fprintf(string(progress))
    end

end

fprintf("\n Complete!\n")

r1_knockouts = survivedKnockouts;
r1_num_knockouts = num_survivedKnockout;

save 'r1_knockouts.mat' r1_knockouts
save 'r1_num_knockouts.mat' r1_num_knockouts

%% Round 2 Knockouts

clear all;
options = optimoptions('linprog','Display','none');

load Ec_iJO1366.mat; load r1_knockouts.mat; load r1_num_knockouts.mat;

survivedKnockouts = [];
num_survivedKnockout = [];
count=0;

for k=1:length(r1_knockouts)
    for numGene=1:length(r1_knockouts)

        count = count+1;
        load Ec_iJO1366.mat

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
            model.lb(rxnList(j)) = 0;
            model.ub(rxnList(j)) = 0;
        end

        [v,vg,exitFlag] = linprog(-model.c,[],[],model.S,zeros([1,1805]),model.lb,model.ub,options);

        if exitFlag==-2
            continue
        end

        if abs(vg)>.0001
            survivedKnockouts = cat(1,survivedKnockouts,[r1_knockouts(k,:),r1_knockouts(numGene)]);
            num_survivedKnockout = cat(1,num_survivedKnockout,[r1_num_knockouts(k,:),r1_num_knockouts(numGene)]);


        end
        if mod(count,100)==0
            progress = count/(length(r1_num_knockouts)*length(r1_num_knockouts));
            fprintf(string(progress))
        end

    end
end

fprintf("\n Complete!\n")

[~, idx] = unique(sort(num_survivedKnockout(:,1:length(survivedKnockouts(1,1))),2),'rows','stable');
num_survivedKnockout = num_survivedKnockout(idx,:);
survivedKnockouts = survivedKnockouts(idx,:);

r2_knockouts = survivedKnockouts;
r2_num_knockouts = num_survivedKnockout;

save 'r2_knockouts.mat' r2_knockouts
save 'r2_num_knockouts.mat' r2_num_knockouts


















