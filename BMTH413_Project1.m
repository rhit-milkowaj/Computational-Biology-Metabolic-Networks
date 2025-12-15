%% BMTH413 Project 1 - Metabolic Models

% Written by Taylor Donen

%% Run Experiment: Glucose 
% Experiment 1
% Requires no change 

% Import data
model = importdata("Ec_iJO1366.mat");
data = importdata("Genes_for_BMTH413.xlsx");

comp_essential = [];
growth_rates = [];

% Computational Results
for numGene = 1:length(model.genes)

   geneVector = numGene;
   rxnList = [];

   for i=1:length(geneVector)
      rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
   end

   % Set lower and upper bounds of rxnList to zero 
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

   % Change upper and lower bounds
   new_ub = model.ub;
   new_lb = model.lb;

   for i=1:length(new_ub)
        if ismember(i, rxnList)
            new_ub(i) = 0;
            new_lb(i) = -10;
        end
   end

   % Solve linear program with program ub and lb
   options = optimoptions('linprog','Display','none');
   [v,g_star] = linprog(-model.c,[],[],model.S,model.b,new_lb,new_ub,options);

   % Check if g_star is < 0, then the gene is essential
   if abs(g_star) < (0.004)
       comp_essential = [comp_essential, model.genes(numGene)];
       growth_rates = [growth_rates, g_star];
   end

end

comp_nonessential = {};
for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ~ismember(gene_name, comp_essential)
        comp_nonessential = [comp_nonessential; gene_name];
    end
end

% Experimental Results
% For Glucose Medium

exp_always_essential = data(4:141,6);

exp_essential = data(4:23,10);
exp_essential = [exp_essential; exp_always_essential];

exp_nonessential = {};
for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ~ismember(gene_name, exp_essential)
        exp_nonessential = [exp_nonessential; gene_name];
    end
end

% Compare Experimental vs. Computational Essential Genes
tp_glucose = {};
fp_glucose = {};
fn_glucose = {};
tn_glucose = {};

for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ismember(gene_name, exp_essential) && ismember(gene_name, comp_essential)
        tp_glucose = [tp_glucose; gene_name];
    elseif ismember(gene_name, exp_nonessential) && ismember(gene_name, comp_essential)
        fp_glucose = [fp_glucose; gene_name];
    elseif ismember(gene_name, exp_essential) && ismember(gene_name, comp_nonessential)
        fn_glucose = [fn_glucose; gene_name];
    elseif ismember(gene_name, exp_nonessential) && ismember(gene_name, comp_nonessential)
        tn_glucose = [tn_glucose; gene_name];
    end
end

tp_number = numel(tp_glucose);
fp_number = numel(fp_glucose);
fn_number = numel(fn_glucose);
tn_number = numel(tn_glucose);


fprintf(['There are %d genes that are true positive (experimentally & ' ...
    'computationally essential) for glucose.\n\n'], tp_number);
fprintf(['There are %d genes that are false positive (experimentally non-essential & ' ...
    'computationally essential) for glucose.\n\n'], fp_number);
fprintf(['There are %d genes that are false negative (experimentally essential & ' ...
    'computationally non-essential) for glucose.\n\n'], fn_number);
fprintf(['There are %d genes that are true negative (experimentally & ' ...
    'computationally non-essential) for glucose.\n\n'], tn_number);


%% Run Experiment: Glycerol 
% Experiment 2
% Glucose lb = 0, glycerol lb = -10

% Import data
model = importdata("Ec_iJO1366.mat");
data = importdata("Genes_for_BMTH413.xlsx");

comp_essential = [];
growth_rates = [];

% Computational Results
for numGene = 1:length(model.genes)

   geneVector = numGene;
   rxnList = [];

   for i=1:length(geneVector)
      rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
   end

   % Set lower and upper bounds of rxnList to zero 
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

   % Change upper and lower bounds
   new_ub = model.ub;
   new_lb = model.lb;

   % Glucose reaction index = 164, knock out glucose reaction
   new_lb(164) = 0;   

   % Glycerol reaction index = 171, allow glycerol uptake
   new_lb(171) = -10;   % allow glycerol uptake 

    % Set environment to not have O2 (lb to zero)
    % 02 reaction index = 252
   % new_ub(252) = 0;
   % new_lb(252) = 0;

    % % Set environment to not have CO2 (lb  to zero)
    % % C02 reaction index = 85
    % new_ub(85) = 0;
    % new_lb(85) = 0;

    % % Set environment to not have H20 (lb to zero)
    % % H20 reaction index = 187
    % new_ub(187) = 0;
    % new_lb(187) = 0;

    % % Set environment to not have H2 (lb to zero)
    % % H2 reaction index = 186
    % new_ub(186) = 0;
    % new_lb(186) = 0;

    % % Set environment to not have K (lb to zero)
    % % K reaction index = 206
    % new_ub(206) = 0;
    % new_lb(206) = 0;

   for i=1:length(new_ub)
        if ismember(i, rxnList)
            new_ub(i) = 0;
            new_lb(i) = 0;
        end
   end

   % Solve linear program with program ub and lb
   options = optimoptions('linprog','Display','none');
   [v,g_star] = linprog(-model.c,[],[],model.S,model.b,new_lb,new_ub,options);

   % Check if g_star is < 0, then the gene is essential
   if abs(g_star) < (0.004)
       comp_essential = [comp_essential, model.genes(numGene)];
       growth_rates = [growth_rates, g_star];
   end

end

comp_nonessential = {};
for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ~ismember(gene_name, comp_essential)
        comp_nonessential = [comp_nonessential; gene_name];
    end
end

% Experimental Results
% For Glycerol Medium

exp_always_essential = data(4:141,6);

exp_essential = data(4:23,12);
exp_essential = [exp_essential; exp_always_essential];

exp_nonessential = {};
for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ~ismember(gene_name, exp_essential)
        exp_nonessential = [exp_nonessential; gene_name];
    end
end

% Compare Experimental vs. Computational Essential Genes

tp_glycerol = {};
fp_glycerol = {};
fn_glycerol = {};
tn_glycerol = {};

for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ismember(gene_name, exp_essential) && ismember(gene_name, comp_essential)
        tp_glycerol = [tp_glycerol; gene_name];
    elseif ismember(gene_name, exp_nonessential) && ismember(gene_name, comp_essential)
        fp_glycerol = [fp_glycerol; gene_name];
    elseif ismember(gene_name, exp_essential) && ismember(gene_name, comp_nonessential)
        fn_glycerol = [fn_glycerol; gene_name];
    elseif ismember(gene_name, exp_nonessential) && ismember(gene_name, comp_nonessential)
        tn_glycerol = [tn_glycerol; gene_name];
    end
end

tp_number = numel(tp_glycerol);
fp_number = numel(fp_glycerol);
fn_number = numel(fn_glycerol);
tn_number = numel(tn_glycerol);


fprintf(['There are %d genes that are true positive (experimentally & ' ...
    'computationally essential) for glycerol.\n\n'], tp_number);
fprintf(['There are %d genes that are false positive (experimentally non-essential & ' ...
    'computationally essential) for glycerol.\n\n'], fp_number);
fprintf(['There are %d genes that are false negative (experimentally essential & ' ...
    'computationally non-essential) for glycerol.\n\n'], fn_number);
fprintf(['There are %d genes that are true negative (experimentally & ' ...
    'computationally non-essential) for glycerol.\n\n'], tn_number);

%% Run Experiment: Rich Media 
% Experiment 3
% Change all import lower bounds to -1000

% Import data
model = importdata("Ec_iJO1366.mat");
data = importdata("Genes_for_BMTH413.xlsx");

comp_essential = [];
growth_rates = [];

% Computational Results
for numGene = 1:length(model.genes)

   geneVector = numGene;
   rxnList = [];

   for i=1:length(geneVector)
      rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
   end

   % Set lower and upper bounds of rxnList to zero 
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

   % Change upper and lower bounds
   new_ub = model.ub;
   new_lb = model.lb;

   
   % Changes all import lower bound to -1000 for rich media
   for i=9:332
       new_lb(i) = -1000;
   end

   for i=1:length(new_ub)
        if ismember(i, rxnList)
            new_ub(i) = 0;
            new_lb(i) = 0;
        end
   end

   % Solve linear program with program ub and lb
   options = optimoptions('linprog','Display','none');
   [v,g_star] = linprog(-model.c,[],[],model.S,model.b,new_lb,new_ub,options);

   % Check if g_star is < 0, then the gene is essential
   if abs(g_star) < (0.004)
       comp_essential = [comp_essential, model.genes(numGene)];
       growth_rates = [growth_rates, g_star];
   end

end

comp_nonessential = {};
for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ~ismember(gene_name, comp_essential)
        comp_nonessential = [comp_nonessential; gene_name];
    end
end

% Experimental Results
% For Rich Medium

exp_essential = data(4:141,6);

exp_nonessential = {};
for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ~ismember(gene_name, exp_essential)
        exp_nonessential = [exp_nonessential; gene_name];
    end
end

% Compare Experimental vs. Computational Essential Genes

tp_rich = {};
fp_rich = {};
fn_rich = {};
tn_rich = {};

for i=1:length(model.genes)
    gene_name = model.genes{i};
    if ismember(gene_name, exp_essential) && ismember(gene_name, comp_essential)
        tp_rich = [tp_rich; gene_name];
    elseif ismember(gene_name, exp_nonessential) && ismember(gene_name, comp_essential)
        fp_rich = [fp_rich; gene_name];
    elseif ismember(gene_name, exp_essential) && ismember(gene_name, comp_nonessential)
        fn_rich = [fn_rich; gene_name];
    elseif ismember(gene_name, exp_nonessential) && ismember(gene_name, comp_nonessential)
        tn_rich = [tn_rich; gene_name];
    end
end

tp_number = numel(tp_rich);
fp_number = numel(fp_rich);
fn_number = numel(fn_rich);
tn_number = numel(tn_rich);


fprintf(['There are %d genes that are true positive (experimentally & ' ...
    'computationally essential) in rich media.\n\n'], tp_number);
fprintf(['There are %d genes that are false positive (experimentally non-essential & ' ...
    'computationally essential) in rich media.\n\n'], fp_number);
fprintf(['There are %d genes that are false negative (experimentally essential & ' ...
    'computationally non-essential) in rich media.\n\n'], fn_number);
fprintf(['There are %d genes that are true negative (experimentally & ' ...
    'computationally non-essential) in rich media.\n\n'], tn_number);

%% Test O2, CO2, H2O, H2, K on Glucose Media

clear all; close all; clc;

% Import data
model = importdata("Ec_iJO1366.mat");
data = importdata("Genes_for_BMTH413.xlsx");

essential_genes = [];
growth_rates = [];

% Computational Results
for numGene = 1:length(model.genes)

   geneVector = numGene;
   rxnList = [];

   for i=1:length(geneVector)
      rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
   end

   % Set lower and upper bounds of rxnList to zero 
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

   % Change upper and lower bounds
   new_ub = model.ub;
   new_lb = model.lb;

   % Set environment to not have O2 (lb to zero)
   % 02 reaction index = 252
   % new_lb(252) = 0;

    % % Set environment to not have CO2 (lb  to zero)
    % % C02 reaction index = 85
    % new_lb(85) = 0;

    % Set environment to not have H20 (lb to zero)
    % H20 reaction index = 187
    % new_lb(187) = 0;

    % % Set environment to not have H2 (lb to zero)
    % % H2 reaction index = 186
    % new_lb(186) = 0;

    % % Set environment to not have K (lb to zero)
    % K reaction index = 206
    % new_lb(206) = 0;

   for i=1:length(new_ub)
        if ismember(i, rxnList)
            new_ub(i) = 0;
            new_lb(i) = 0;
        end
   end

   % Solve linear program with program ub and lb
   options = optimoptions('linprog','Display','none');
   [v,g_star] = linprog(-model.c,[],[],model.S,model.b,new_lb,new_ub,options);

   % Check if g_star is < 0, then the gene is essential
   if abs(g_star) < (0.004)
       essential_genes = [essential_genes, model.genes(numGene)];
       growth_rates = [growth_rates, g_star];
   end

end


%% Normal Distribution 

model = importdata("Ec_iJO1366.mat");
alpha = 0.1;

% 8th column: reaction that measures growth
mu = full(model.S(:,8));
muInd = find(mu~=0);
mu = mu(muInd);

std = (0.1*mu(:))/3;

X = zeros(100,72);
y = zeros(100,1);

% Iterate through 100 times 
for j=1:100

    % randn changes each interation
    diff_mu = mu + std(:).*randn(length(mu),1);
    X(j,:) = diff_mu;

    Stmp = model.S;
    Stmp(muInd,8) = diff_mu;

    comp_essential = [];
    all_growth_rates = zeros(1367,1);

    % Find new growth rate 
    for numGene = 1:length(model.genes)

        geneVector = numGene;
        rxnList = [];

        for i=1:length(geneVector)
            rxnList = union(rxnList,find(model.rxnGeneMat(:,geneVector(i))==1));
        end

        % Set lower and upper bounds of rxnList to zero 
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

        % Change upper and lower bounds
        new_ub = model.ub;
        new_lb = model.lb;

        for i=1:length(new_ub)
            if ismember(i, rxnList)
                new_ub(i) = 0;
                new_lb(i) = 0;
            end
        end

        % Solve linear program with program ub and lb
        options = optimoptions('linprog','Display','none');
        [v,g_star] = linprog(-model.c,[],[],Stmp,model.b,new_lb,new_ub,options);

        if ~isempty(g_star)
            all_growth_rates(numGene) = g_star;
        end
           
        % Check if g_star is < 0, then the gene is essential
        if abs(g_star) < (0.004)
            comp_essential = [comp_essential, model.genes(numGene)];
        end

    end

    y(j) = mean(all_growth_rates);
end


%% Use linear regression code 

[a,Rsqr,pVal,yHat,cInt,pInt] = Donen_LinearRegression_V1(X,y);

plot(y, yHat, 'k.', 'MarkerSize', 15, 'Color', 'Blue')
xlabel('True Mean Growth')
ylabel('Predicted Mean Growth')
