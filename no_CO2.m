function model = no_CO2(model)
% changes lower bound of net carbon dioxide transport from -1000 to 0, 
% essentially it goes from being able to import lots of CO2 to only 
% being able to export it

model.lb(85) = 0;

end