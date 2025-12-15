function model = no_oxygen(model)
% changes lower bound of net oxygen transport from -1000 to 0, essentially
% it goes from being able to import lots of oxygen to only being able to
% export it

model.lb(252) = 0;

end