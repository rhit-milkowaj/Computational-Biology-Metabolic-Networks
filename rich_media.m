function model = rich_media(model)
% changes all import lower bound to -1000, essentially the cell can import
% as much of whatever it wants
for i=9:332
    model.lb(i) = -1000;

end