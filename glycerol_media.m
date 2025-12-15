function model = glycerol_media(model)
% Changes the media from a glucose media to a glycerol media

model.lb(164) = 0; % sets glucose transport into cell to 0
model.lb(174) = -10; % sets glycerol transport into cell max to 10


end