function [h, m, s] = sec2hms(t)
    % converts seconds to hours, minutes, and seconds
    h = floor(t / 3600);
    t = t - h * 3600;
    m = floor(t / 60);
    s = t - m * 60;
end