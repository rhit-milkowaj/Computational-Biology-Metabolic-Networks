function updateProgressDisplay(i, totalIterations, startTime, barLength)
% updateProgressDisplay Updates the command window with progress information.
%   Inputs:
%     i              - Current iteration number.
%     totalIterations - Total number of iterations expected.
%     startTime      - The time returned by 'tic' at the start of the process.
%     barLength      - (Optional) Length of the text-based progress bar. Default is 30.

if nargin < 4
    barLength = 30; % Default bar length
end

elapsedTime = toc(startTime); % get time elapsed since start time (seconds)
progress = i / totalIterations; % fraction of total progress made

estimatedTotalTime = elapsedTime / progress;
remainingTime = estimatedTotalTime - elapsedTime;

% Convert remaining time (in seconds) to hours, minutes, seconds 
[h, m, s] = sec2hms(remainingTime);

% gets current time to estimate ending time
t_now = datetime('now','TimeZone','local','Format','HH:mm:ss');
t_est = t_now + hours(h) + minutes(m) + seconds(s);
t_est_str = datestr(t_est, 'HH:MM:SS');
print_str = append('Est. Finish Time: ',  t_est_str,'\r');

% Display update 
fprintf('Processing: %s | Progress: %.1f%% | ', ...
    getProgressBar(progress, barLength), progress * 100);
fprintf(print_str)

end



