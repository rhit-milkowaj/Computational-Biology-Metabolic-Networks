function progressBar = getProgressBar(progress, barLength)
    % getProgressBar Creates a text-based progress bar
    filledLength = round(barLength * progress);
    emptyLength = barLength - filledLength;
    progressBar = ['[' repmat('#', 1, filledLength) repmat('-', 1, emptyLength) ']'];
end