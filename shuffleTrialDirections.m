%% function to shuffle directions
% inputs:
% directions: 1d vector with directions
% numTrialsPerBlock: number of trial per phase (e.g.: 8 directions)
% numBlocks: number of blocks (or repeats; e.g.: 3)
%
% Make sure that length(directions) = numTrialsPerBlock * numBlocks


function [shuffledDirs] = shuffleTrialDirections(directions,numTrialsPerBlock, numBlocks)

shuffledInd = zeros(numTrialsPerBlock*numBlocks,1);

for i = 1:numBlocks

    shuffledInd(i*numTrialsPerBlock-numTrialsPerBlock+1:i*numTrialsPerBlock) = randperm(numTrialsPerBlock);
    
end

shuffledDirs = directions(shuffledInd);