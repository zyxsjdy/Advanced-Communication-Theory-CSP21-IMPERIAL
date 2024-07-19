function [delayEst] = fChannelEstimation(symbolsOut, goldSeq, nPaths)
    % Function:
    %   - perform channel estimation for the desired source using the received
    %  signal
    %
    % InputArg(s):
    %   - symbolsOut (L*1): channel symbol chips received
    %   - goldSeq(W*1): gold sequence used in the modulation process
    %   - nPaths: number of paths for each source
    %
    % OutputArg(s):
    %   - delayEst: estimated delay of the signal paths
    %
    % Comments:
    %   - number of paths for each source should be known
    %
    % Author & Date: Yang (i@snowztail.com) - 21 Dec 18

    % obtain the maximum possible relative delay and number of signals
    [nDelays, nSignals] = size(goldSeq);
    % assume the max actual delay can be large (improve accuracy with cost of
    % complexity)
    delayMax = length(symbolsOut) - nDelays;
    % estimated delay for all paths
    delayEst = zeros(sum(nPaths), 1);
    % all possible correlation functions
    corFun = zeros(delayMax, nSignals);
    % initialise path counter
    pathCounter = 1;
    for iDelay = 1: delayMax
        % calculate the correlation functions for all possible delays
        corFun(iDelay, :) = abs(symbolsOut(iDelay: iDelay + nDelays - 1).' * goldSeq);
    end
    
    for iSignal = 1: nSignals
        % sort the correlation indexes
        [~, delayIndex] = sort(corFun(:, iSignal), 'descend');
        % find the residues of indexes and delete the repeated values
        uniqueResidue = unique(mod(delayIndex, nDelays), 'stable');
        % the first few residues can suggest the delays of paths of the
        % corresponding source
        delayEst(pathCounter: pathCounter + nPaths(iSignal) - 1) = sort(uniqueResidue(1: nPaths(iSignal)))-1;
        % update the path counter
        pathCounter = pathCounter + nPaths(iSignal);
    end
    % set the invalid estimation as zero
    delayEst(delayEst < 0) = length(goldSeq)-1;
end

