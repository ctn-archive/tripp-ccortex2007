% function spikeTrain = space(spikeTrain, RT) spaces out spikes in a single spike train 
% according to an absolute refractory time. 
% 
% The options are to drop or shift them.  We pretend that this represents
% a coincidence of current from correlated and uncorrelated sources, and
% drop them.  If you change this, be sure to shift the time of the
% offending spike rather than just increasing the ISI, as the latter
% would shift subsequent correlated activity.    
function spikeTrain = space(spikeTrain, RT)
    done = 0;
    while (~done) %have to go one at a time, because removing one spike effects previous and following ISI
        spikeTrain = sort(spikeTrain);
        ISI = diff(spikeTrain);
        closeIndices = find(diff(spikeTrain) < RT & spikeTrain(2:end) > 0) + 1;
        if (length(closeIndices) == 0)
            done = 1;
        else 
            spikeTrain(closeIndices(1)) = 0;
        end
    end
    spikeTrain = sort(spikeTrain);