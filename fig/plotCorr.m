% function plotCorr(datafile, jitterIndex) plots a summary of MSE vs
% correlation. 
% 
% dataFile: filename of saved data from a exp_correlation.m
% jitterIndex: index of jitter level in this file 

function plotCorr(datafile, jitterIndex)
    data = load(datafile); 
    
    hold on
    plotData(data, data.epochRateIndices, 'k-', jitterIndex);
    plotData(data, data.thresholdIndices10, 'k:', jitterIndex);
    plotData(data, data.thresholdIndices22, 'k:', jitterIndex);
    plotData(data, data.thresholdIndices55, 'k:', jitterIndex);
    xlabel('Pairwise correlation')
    ylabel('MSE')
    
    addLabel(data, data.thresholdIndices10(1), jitterIndex, 'a');
    addLabel(data, data.thresholdIndices22(1), jitterIndex, 'b');
    addLabel(data, data.thresholdIndices55(1), jitterIndex, 'g');

function addLabel(data, caseIndex, jitterIndex, labelText)
    h = text(data.meanPeakCorrelation(caseIndex)+.005, data.meanE(caseIndex,jitterIndex), labelText)
    set(h, 'FontName', 'Symbol')  

function plotData(data, caseIndices, symbol, jitterIndex)
    x = data.meanPeakCorrelation(caseIndices);
    y = squeeze(data.meanE(caseIndices,jitterIndex,:));
    meanY = mean(y, 2);
    sdY = std(y, 0, 2);
 	plot(x, meanY, symbol);
 
	hold on    
	for i = 1:length(x)
        plot([x(i) x(i)], [meanY(i)-sdY(i) meanY(i)+sdY(i)], symbol) 
	end
