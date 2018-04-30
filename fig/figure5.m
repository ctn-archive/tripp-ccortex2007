% this script produces figure 5 of the article 

popdata = load('data_populationSizeRate.mat');
popMeanErr = mean(popdata.meanE,3);
popSdErr = std(popdata.meanE,0,3);

rateindices = 1:8; 
figure(4), semilogx(popdata.cases(rateindices,2), popMeanErr(rateindices,jitterindices), 'k');
set(gca, 'NextPlot', 'add')
for i = rateindices
    for j = jitterindices
        plot([popdata.cases(i,2) popdata.cases(i,2)], [popMeanErr(i,j)-popSdErr(i,j) popMeanErr(i,j)+popSdErr(i,j)], 'k')
    end
end
set(gcf, 'Position', [360 669 308 265])
xlabel('Mean firing rate (Hz)')
ylabel('MSE')
for i = jitterindices
    h = text(popdata.cases(8,2)+50, popMeanErr(8,i), sprintf('%i ms', 1000*popdata.jitter(i)))
end

poissonrateindices = 9:16; 
semilogx(popdata.cases(poissonrateindices,2), popMeanErr(poissonrateindices,jitterindices), 'k:');
