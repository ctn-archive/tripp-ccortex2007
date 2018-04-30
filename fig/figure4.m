% This script produces figure 4 of the article from saved data.  

popdata = load('data_populationSizeRate.mat');
popMeanErr = mean(popdata.meanE,3);
popSdErr = std(popdata.meanE,0,3);

popindices = [17 5 18 19 20 21];
jitterindices = 1:5; 
popsize = popdata.cases(popindices,1);
figure, hold on, plot(popsize, popMeanErr(popindices,jitterindices), 'k');
for i = popindices
    for j = jitterindices
        plot([popdata.cases(i,1) popdata.cases(i,1)], [popMeanErr(i,j)-popSdErr(i,j) popMeanErr(i,j)+popSdErr(i,j)], 'k')
    end
end
set(gcf, 'Position', [360 669 308 265])
xlabel('Population size (# neurons)')
ylabel('MSE')
for i = jitterindices
    h = text(40, popMeanErr(popindices(1),i), sprintf('%i ms', 1000*popdata.jitter(i)))
end

