% this script produces figure 7 of the article 

load data_power.mat

figure, hold on 
plot(mistimeRate/lambda*100, trials(1,:), 'k')
plot(mistimeRate/lambda*100, trials(2,:), 'k--')
set(gca, 'XLim', [0 70])
set(gca, 'YLim', [0 1500])
set(gcf, 'Position', [360 669 308 265])
xlabel('% Mistimed Post-Synaptic Spikes')
ylabel('Trials')

clear all
load data_subtleDrive.mat

%The choice of balancing inhibition below (i.e. columns 5 and 15 from layer2a and 
%layer2b, respectively) is from subjective judgement of what level of balancing 
%inhibition produces the most stereotyped firing post-synaptically. Compare
%others by changing these indices as desired. 

%Note: we histogram two copies to convert from spikes/(100 trials) to spikes/s
foo = reshape(layer1(:,3,:), size(layer1,1)*size(layer1,3), 1); foo = foo - .1; foo = foo(find(foo > 0 & foo < .5)); figure, hist([foo; foo], 100)
hold on, plot([.05 .05 .1], [29 34 34 ], 'k'), set(gcf, 'Position', [405 711 506 160])
foo = reshape(layer2b(:,15,:), size(layer2b,1)*size(layer2b,3), 1); foo = foo - .1; foo = foo(find(foo > 0 & foo < .5)); figure, hist([foo; foo], 100)
hold on, plot([.05 .05], [29 34], 'k'), set(gcf, 'Position', [405 711 506 160])
foo = reshape(layer2a(:,5,:), size(layer2a,1)*size(layer2a,3), 1); foo = foo - .1; foo = foo(find(foo > 0 & foo < .5)); figure, hist([foo; foo], 100)
hold on, plot([.05 .05], [29 34], 'k'), set(gcf, 'Position', [405 711 506 160])
