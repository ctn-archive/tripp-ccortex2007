% This script produces the figure 2 panels that rely on data produced by
% experiment scripts (panels F and G). 

% plot error vs. coefficient of variation with spike jitter 
covdata = load('data_COV.mat');
meanErr = mean(covdata.meanErr,3);
sdErr = std(covdata.meanErr,0,3);

figure(1)
loglog(covdata.meanCOV, meanErr', 'k')
set(gcf, 'Position', [360 669 308 265])
set(gca, 'NextPlot', 'add')
for i = 1:size(meanErr, 1)
    for j = [1 6] 
        loglog([covdata.meanCOV(i) covdata.meanCOV(i)], [meanErr(i,j)-sdErr(i,j) meanErr(i,j)+sdErr(i,j)], 'k')
    end
end
for i = 1:size(covdata.noise,1)
    h = text(covdata.meanCOV(3), meanErr(3,i), sprintf('%i ms', 1000*covdata.noise(i,1)))
end
% izdata = load('data_Izhikevich.mat')
% izErr = mean(izdata.meanErr,3);
% plot(izdata.meanCOV(1), izErr(1,:), 'ko');
% plot(izdata.meanCOV(2), izErr(2,:), 'kx');
% xlabel('Mean COV')
% ylabel('MSE')
% set(gca, 'XLim', [.1 2])
% 
% % plot error vs. coefficient of variation with noise spikes 
% j0 = 1:3; % indices of noise conditions with jitter = 0ms
% j4 = 4:6; % indices of noise conditions with jitter = 4ms
% 
% covdata = load('data_COVNS.mat');
% meanErr = mean(covdata.meanErr,3);
% sdErr = std(covdata.meanErr,0,3);
% 
% figure(2)
% loglog(covdata.meanCOV, meanErr(:,j0)', 'k')
% set(gca, 'NextPlot', 'add')
% loglog(covdata.meanCOV, meanErr(:,j4)', 'k:')
% set(gcf, 'Position', [360 669 308 265])
% 
% for i = 1:size(meanErr, 1)
%     for j = [j0(1) j0(end)] 
%         loglog([covdata.meanCOV(i) covdata.meanCOV(i)], [meanErr(i,j)-sdErr(i,j) meanErr(i,j)+sdErr(i,j)], 'k')
%     end
% end
% for i = j0
%     h = text(covdata.meanCOV(3), meanErr(3,i), sprintf('%i %%', 100*covdata.noise(i,2)))
% end
% xlabel('Mean COV')
% ylabel('MSE')
% set(gca, 'XLim', [.1 2])
