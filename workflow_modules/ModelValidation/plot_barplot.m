%% Oliver Grimm and Lucy Hutchinson 2021
%
%% Make a bar plot of the condensed results
% Note this plot is customised for the results from the manuscript

function plot_barplot(res, pats_chosen,N, filename_save)
res = 100*res; % Get a percentage
mean_ctrl = mean(res(:,2:N),2);
std_ctrl = std(res(:,2:N),0,2);

figure('units','normalized','position',[0.1,0.1,0.7,0.5])
b2=bar(res,'LineWidth',1.5)
b2(1).FaceColor='Flat';
for ii = 2:N
    b2(ii).FaceColor='Flat';
    b2(ii).CData=[1,1,1];
end
b2(1).CData=repmat([0 0 0],length(pats_chosen),1);

for ii = 1:length(pats_chosen)
patlabels{ii}=['p',num2str(pats_chosen(ii))];
end
xticklabels(patlabels)
set(gca,'Fontsize',13)
ylabel('% of simulations accepted by SAM')
title('Accuracy of Model Predictions')
legend([b2(1),b2(2)],{'Population Parameters','Control Parameters'},'location','NorthEastOutside')

print(filename_save,'-dpng')