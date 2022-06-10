%% Oliver Grimm and Lucy Hutchinson 2021
%
function plot_model_scenarios(pats_test, path, savedir)

figure('units','normalized','position',[0.1,0.1,0.85,0.85])
colours = [51,153,255; 255, 153, 51; 255,0,127]/255;

for p = 1:length(pats_test)
for scen = 1:3
path_curr = [path,num2str(scen),'_pat',num2str(pats_test(p))];
scenario_curr = load(path_curr);

    subplot(2,4,p)
    CD8_timecourse_mean_curr = nanmean(scenario_curr.CD8_timecourse_all,3);
    CD8_timecourse_mean_curr=squeeze(CD8_timecourse_mean_curr);
    % Stack the tiles and parsets on top of each other to reach a 2d matrix
    temp = CD8_timecourse_mean_curr(:,:,1);
    
    for ll = 2:size(CD8_timecourse_mean_curr,3)
        temp = vertcat(temp,CD8_timecourse_mean_curr(:,:,ll));
    end
    
    CD8_timecourse_all_curr  = temp;
    
    plot(nanmedian(CD8_timecourse_all_curr),'color',colours(scen,:),'linewidth',2)
    
    % Use quantiles and fill
    hold on
    quantiles_curr25 = quantile(CD8_timecourse_all_curr,[0.25,0.75]);
    
    % Put into a long vector for fill
    if sum(sum(isnan(quantiles_curr25)))
        quantiles25_fill=[quantiles_curr25(1,~isnan(quantiles_curr25(1,:))),fliplr(quantiles_curr25(2,~isnan(quantiles_curr25(2,:))))];
    else
        quantiles25_fill = [quantiles_curr25(1,:),fliplr(quantiles_curr25(2,:))];
    end
    
    fill([1:length(quantiles25_fill)/2,length(quantiles25_fill)/2:-1:1]',quantiles25_fill',colours(scen,:),'facealpha',0.2,'linestyle','none')
    
    xlabel('time (days)')
    ylabel('CD8 number per tile')
    xlim([0,50])
    ylim([0 2000])
    title(['Patient ',num2str(pats_test(p))])
end
end
mkdir(savedir)
print([savedir,'/TimecourseCD8_median_iqr_scenarios_allpats'],'-dpng')
