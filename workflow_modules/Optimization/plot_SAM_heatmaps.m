%% Oliver Grimm and Lucy Hutchinson 2021
%
function [] = plot_SAM_heatmaps(pats_chosen, SAM_matrix_NANIndicator, VarSAM_scale, pars_all, par_names, filename)

% Scale the parameters between 1 and 0
transform_pars_all = [log10(pars_all(:,1)),log10(pars_all(:,2)),log10(pars_all(:,3)),pars_all(:,4)];
pars_all_scaled = (transform_pars_all-min(transform_pars_all))./range(transform_pars_all);

% Sort by mean number points within bounds across all pats (rows)
SAM_matrix_NANIndicator_mean = nanmean(SAM_matrix_NANIndicator,2);
SAM_matrix_NANIndicator_mean_minus1 = SAM_matrix_NANIndicator_mean;
SAM_matrix_NANIndicator_mean_minus1(isnan(SAM_matrix_NANIndicator_mean))=-1;
[~,ind_SAM_matrix_NANIndicator_mean] = sort(SAM_matrix_NANIndicator_mean_minus1);
pars_sort = pars_all_scaled(ind_SAM_matrix_NANIndicator_mean,:);
SAM_matrix_NANIndicator_sort = SAM_matrix_NANIndicator(ind_SAM_matrix_NANIndicator_mean,:);

% Sort VarSAM condition by same order as SAM matrix
VarSAM_sort = VarSAM_scale(ind_SAM_matrix_NANIndicator_mean,:);

% Add NaNs to the VarSAM matrix to match the SAM matrix with NaN indicator
VarSAM_sort(isnan(SAM_matrix_NANIndicator_sort))=nan;

figure('units','normalized','position',[0 0 1 1],'inverthardcopy','off')
subplot(1,5,1)
imagesc(pars_sort)
ax=gca;ax.YGrid='on';ax.YMinorGrid='on'; ax.GridColor = 'k'; ax.GridAlpha=1;
ax.FontSize=12;
% Move to the left
pos = get(ax,'Position');
ax.Position = [pos(1)-0.05,pos(2),pos(3),pos(4)];
xticks([1,2,3,4]);
xticklabels(par_names)
xtickangle(45);
colormap(ax,'autumn')
c1=colorbar;
c1.FontSize=12;
c1.Label.String='scaled parameter value';
c1.Label.Rotation = 270;
c1.Label.FontSize = 14;
pos = get(c1,'Position');
c1.Label.Position = [pos(1)+5 pos(2)+0.4]; % to change its position
ylabel('Parameter set ID (sorted)','Fontsize',14)

title('Par sets sorted by SAM','fontsize',16)

subplot(1,5,[2,3])
h=imagesc(SAM_matrix_NANIndicator_sort);
ax1=gca;ax1.YGrid='on';ax1.YMinorGrid='on'; ax1.GridColor = 'w'; ax1.GridAlpha=1;
ax1.FontSize=9;
xticks(1:length(pats_chosen))
xticklabels(pats_chosen);
xtickangle(90);
colormap(ax1,'summer')
xlabel('Patient','Fontsize',14)
title('SAM: Optimization parameters','fontsize',16)
ylabel('Parameter set ID (sorted)','Fontsize',14)

% Set the nans to transparent then black
set(h,'alphadata',~isnan(SAM_matrix_NANIndicator_sort))
set(gca,'color','black')


subplot(1,5,[4,5])
h1=imagesc(VarSAM_sort);
ax2=gca;ax2.YGrid='on';ax2.YMinorGrid='on'; ax2.GridColor = 'w'; ax2.GridAlpha=1;
ax2.FontSize=9;
% Move to the right
pos = get(ax2,'Position');
ax2.Position = [pos(1)+0.05,pos(2),pos(3),pos(4)];

xticks(1:length(pats_chosen))
xticklabels(pats_chosen);
xtickangle(90);
colormap(ax2,'summer')
xlabel('Patient','Fontsize',14)
ylabel('Parameter set ID (sorted)','Fontsize',14)
title('VarSAM: Optimization parameters','fontsize',16)
c=colorbar;
c.FontSize=12;
c.Label.String='SAM/VarSAM value';
c.Label.Rotation = 270;
pos = get(c,'Position');
c.Label.Position = [pos(1)+3 pos(2)+0.4]; % to change its position
c1.Label.FontSize = 14;

% Set the nans to transparent then black
set(h1,'alphadata',~isnan(VarSAM_sort))
set(gcf,'color','w')
set(gca,'color','black')%

% Save the plots
 print(filename,'-dpng')
  print(filename,'-dtiff')

