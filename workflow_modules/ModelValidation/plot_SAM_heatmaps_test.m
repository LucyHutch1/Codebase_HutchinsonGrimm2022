%% Oliver Grimm and Lucy Hutchinson 2021
%
function [] = plot_SAM_heatmaps_test(pats_chosen, accepted1cont,v2_accept_matrix, pars_all, par_names, total_sets,filename)


n_pars = size(pars_all,1);
n_pats = length(pats_chosen);

% Scale the parameters between 1 and 0
pars_all_scaled = (pars_all-min(pars_all))./(max((pars_all-min(pars_all))));

% Reshape the matrices to have the test columns first, followed by all
% control columns
accepted1cont_all = reshape(accepted1cont,n_pars,n_pats*total_sets);
v2accept_all = reshape(v2_accept_matrix,n_pars,n_pats*total_sets);
v2accept_all_scaled = v2accept_all;
v2accept_all_scaled(v2accept_all>1)=1./(v2accept_all(v2accept_all>1)); % Scale so everything is below 1

% Add nans where there are nans in the accept1contmatrix
v2accept_all_scaled_wFlags = v2accept_all_scaled;
v2accept_all_scaled_wFlags(isnan(accepted1cont_all))=nan;

figure('units','normalized','position',[0 0 1 1],'inverthardcopy','off')
subplot(1,8,1)
imagesc(pars_all_scaled)
ax=gca;ax.YGrid='on';ax.YMinorGrid='on'; ax.GridColor = 'k'; ax.GridAlpha=1;
ax.FontSize=12;
% Move to the left
pos = get(ax,'Position');
ax.Position = [pos(1)-0.05,pos(2),pos(3),pos(4)];

colormap(ax,'autumn')
xticks([1,2,3,4]);
xticklabels(par_names)
xtickangle(45);
ylabel('Population parameter set ID','fontsize',14)
c=colorbar;
c.FontSize=14;
c.Label.String='scaled paramter value';
c.Label.Rotation = 270;
c.Label.FontSize = 14;
pos = get(c,'Position');
c.Label.Position = [pos(1)+7 pos(2)+0.4]; % to change its position

title('Final par sets','fontsize',14)

subplot(1,9,[2:5])
h=imagesc(accepted1cont_all);
ax1=gca;ax1.YGrid='on';ax1.YMinorGrid='on'; ax1.GridColor = 'w'; ax1.GridAlpha=1;
ax1.FontSize=12;
colormap(ax1,'summer')
xticks(1:size(accepted1cont_all,2))
xticklabels(repmat(pats_chosen',1,7));
xtickangle(90);
xlabel('Patient','fontsize',14)
ylabel('Population parameter set ID','fontsize',14)

% Set the nans to transparent then black
set(h,'alphadata',~isnan(accepted1cont_all))
set(gca,'color','black')

hold on
% Add vertical lines and subdivisions
for ii = 1:total_sets
plot([n_pats,n_pats]+[0.5,0.5]+(ii-1)*n_pats,[0,62.5],'w-','linewidth',2)
if ii == 1
str_temp = 'Test';
else
    str_temp = ['Control ',num2str(ii-1)];
end
text(n_pats/2+(ii-1)*n_pats-2,0.3,str_temp,'fontsize',14)
end
annotation('textbox',[0.3,0.955,0.3,0.03],'string','SAM: Test and Control results','fitboxtotext','on','fontsize',16,'LineStyle','none','FontWeight','bold')


subplot(1,9,[6:9])
h1=imagesc(v2accept_all_scaled_wFlags);
ax2=gca;ax2.YGrid='on';ax2.YMinorGrid='on'; ax2.GridColor = 'w'; ax2.GridAlpha=1;
ax2.FontSize=12;
% Move to the right
pos = get(ax2,'Position');
ax2.Position = [pos(1)+0.05,pos(2),pos(3),pos(4)];

colormap(ax2,'summer')
xticks(1:size(accepted1cont_all,2))
xticklabels(repmat(pats_chosen',1,7));
xtickangle(90);
xlabel('Patient','fontsize',14)
c1=colorbar;
c1.FontSize=14;
c1.Label.String='SAM/VarSAM value';
c1.Label.Rotation = 270;
c1.Label.FontSize = 14;
pos = get(c1,'Position');
c1.Label.Position = [pos(1)+3.5 pos(2)+0.4]; % to change its position
% Set the nans to transparent then black
set(h1,'alphadata',~isnan(v2accept_all_scaled_wFlags))
set(gcf,'color','white')
set(gca,'color','black')


ylabel('Population parameter set ID','fontsize',14)

hold on
% Add vertical lines and subdivisions
for ii = 1:total_sets
plot([n_pats,n_pats]+[0.5,0.5]+(ii-1)*n_pats,[0,62.5],'w-','linewidth',2)
if ii == 1
str_temp = 'Test';
else
    
    str_temp = ['Control ',num2str(ii-1)];
end
text(n_pats/2+(ii-1)*n_pats-2,0.3,str_temp,'fontsize',14)
end
annotation('textbox',[0.7,0.955,0.3,0.03],'string','VarSAM: Test and Control results','fitboxtotext','on','fontsize',16,'LineStyle','none','FontWeight','bold')


print(filename,'-dpng')
