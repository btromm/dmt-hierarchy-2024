%% gec analysis
load('/Users/myco/computer_magic/nifty_utils/dbs2yeo/dk80toyeo7.mat'); %32-50 are subcortical structures, not part of RSNs!
%load your data

NSUB=16;

%G_totalCT=G_totalOccasionalSober;
%G_totalSZ=G_totalOccasionalCannabis;
G_totalCT_Yeo=zeros(NSUB,7);
G_totalSZ_Yeo=zeros(NSUB,7);

%G_totalCT_Yeo=decompose_rsn(G_totalCT);
%G_totalSZ_Yeo=decompose_rsn(G_totalSZ);


% total amount of orchestration by network
for nsub=1:NSUB
    for i=1:80
        if av(i)~=0
            G_totalCT_Yeo(nsub,av(i))=G_totalCT_Yeo(nsub,av(i))+G_totalCT(nsub,i);
            G_totalSZ_Yeo(nsub,av(i))=G_totalSZ_Yeo(nsub,av(i))+G_totalSZ(nsub,i);
        end
    end
end



%% G_total

for n=1:7 % for each RSN
    tmpCT=G_totalCT_Yeo(:,n)';
    tmpSZ=G_totalSZ_Yeo(:,n)';
    stats=permutation_htest2_np([tmpCT tmpSZ],[ones(1,numel(tmpCT)) 2*ones(1,numel(tmpSZ))],5000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
    t_idv(n)=stats.tvals;
    diffs_idv(n)=stats.diffs;
end
nsig=FDR_benjHoch(pp,0.05);

figure;
Y=[nanmean(G_totalCT_Yeo)' nanmean(G_totalSZ_Yeo)']; 
Error=[(std(G_totalCT_Yeo)./sqrt(length(G_totalCT_Yeo)))' (std(G_totalSZ_Yeo)./sqrt(length(G_totalSZ_Yeo)))']; %standard error
b=bar(Y,'grouped')
hold on;
% Error bars
ngroups = size(Y, 1);
nbars = size(Y, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Loop through each group and add error bars
for i = 1:nbars
    % Calculate the center positions of the bars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    % Create error bars
    errorbar(x, Y(:,i), Error(:,i), 'k', 'linestyle', 'none');
end

ylim([0 30])
%ylabel('In-degree')
legend({'Baseline','Ayahuasca'},'Location','Northwest')
xticklabels({'VIS','SM','DAN','SAL','LIM','FP','DMN'})
map=cbrewer2('Pastel1',2);
b(1).FaceColor=map(2,:);
b(2).FaceColor=map(1,:);
set(gca,'fontname','SF Pro');
set(gca,'fontsize',18);

%% same thing with hl's
%load your data

hl_totalCT_Yeo=decompose_rsn(hierarchicallevelsOccasionalSober);
hl_totalSZ_Yeo=decompose_rsn(hierarchicallevelsOccasionalCannabis);

for n=1:7 % for each RSN
    tmpCT=hl_totalCT_Yeo(:,n)';
    tmpSZ=hl_totalSZ_Yeo(:,n)';
    stats=permutation_htest2_np([tmpCT tmpSZ],[ones(1,numel(tmpCT)) 2*ones(1,numel(tmpSZ))],5000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
    t_idv(n)=stats.tvals;
    diffs_idv(n)=stats.diffs;
end
nsig=FDR_benjHoch(pp,0.05);

figure;
Y=[nanmean(hl_totalCT_Yeo)' nanmean(hl_totalSZ_Yeo)']
Error=[(std(hl_totalCT_Yeo)./sqrt(length(hl_totalCT_Yeo)))' (std(hl_totalSZ_Yeo)./sqrt(length(hl_totalSZ_Yeo)))']
b=bar(Y,'grouped')
hold on;
% Error bars
ngroups = size(Y, 1);
nbars = size(Y, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));

% Loop through each group and add error bars
for i = 1:nbars
    % Calculate the center positions of the bars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    % Create error bars
    errorbar(x, Y(:,i), Error(:,i), 'k', 'linestyle', 'none');
end

ylim([0 0.75])
%ylabel('In-degree')
legend({'Baseline','Cannabis'},'Location','Northwest')
xticklabels({'VIS','SM','DAN','SAL','LIM','FP','DMN'})
map=cbrewer2('Pastel1',2);
b(1).FaceColor=map(2,:);
b(2).FaceColor=map(1,:);
set(gca,'fontname','SF Pro');
set(gca,'fontsize',18);