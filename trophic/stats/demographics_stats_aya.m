%% AYAHUASCA STATS
%SZ=AYA, CT=PCB

%% Loading stuff
clear all;
path2data='/Users/myco/human_brains/ayahuasca';
cd(path2data);
load results_Ceff_aya;
load results_aya_modelfree.mat
load demographics.mat;
numSub=24;
results_path='/Users/myco/human_brains/INSIDEOUT_THESIS';
load '/Users/myco/computer_magic/nifty_utils/dbs2yeo/dk80toyeo7.mat'; %rsns

%% Recency analysis

pre=trophiccoherenceCT;
post=trophiccoherenceSZ;
change=(post-pre)';
subj=(1:24)';
X=table(post',pre',recency,ceremonies,dose,subj,'VariableNames',{'post','pre','recency','ceremonies','dose','subj'});
X=rmmissing(X);
mdl=fitlme(X,'post ~ recency*ceremonies + dose + pre + (1|subj)')

hfig=figure;
scatter(ceremonies,post,50,'blue','filled','o','MarkerEdgeColor','k');
hold on;
p=polyfit(ceremonies,post,1);
yfit=polyval(p,ceremonies);
plot(ceremonies,yfit,'r-','LineWidth',2);
xlabel('Number of ceremonies');
ylabel('Directedness');
R=corrcoef(ceremonies,post);
R2=R(1,2)^2;
text(mean(ceremonies),max(yfit),sprintf('R^2 = %.3f', R2),'FontSize',12)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5, 'FontSize', 12);

%% recency, post, recency
pre=trophiccoherenceCT;
post=trophiccoherenceSZ;
change=(post-pre)';
subj=(1:24)';
X=table(post',pre',recency,ceremonies,dose,subj,'VariableNames',{'post','pre','recency','ceremonies','dose','subj'});
X=rmmissing(X);
mdl=fitlme(X,'post ~ recency*ceremonies + dose + pre + (1|subj)')

hfig=figure;
scatter(recency,post,50,'blue','filled','o','MarkerEdgeColor','k');
hold on;
p=polyfit(recency,post,1);
yfit=polyval(p,recency);
plot(recency,yfit,'r-','LineWidth',2);
xlabel('Number of ceremonies');
ylabel('Directedness');
R=corrcoef(recency,post);
R2=R(1,2)^2;
text(mean(recency),max(yfit),sprintf('R^2 = %.3f', R2),'FontSize',12)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5, 'FontSize', 12);

%% recency, ceremonies (pre)
pre=trophiccoherenceCT;
post=trophiccoherenceSZ;
change=(post-pre)';
subj=(1:24)';
X=table(post',pre',recency,ceremonies,dose,subj,'VariableNames',{'post','pre','recency','ceremonies','dose','subj'});
X=rmmissing(X);
mdl=fitlme(X,'pre ~ recency*ceremonies + dose + (1|subj)')


hfig=figure;
scatter(X.ceremonies,X.pre,50,'blue','filled','o','MarkerEdgeColor','k');
hold on;
p=polyfit(X.ceremonies,X.pre,1);
yfit=polyval(p,X.ceremonies);
plot(X.ceremonies,yfit,'r-','LineWidth',2);
xlabel('Number of ceremonies');
ylabel('Directedness (baseline)');
R=corrcoef(X.ceremonies,X.pre);
R2=R(1,2)^2;
text(mean(X.ceremonies),max(yfit),sprintf('R^2 = %.3f', R2),'FontSize',12)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5, 'FontSize', 12);

%% recency, recency (pre)
pre=trophiccoherenceCT;
post=trophiccoherenceSZ;
change=(post-pre)';
subj=(1:24)';
X=table(post',pre',recency,ceremonies,dose,subj,'VariableNames',{'post','pre','recency','ceremonies','dose','subj'});
X=rmmissing(X);
mdl=fitlme(X,'pre ~ recency*ceremonies + dose + (1|subj)')

hfig=figure;
scatter(X.recency,X.pre,50,'blue','filled','o','MarkerEdgeColor','k');
hold on;
p=polyfit(X.recency,X.pre,1);
yfit=polyval(p,X.recency);
plot(X.recency,yfit,'r-','LineWidth',2);
xlabel('Time since last ceremony (days)');
ylabel('Directedness (baseline)');
R=corrcoef(X.recency,X.pre);
R2=R(1,2)^2;
text(mean(X.recency),max(yfit),sprintf('R^2 = %.3f', R2),'FontSize',12)
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.5, 'FontSize', 12);