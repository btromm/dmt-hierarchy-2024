%% trophic_stats_cannabis
clear all;
path2data='/Users/myco/human_brains/cannab';
cd(path2data)
load results_Ceff_cannabis.mat;
%% render brains

%chronic
rendersurface_dbs80(normalize(mean(G_totalChronicCannabis),'range'),0,1,0);

%occasional
rendersurface_dbs80(normalize(mean(G_totalOccasionalCannabis),'range'),0,1,0);

% diff

ccnorm=normalize(nanmean(G_totalChronicCannabis),'range');
ocnorm=normalize(nanmean(G_totalOccasionalCannabis),'range');
diff=normalize(ccnorm-ocnorm,'range',[-1 1]);
rendersurface_dbs80(diff,-1,1,1,'RdBu11')

%hl
mnhlCC=nanmean(hierarchicallevelsChronicCannabis);
mnhlOC=nanmean(hierarchicallevelsOccasionalCannabis);

rendersurface_dbs80(normalize(mnhlCC,'range'),0,1,0)
rendersurface_dbs80(normalize(mnhlOC,'range'),0,1,0)

diff=normalize(mnhlCC-mnhlOC,'range',[-1 1]);
rendersurface_dbs80(diff,-1,1,1,'RdBu11');


%% GRAPH EXPORT TO R
gCC=squeeze(mean(CeffChronicCannabis))';
gCC(find(gCC<0.015))=0;
gCC(:,32:49)=[];
gCC(32:49,:)=[];
gCC=digraph(gCC);
writetable(gCC.Edges,'CCEdges.csv')
hlCC=nanmean(hierarchicallevelsChronicCannabis)';
hlCC(32:49)=[];
writematrix(hlCC,'hlCC.csv');

gCS=squeeze(mean(CeffChronicSober))';
gCS(find(gCS<0.015))=0;
gCS(:,32:49)=[];
gCS(32:49,:)=[];
gCS=digraph(gCS);
writetable(gCS.Edges,'CSEdges.csv');
hlCS=nanmean(hierarchicallevelsChronicSober)';
hlCS(32:49)=[];
writematrix(hlCS,'hlCS.csv');

gOC=squeeze(mean(CeffOccasionalCannabis))';
gOC(find(gOC<0.015))=0;
gOC(:,32:49)=[];
gOC(32:49,:)=[];
gOC=digraph(gOC);
writetable(gOC.Edges,'OCEdges.csv')
hlOC=nanmean(hierarchicallevelsOccasionalCannabis)';
hlOC(32:49)=[];
writematrix(hlOC,'hlOC.csv');

gOS=squeeze(mean(CeffOccasionalSober))';
gOS(find(gOS<0.015))=0;
gOS(:,32:49)=[];
gOS(32:49,:)=[];
gOS=digraph(gOS);
writetable(gOS.Edges,'OSEdges.csv');
hlOS=nanmean(hierarchicallevelsOccasionalSober)';
hlOS(32:49)=[];
writematrix(hlOS,'hlOS.csv');


%% svm classification for trophic levels
figure(3);

grp_cannabis=[repelem("Occasional/Baseline",80,1);repelem("Chronic/Baseline",80,1);repelem("Occasional/Cannabis",80,1);repelem("Chronic/Cannabis",80,1)];
violinplot([mean(hierarchicallevelsSoberOccasional)' mean(hierarchicallevelsSoberChronic)' mean(hierarchicallevelsCannabisOccasional)' mean(hierarchicallevelsCannabisChronic)'],grp_cannabis);
concatlevels = [hierarchicallevelsSoberOccasional; hierarchicallevelsSoberChronic; hierarchicallevelsCannabisOccasional; hierarchicallevelsCannabisChronic];
mnSoberOccasional=mean(hierarchicallevelsSoberOccasional);
mnSoberChronic=mean(hierarchicallevelsSoberChronic);
mnCannabisOccasional=mean(hierarchicallevelsCannabisOccasional);
mnCannabisChronic=mean(hierarchicallevelsCannabisChronic);
[~,top5SO]=maxk(mnSoberOccasional,5);
[~,top5SC]=maxk(mnSoberChronic,5);
[~,top5CO]=maxk(mnCannabisOccasional,5);
[~,top5CC]=maxk(mnCannabisChronic,5);

labels=[zeros(39,1); ones(27,1); ones(38,1).*2;ones(33,1).*3];
numObs=length(labels);
folds=floor(5);
c=cvpartition(labels,'KFold',folds);

for i=1:folds
    train_idx=c.training(i);
    test_idx=c.test(i);

    train_data=concatlevels(train_idx,:);
    train_labels=labels(train_idx);

    test_data=concatlevels(test_idx,:);
    test_labels=labels(test_idx);

    svm_model=fitcecoc(train_data,train_labels);

    [predictions,scores]=predict(svm_model,test_data);
    if ~exist('prediction','var')
        prediction = predictions;
    else
        prediction=[prediction; predictions];    
    end
   

    accuracy=sum(predictions==test_labels)/numel(test_labels);
    accuracies(i)=accuracy;
end

avg_accuracy=mean(accuracies);
disp(['Average Accuracy (',num2str(folds),'-fold cross-validation): ', num2str(avg_accuracy)])

% Compute confusion matrix
confMat = confusionmat(test_labels, predictions);

% Display the confusion matrix as a heatmap
figure(4);
confusionchart(confMat, {'Sober/Occasional', 'Sober/Chronic','Cannabis/Occasional','Cannabis/Chronic'});
title('Confusion Matrix (Heatmap)');

% You want to use SVM on both the GEC and then the trophic coherence. This
% gives more detail through trophic coherence, but also may show that
% classifying conditions is more accurate using trophic than GEC alone


%% regional level t-tests
%hierarchicallevels_baseline
%{
alpha=0.05/80;
for i=1:80
    [h,p]=ttest(hierarchicallevelsCT(:,i),hierarchicallevelsSZ(:,i));
    corrs(i,1)=corr(hierarchicallevelsCT(:,i),hierarchicallevelsSZ(:,i));
    tests(i,2)=p;
    if tests(i,2)<alpha
        tests(i,1)=1;
    else
        tests(i,1)=0;
    end
end
%}

%% Appendix
%% show fitted matrices
figure(2);
subplot(1,2,1);
boxplot([fittFC_SoberOccasional'; fittFC_CannabisOccasional'; fittFC_SoberChronic'; fittFC_CannabisChronic'],grp_cannabis);
title('Fitted FC Matrix');
subplot(1,2,2);
boxplot([fittCVtau_SoberOccasional'; fittCVtau_CannabisOccasional'; fittCVtau_SoberChronic'; fittCVtau_CannabisChronic'],grp_cannabis);
title('Fitted COV matrix');