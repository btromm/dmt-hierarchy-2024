%% DMT
clear all;
path2data='/Users/myco/human_brains/dmt_update';
cd(path2data);
load results_Ceff_dmt_update;
load results_DMT_modelfree.mat
numSub=17;
results_path='/Users/myco/human_brains/INSIDEOUT_THESIS/DMT_update';

%% render brains
rendersurface_dbs80(mean(G_inPCB_Post),0,2,0)
rendersurface_dbs80(mean(G_outPCB_Post),0,2,0)
rendersurface_dbs80(mean(G_totalPCB_Post),0,3.5,0)

rendersurface_dbs80(mean(G_inDMT_Post),0,2,0)
rendersurface_dbs80(mean(G_outDMT_Post),0,2,0)
rendersurface_dbs80(normalize(mean(G_totalDMT_Post),'range'),0,1,0)

%differences
rendersurface_dbs80(mean(G_totalDMT_Post-G_totalPCB_Post),-1.75,1.75,1,'RdBu11')

mnhlCT=nanmean(hierarchicallevelsPCB_Post);
mnhlSZ=normalize(nanmean(hierarchicallevelsDMT_Post),'range');
rendersurface_dbs80(mnhlSZ,0,1,0)
%% GRAPH EXPORT
gCT=squeeze(mean(CeffPCB_Post))';
gCT(find(gCT<0.015))=0;
gCT(:,32:49)=[];
gCT(32:49,:)=[];
gCT=digraph(gCT);
writetable(gCT.Edges,'PCBEdges.csv')
hlCT=mean(hierarchicallevelsPCB_Post)';
hlCT(32:49)=[];
writematrix(hlCT,'hlPCB.csv');

gSZ=squeeze(mean(CeffDMT_Post))';
gSZ(find(gSZ<0.015))=0;
gSZ(:,32:49)=[];
gSZ(32:49,:)=[];
gSZ=digraph(gSZ);
writetable(gSZ.Edges,'DMTEdges.csv');
hlSZ=mean(hierarchicallevelsDMT_Post)';
hlSZ(32:49)=[];
writematrix(hlSZ,'hlDMT.csv');

%% modelfree comparison (*)

pre=[FowRev_pcbpre'; FowRev_dmtpre'];
post=[FowRev_pcbpost'; FowRev_dmtpost'];
group=categorical([repelem("Placebo",17,1); repelem("DMT",17,1)]);
subj=repmat(1:17,1,2)';
X=table(pre,post,group,subj,'VariableNames',{'pre','post','group','subj'});
X.group=categorical(X.group);
mdl=fitlme(X,'post ~ group+pre + (1 | subj)')

%% tc comparison (*)

pre=[trophiccoherencePCB_Pre'; trophiccoherenceDMT_Pre'];
post=[trophiccoherencePCB_Post'; trophiccoherenceDMT_Post'];
group=categorical([repelem("Placebo",17,1); repelem("DMT",17,1)]);
subj=repmat(1:17,1,2)';
X=table(pre,post,group,subj,'VariableNames',{'pre','post','group','subj'});
X.group=categorical(X.group);
mdl=fitlme(X,'post ~ group+pre + (1 | subj)')


%% plot

% Plotting
figure; % create a new figure
hold on; % allow multiple plots on the same figure

% Generate x-values for pre and post measurements
xPre = repmat(1, numSub, 1);
xPost = repmat(2, numSub, 1);

for i = 1:numSub
    % Plotting for Placebo
    plot([xPre(i), xPost(i)], [trophiccoherencePCB_Pre(i), trophiccoherencePCB_Post(i)], '-o', 'Color', [0, 0.5, 0.5], 'DisplayName', ['Subject ' num2str(i) ' Placebo']);
    
    % Plotting for DMT
    plot([xPre(i), xPost(i)], [trophiccoherenceDMT_Pre(i), trophiccoherenceDMT_Post(i)], '-x', 'Color', [0.5, 0, 0.5], 'DisplayName', ['Subject ' num2str(i) ' DMT']);
end

% Formatting the plot
xticks([1, 2]);
xticklabels({'Pre', 'Post'});
xlabel('Treatment Timepoint');
ylabel('Trophic Coherence');
legend('off'); % To avoid overcrowding of legend items
title('Trophic Coherence Before and After Treatment');
grid on; % adds a grid for better readability
hold off; % stop plotting on the same figure




%% Hierarchical Levels RSN
%{
yeo7{1}='Vis';
yeo7{2}='SomMot';
yeo7{3}='DorsAttn';
yeo7{4}='SalVentAttn';
yeo7{5}='Limbic';
yeo7{6}='FPN';
yeo7{7}='Default';
%}

yeohl_CT = zeros(numSub,7);
yeohl_CT_tally=zeros(numSub,7);
yeohl_SZ = zeros(numSub,7);
yeohl_SZ_tally=zeros(numSub,7);

for i=1:numSub
    for j=1:80
        if av(j) ~=0 % make sure it's part of an RSN
            yeohl_CT(i,av(j))=yeohl_CT(i,av(j))+hierarchicallevelsCT(i,j);
            yeohl_CT_tally(i,av(j))=yeohl_CT_tally(i,av(j))+1;

            yeohl_SZ(i,av(j))=yeohl_SZ(i,av(j))+hierarchicallevelsSZ(i,j);
            yeohl_SZ_tally(i,av(j))=yeohl_SZ_tally(i,av(j))+1;
        end
    end
end

yeohl_CT = yeohl_CT ./ yeohl_CT_tally; %normalize to mean
yeohl_SZ = yeohl_CT ./ yeohl_SZ_tally;

hierarchicalRSN=[yeohl_CT yeohl_SZ]; % to R!
writematrix(hierarchicalRSN);

%% top 5 regions
mnCT=nanmean(hierarchicallevelsCT);
mnSZ=nanmean(hierarchicallevelsSZ);
% find top 5 regions
[~,top5CT]=maxk(mnCT,5);
[~,top5SZ]=maxk(mnSZ,5);

% find top 5 for each subject
for i=1:24
    [~,top5CT(i,:)]=maxk(hierarchicallevelsCT(i,:),5);
    [~,top5SZ(i,:)]=maxk(hierarchicallevelsSZ(i,:),5);
end


%% svm
concatlevels = [hierarchicallevelsCT; hierarchicallevelsSZ];
labels=[zeros(24,1); ones(24,1)];
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

    svm_model=fitcsvm(train_data,train_labels);

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
figure;
confusionchart(confMat, {'Group 1', 'Group 2'});
title('Confusion Matrix (Heatmap)');


save results_Ceff_aya.mat CeffCT CeffSZ  ...
     hierarchicallevelsCT hierarchicallevelsSZ trophiccoherenceSZ trophiccoherenceCT;


% You want to use SVM on both the GEC and then the trophic coherence. This
% gives more detail through trophic coherence, but also may show that
% classifying conditions is more accurate using trophic than GEC alone

%% appendix
%% show fitted matrices -> measure of gof to empirical
figure(1);
subplot(1,2,1);
boxplot([fittFC_PCB_Post' fittFC_DMT_Post']);
title('Fitted FC Matrix');
subplot(1,2,2);
boxplot([fittCVtau_PCB_Post' fittCVtau_DMT_Post']);
title('Fitted COV matrix');