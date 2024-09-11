%% GROUP RENDERS
clear all;
path2data='/Users/myco/human_brains/ayahuasca';
cd(path2data);
load results_Ceff_aya;

G_totalBSL=G_totalCT;
G_totalAYA=G_totalSZ;

path2data='/Users/myco/human_brains/dmt_update';
cd(path2data);
load results_Ceff_DMT;

G_totalPCB=G_totalCT;
G_totalDMT=G_totalSZ;

path2data='/Users/myco/human_brains/cannab';
cd(path2data);
load results_Ceff_cannabis.mat;



%% aya vs dmt
% diffs
ayanorm=normalize(mean(G_totalAYA),'range');

dmtnorm=normalize(mean(G_totalDMT),'range');
diff_dmt=ayanorm-dmtnorm;
rendersurface_dbs80(diff_dmt,-1,1,1,'RdBu4');

% similarities (intersect)
[~,idx1]=maxk(ayanorm,20);
[~,idx2]=maxk(dmtnorm,20);
intr=intersect(idx1,idx2); %indexes that intersect

ayanorm(setdiff(1:end,intr))=0;
dmtnorm(setdiff(1:end,intr))=0;
sim_aya=mean([ayanorm; dmtnorm]);
rendersurface_dbs80(normalize(sim_aya,'range',[0 1]),0,1,0);

%% chronic vs. occasional
ccnorm=normalize(nanmean(G_totalChronicCannabis),'range');
ocnorm=normalize(nanmean(G_totalOccasionalCannabis),'range');
diff_cannab=normalize(ccnorm-ocnorm,'range',[-1 1]);
rendersurface_dbs80(diff_cannab,-1,1,1,'RdBu11')

% similarities (intersect)
[~,idx1]=maxk(ccnorm,20);
[~,idx2]=maxk(ocnorm,20);
intr=intersect(idx1,idx2); %indexes that intersect

ccnorm(setdiff(1:end,intr))=0;
ocnorm(setdiff(1:end,intr))=0;
sim_cannab=mean([ccnorm; ocnorm]);
rendersurface_dbs80(normalize(sim_cannab,'range',[0 1]),0,1,0);

%% (aya v dmt) v (chronic v occ)
diff=normalize(diff_dmt-diff_cannab,'range',[-1 1]);
rendersurface_dbs80(diff,-1,1,1,'RdBu11');

% intersect
[~,idx1]=maxk(sim_aya,80);
[~,idx2]=maxk(sim_cannab,80);
intr=intersect(idx1,idx2);

sim_aya(setdiff(1:end,intr))=0;
sim_cannab(setdiff(1:end,intr))=0;
sim=mean([sim_aya; sim_cannab]);
rendersurface_dbs80(normalize(sim,'range',[0 1]),0,1,0);

%% HL
clear all;
path2data='/Users/myco/human_brains/ayahuasca';
cd(path2data);
load results_Ceff_aya;
load /Users/myco/human_brains/ayahuasca/results_hlRSN_dmt.mat;
NSIG_AYA=nsig;


HLBSL=hierarchicallevelsCT;
HLAYA=hierarchicallevelsSZ;

path2data='/Users/myco/human_brains/dmt_update';
cd(path2data);
load results_Ceff_dmt_update;
load /Users/myco/human_brains/dmt_update/results_hlRSN_dmt.mat;
NSIG_DMT=nsig;


HLPCB=hierarchicallevelsCT;
HLDMT=hierarchicallevelsSZ;

path2data='/Users/myco/human_brains/cannab';
cd(path2data);
load results_Ceff_cannabis.mat;
load results_hltest_chr.mat;
NSIG_CHR=nsig;
load results_hltest_occ.mat;
NSIG_OCC=nsig;

%% aya vs dmt
% diffs
ayanorm=normalize(nanmean(HLAYA),'range');
dmtnorm=normalize(nanmean(HLDMT),'range');
diff_dmt=ayanorm-dmtnorm;
diff_dmt(setdiff(1:80,NSIG_DMT))=0;
diff_dmt(setdiff(1:80,NSIG_AYA))=0;

%diff_dmt=nanmean(HLAYA)-nanmean(HLDMT);
rendersurface_dbs80(diff_dmt,-1,1,1,'RdBu4');

% similarities (intersect)
ayanorm=normalize(nanmean(HLAYA),'range');
dmtnorm=normalize(nanmean(HLDMT),'range');
[~,idx1]=maxk(ayanorm,40);
[~,idx2]=maxk(dmtnorm,40);
intr=intersect(idx1,idx2); %indexes that intersect

ayanorm(setdiff(1:end,intr))=0;
dmtnorm(setdiff(1:end,intr))=0;
sim_aya=mean([ayanorm; dmtnorm]);
rendersurface_dbs80(normalize(sim_aya,'range',[0 1]),0,1,0);

%% chronic vs. occasional
ccnorm=normalize(nanmean(hierarchicallevelsChronicCannabis),'range');
ocnorm=normalize(nanmean(hierarchicallevelsOccasionalCannabis),'range');
diff_cannab=normalize(ccnorm-ocnorm,'range',[-1 1]);
rendersurface_dbs80(diff_cannab,-1,1,1,'RdBu11')

% similarities (intersect)
ccnorm=normalize(nanmean(hierarchicallevelsChronicCannabis),'range');
ocnorm=normalize(nanmean(hierarchicallevelsOccasionalCannabis),'range');
[~,idx1]=maxk(ccnorm,40);
[~,idx2]=maxk(ocnorm,40);
intr=intersect(idx1,idx2); %indexes that intersect

ccnorm(setdiff(1:end,intr))=0;
ocnorm(setdiff(1:end,intr))=0;
sim_cannab=mean([ccnorm; ocnorm]);
rendersurface_dbs80(normalize(sim_cannab,'range',[0 1]),0,1,0);

%% (aya v dmt) v (chronic v occ)
diff=normalize(diff_dmt-diff_cannab,'range',[-1 1]);
rendersurface_dbs80(diff,-1,1,1,'RdBu11');

% intersect
[~,idx1]=maxk(sim_aya,80);
[~,idx2]=maxk(sim_cannab,80);
intr=intersect(idx1,idx2);

sim_aya(setdiff(1:end,intr))=0;
sim_cannab(setdiff(1:end,intr))=0;
sim=mean([sim_aya; sim_cannab]);
rendersurface_dbs80(normalize(sim,'range',[0 1]),0,1,0);

%% regional changes in hierarchy
clear all;
path2yeo='/Users/myco/computer_magic/nifty_utils/dbs2yeo';
addpath(path2yeo);
load dk80toyeo7.mat;
av(32:49)=8;

path2data='/Users/myco/human_brains/ayahuasca';
cd(path2data);
load results_Ceff_aya;

hlBSL=hierarchicallevelsCT;
hlAYA=hierarchicallevelsSZ;

path2data='/Users/myco/human_brains/dmt_update';
cd(path2data);
load results_Ceff_DMT;

hlPCB=hierarchicallevelsCT;
hlDMT=hierarchicallevelsSZ;

path2data='/Users/myco/human_brains/cannab';
cd(path2data);
load results_Ceff_cannabis.mat;

hlBSL=decompose_rsn(hlBSL,av);
hlAYA=decompose_rsn(hlAYA,av);
hlPCB=decompose_rsn(hlPCB,av);
hlDMT=decompose_rsn(hlDMT,av);
hlCC=decompose_rsn(hierarchicallevelsChronicCannabis,av);
hlCS=decompose_rsn(hierarchicallevelsChronicSober,av);
hlOC=decompose_rsn(hierarchicallevelsOccasionalCannabis,av);
hlOS=decompose_rsn(hierarchicallevelsOccasionalSober,av);

rsn_aya=mean(hlAYA-hlBSL);
rsn_dmt=mean(hlDMT-hlPCB);
rsn_chronic=mean(hlCC)-mean(hlCS);
rsn_occ=mean(hlOC-hlOS);

writematrix(normalize([rsn_aya; rsn_dmt; rsn_chronic; rsn_occ],'range',[-1 1]),'/Users/myco/human_brains/LSDAyaCannab/rsnhlsigchanges.csv');
%% SVM HL
% load data
clear all;
path2data='/Users/myco/human_brains/msc_thesis/ayahuasca';
cd(path2data);
load results_Ceff_aya;

hlBSL=normalize(hierarchicallevelsCT,'range');
hlAYA=normalize(hierarchicallevelsSZ,'range');

path2data='/Users/myco/human_brains/dmt_update';
cd(path2data);
load results_Ceff_DMT;

hlPCB=normalize(hierarchicallevelsCT,'range');
hlDMT=normalize(hierarchicallevelsSZ,'range');

path2data='/Users/myco/human_brains/cannab';
cd(path2data);
load results_Ceff_cannabis.mat;

hierarchicallevelsChronicCannabis=normalize(hierarchicallevelsChronicCannabis,'range');
hierarchicallevelsChronicSober=normalize(hierarchicallevelsChronicSober,'range');
hierarchicallevelsOccasionalCannabis=normalize(hierarchicallevelsOccasionalCannabis,'range');
hierarchicallevelsOccasionalSober=normalize(hierarchicallevelsOccasionalSober,'range');

% organize data
concatlevels=[hlBSL; hlAYA; hlPCB; hlDMT; hierarchicallevelsChronicCannabis; hierarchicallevelsChronicSober; hierarchicallevelsOccasionalCannabis; hierarchicallevelsOccasionalSober];
labels=[repelem("Baseline",24,1); repelem("Ayahuasca",24,1); repelem("Placebo",16,1); repelem("DMT",16,1);...
    repelem("Chronic/Cannabis",33,1); repelem("Chronic/Sober",28,1); repelem("Occasional/Cannabis",41,1); repelem("Occasional/Sober",41,1)];

%fit svm
rng default
Mdl = fitcecoc(concatlevels,labels,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))

p=crossval(Mdl,'KFold',5);
genError = kfoldLoss(p);
accuracy=1-genError
l=kfoldPredict(p);

confMat = confusionmat(l, labels);

% Display the confusion matrix as a heatmap
figure;
confusionchart(confMat,{'Baseline/Ayahuasca','Ayahuasca','Placebo','DMT','Chronic/Cannabis','Chronic/Sober','Occasional/Cannabis','Occasional/Sober'});
title('Confusion Matrix (Heatmap)');

% You want to use SVM on both the GEC and then the trophic coherence. This
% gives more detail through trophic coherence, but also may show that
% classifying conditions is more accurate using trophic than GEC alone

%% CEFF SVM -> all bounded b/w 0 & 0.200
clear all;
path2data='/Users/myco/human_brains/msc_thesis/ayahuasca';
cd(path2data);
load results_aya_modelfree;

for i=1:24
    CeffBSL(i,:)=Tenet2_R(i,:);
    CeffAYA(i,:)=Tenet2_A(i,:);
end
CeffBSL=normalize(CeffBSL,'range');
CeffAYA=normalize(CeffAYA,'range');
path2data='/Users/myco/human_brains/dmt_update';
cd(path2data);
load results_dmt_modelfree;
for i=1:16
    CeffPCB(i,:)=Tenet2_R(i,:);
    CeffDMT(i,:)=Tenet2_A(i,:);
end
CeffPCB=normalize(CeffPCB,'range');
CeffDMT=normalize(CeffDMT,'range');
path2data='/Users/myco/human_brains/cannab';
cd(path2data);
load results_cannabischr_modelfree;

for i=1:28
    CeffCS(i,:)=Tenet2_R(i,:);
end
for i=1:33
    CeffCC(i,:)=Tenet2_A(i,:);
end
CeffCS=normalize(CeffCS,'range');
CeffCC=normalize(CeffCC,'range');
load results_cannabisocc_modelfree
for i=1:41
    CeffOS(i,:)=Tenet2_R(i,:);
    CeffOC(i,:)=Tenet2_A(i,:);
end
CeffOS=normalize(CeffOS,'range');
CeffOC=normalize(CeffOC,'range');

concatlevels = [CeffBSL; CeffAYA; CeffPCB; CeffDMT; CeffCC; CeffCS; CeffOC; CeffOS];
labels=[repelem("Baseline",24,1); repelem("Ayahuasca",24,1); repelem("Placebo",16,1); repelem("DMT",16,1);...
    repelem("Chronic/Cannabis",33,1); repelem("Chronic/Sober",28,1); repelem("Occasional/Cannabis",41,1); repelem("Occasional/Sober",41,1)];

rng default
Mdl = fitcecoc(concatlevels,labels,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus'))

p=crossval(Mdl,'KFold',5);
genError = kfoldLoss(p);
accuracy=1-genError
l=kfoldPredict(p);

confMat = confusionmat(l, labels);

% Display the confusion matrix as a heatmap
figure;
confusionchart(confMat,{'Baseline/Ayahuasca','Ayahuasca','Placebo','DMT','Chronic/Cannabis','Chronic/Sober','Occasional/Cannabis','Occasional/Sober'});
title('Confusion Matrix (Heatmap)');
set(gca,'fontname','SF Pro');
set(gca,'fontsize',18);