%% this part is only for figure 1
%% effective connectivity (for figures)

%% Loading stuff
clear all;
path2data='/Users/myco/human_brains/ayahuasca';
cd(path2data);
load results_Ceff_aya;
load results_aya_modelfree.mat
load demographics.mat;
numSub=24;
results_path='/Users/myco/human_brains/INSIDEOUT_THESIS';
load '/Users/myco/Code/nifty_utils/dbs2yeo/dk80toyeo7.mat'; %rsns

mnCeffCT=squeeze(mean(CeffCT));
mnCeffSZ=squeeze(mean(CeffSZ));
figure(1);
imagesc(mnCeffCT);
figure(2);
imagesc(mnCeffSZ);

%% NR matrices
mnNRCT=squeeze(mean(Tenet2_R));
mnNRSZ=squeeze(mean(Tenet2_A));
figure(3);
imagesc(mnNRCT);
figure(4);
imagesc(mnNRSZ);

%% FC Matrices
load aya_ses01_hopf.mat
figure(5);
imagesc(FCemp);
load aya_ses02_hopf.mat
figure(6);
imagesc(FCemp);

%% main part
%% render brains
rendersurface_dbs80(mean(G_inCT),0,1.5,0)
rendersurface_dbs80(mean(G_outCT),0,1.5,0)
rendersurface_dbs80(mean(G_totalCT),0,2,0)

rendersurface_dbs80(mean(G_inSZ),0,1.5,0)
rendersurface_dbs80(mean(G_outSZ),0,1.5,0)
rendersurface_dbs80(normalize(mean(G_totalSZ),'range'),0,1,0)


rendersurface_dbs80(mean(G_totalSZ-G_totalCT),-.75,.75,1,'RdBu11')

mnhlCT=nanmean(hierarchicallevelsCT);
mnhlSZ=normalize(nanmean(hierarchicallevelsSZ),'range');
diff=normalize(mnhlSZ-mnhlCT,'range',[-1 1]);
rendersurface_dbs80(mnhlSZ,0,1,0);

%% Hierarchical Levels RSN
%{
yeo7{1}='Vis';
yeo7{2}='SomMot';
yeo7{3}='DorsAttn';
yeo7{4}='SalVentAttn';
yeo7{5}='Limbic';
yeo7{6}='Cont';
yeo7{7}='Default';
%}

yeoCT=decompose_rsn(hierarchicallevelsCT,av);
yeoSZ=decompose_rsn(hierarchicallevelsSZ,av);

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


%% appendix
%% show fitted matrices -> measure of gof to empirical
figure(1);
subplot(1,2,1);
boxplot([fittFC_CT' fittFC_SZ']);
title('Fitted FC Matrix');
subplot(1,2,2);
boxplot([fittCVtau_CT' fittCVtau_SZ']);
title('Fitted COV matrix');