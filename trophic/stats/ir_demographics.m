%% trophic coherence
% correlations b/w changes in irreversibility and demographics...
    % @kat: this script failed entirely but maybe u can make it work :3

load demographics.mat;
load results_aya_modelfree.mat;

ir=[hierarchy_sober';hierarchy_ayahuasca'];
subj=repmat(1:24,1,2)';
dose2=repmat(dose,2,1);
group=[ones(24,1);2*ones(24,1)];
X=table(ir,dose2,subj,group,'VariableNames',{'ir','dose','subj','group'});
X=rmmissing(X);

mdl = fitlme(X,'ir ~ dose + group + (1|subj)');

%% trophic coherence
load results_Ceff_aya.mat;

tc=[trophiccoherenceCT';trophiccoherenceSZ'];
subj=repmat(1:24,1,2)';
dose2=repmat(dose,2,1);
group=[ones(24,1);2*ones(24,1)];
X2=table(tc,dose2,subj,group,'VariableNames',{'tc','dose','subj','group'});

mdl2 = fitlme(X2,'tc ~ dose + group + (1|subj)');

%% HL

a=mean(nodelevelCT,2)';
b=mean(nodelevelSZ,2)';
HL=[a;b];
subj=repmat(1:24,1,2)';
dose2=repmat(dose,2,1);
group=[ones(24,1);2*ones(24,1)];
X3=table(HL,dose2,subj,group,'VariableNames',{'HL','dose','subj','group'});
mdl3 = fitlme(X3,'HL ~ dose + group + (1|subj)');