dbs80=importdata('/Users/myco/data/atlases/dbs80/dbs80symm_labels.csv');
path2yeo='/Users/myco/code/nifty_utils/dbs2yeo';
addpath(path2yeo);
load dk80toyeo7.mat;
av(32:49)=8;

%% group change of Hierarchy
N=7;
NSUB=24;

Ceff1=CeffCT;
Ceff2=CeffSZ;

% Calculate Hierarchical levels
for nsub=1:NSUB
    Ceff=squeeze(Ceff1(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub_CT(nsub,:)=gamma';
end
for nsub=1:NSUB
    Ceff=squeeze(Ceff2(nsub,:,:));
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;
    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevels_sub_SZ(nsub,:)=gamma';
end

hlyeobsl=decompose_rsn(hierarchicallevels_sub_CT,av);
hlyeodrug=decompose_rsn(hierarchicallevels_sub_SZ,av);



%% stats
for n=1:N % for each region
    SSND(n)=abs(mean(hlyeodrug(:,n))-mean(hlyeobsl(:,n))) ...
        /sqrt(var(hlyeodrug(:,n))+var(hlyeobsl(:,n))); %absolute value of mean diff between groups, normalized by the square root of the total variance of each group
end
[soSSND index]=sort(SSND,'descend');

% permutation test for each region
for n=1:N %for each region
    clear nodelevelCT nodelevelSZ;
    for nsub=1:NSUB
        gamma=hlyeobsl(nsub,n);
        nodelevelCT(nsub,:)=gamma;
    end
    for nsub=1:NSUB
        gamma=hlyeodrug(nsub,n);
        nodelevelSZ(nsub,:)=gamma;
    end

    a=nanmean(nodelevelCT,2)';
    b=nanmean(nodelevelSZ,2)';
    stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
    t_idv(n)=stats.tvals;
    diffs_idv(n)=stats.diffs;
end
% pps=find(pp<0.01);
nsig=FDR_benjHoch(pp,0.2); % find significant regions after FDR correction
ps=mafdr(pp)
alpha=0.05/7;
H=find(pp<alpha)

% minn=pps(end);
% n=minn
clear nodelevelCT nodelevelSZ;
for nsub=1:NSUB
    gamma=hlyeobsl(nsub,nsig);
    nodelevelCT(nsub,:)=gamma;
end
for nsub=1:NSUB
    gamma=hlyeodrug(nsub,nsig);
    nodelevelSZ(nsub,:)=gamma;
end

% permutation test over all significant regions (total diff)
a=nanmean(nodelevelCT,2)';
b=nanmean(nodelevelSZ,2)';
stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.01,'ranksum'); %10k permutations, p < 0.01, non-parametric Wilcoxon r
pvalue=min(stats.pvals)

figure(1)
boxplot([a b],[zeros(1,length(a)) ones(1,length(b))]);

Region=repmat(nsig,2,1);
Condition=[repelem("Baseline",length(nsig),1); repelem("Ayahuasca",length(nsig),1)];
HL=[mean(nodelevelCT)'; mean(nodelevelSZ)'];
t=table(Region,Condition,HL)

save results_hlRSN_aya.mat t stats t_idv

%writematrix(repmat(nsig,2,1),'nsigrsn.csv')
%cond=[repelem("Baseline",7,1); repelem("Ayahuasca",7,1)];
%writematrix(cond,'conditionrsn.csv')
%writematrix([mean(nodelevelCT)'; mean(nodelevelSZ)'],'sighlrsn.csv')


