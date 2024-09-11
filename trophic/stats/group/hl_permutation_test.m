%clear all;

%load results_Ceff_dmt.mat;
dbs80=importdata('/Users/myco/human_brains/atlases/dbs80/dbs80symm_labels.csv');

%% group change of hierarchy
N=80;
NSUB=24;

%eff1=CeffOccasionalSober;
%Ceff2=CeffOccasionalCannabis;
Ceff1=CeffCT;
Ceff2=CeffSZ;

%NSUB=28;
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
%NSUB=33;
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

for n=1:N % for each region
    SSND(n)=abs(mean(hierarchicallevels_sub_SZ(:,n))-mean(hierarchicallevels_sub_CT(:,n))) ...
        /sqrt(var(hierarchicallevels_sub_SZ(:,n))+var(hierarchicallevels_sub_CT(:,n))); %absolute value of mean diff between groups, normalized by the square root of the total variance of each group
end
[soSSND index]=sort(SSND,'descend');

for n=1:N %for each region
    clear nodelevelCT nodelevelSZ;
 %   NSUB=28;
    for nsub=1:NSUB
        gamma=hierarchicallevels_sub_CT(nsub,n);
        nodelevelCT(nsub,:)=gamma;
    end
  %  NSUB=33;
    for nsub=1:NSUB
        gamma=hierarchicallevels_sub_SZ(nsub,n);
        nodelevelSZ(nsub,:)=gamma;
    end

    a=mean(nodelevelCT,2)';
    b=mean(nodelevelSZ,2)';
    stats=permutation_htest2_np([b,a],[ones(1,numel(b)) 2*ones(1,numel(a))],10000,0.01,'ranksum');
    pp(n)=min(stats.pvals);
    t_idv(n)=stats.tvals;
    diffs_idv(n)=stats.diffs;
end
% pps=find(pp<0.01);
nsig=FDR_benjHoch(pp,0.2);

% minn=pps(end);
% n=minn
clear nodelevelCT nodelevelSZ;
%SUB=28;
for nsub=1:NSUB
    gamma=hierarchicallevels_sub_CT(nsub,nsig);
    nodelevelCT(nsub,:)=gamma;
end
%NSUB=33;
for nsub=1:NSUB
    gamma=hierarchicallevels_sub_SZ(nsub,nsig);
    nodelevelSZ(nsub,:)=gamma;
end

a=mean(nodelevelCT,2)';
b=mean(nodelevelSZ,2)';
stats=permutation_htest2_np([b,a],[ones(1,numel(b)) 2*ones(1,numel(a))],10000,0.01,'ranksum');
pvalue=min(stats.pvals)

figure(1)
boxplot([a b],[zeros(1,length(a)) ones(1,length(b))]);

Region=repmat(nsig,2,1);
Condition=[repelem("Baseline",length(nsig),1); repelem("DMT",length(nsig),1)];
HL=[mean(nodelevelCT)'; mean(nodelevelSZ)'];
t=table(Region,Condition,HL)

%save results_hlRSN_chr.mat t stats t_idv diffs_idv nsig

%% render

notsig=setdiff(1:80,nsig);
hlChange=t_idv;
%hlChange=nanmean(hierarchicallevels_sub_SZ)-nanmean(hierarchicallevels_sub_CT);
hlChange(notsig)=0;
rendersurface_dbs80(hlChange,-5,5,1,'RdBu4')


