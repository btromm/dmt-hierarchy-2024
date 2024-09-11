clear all;
path2data='/Users/myco/human_brains/msc_thesis/ayahuasca';
path2code='/Users/myco/computer_magic/trophic_coherence';
addpath(path2code);
cd(path2data);
addpath(path2data);
N=80; %regions

indexN=1:N;
Tau=4;
sigma=0.01; %noise parameter

% learning parameters - control rate of convergence & magnitude of updates
% to connectivity matrix during optimization process.
epsFC=0.0004; %learning rate for adjusting FC based on diff b/w empirical and functional connectivity
epsFCtau=0.0001; %controls learning for adjusting connectivity matrix based on diff b/w empirical and simulated covariance at time delay
maxC=0.2; %maximum value of the connectivity matrix (C) and initial anatomical connectivity

load SC_dbs80HARDIFULL.mat;

Isubdiag = find(tril(ones(N),-1)); %all indexes of lower triangular portion

% Parameters of the data
TR=1.4;  % Repetition Time (seconds)

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

% Load structural connectivity from anatomical DTI tractography
C = SC_dbs80HARDI;
C = C/max(max(C))*maxC; %scales C by 0 to maxC

%%
%% AYA_SES01 (CT) BASELINE
load aya_ses01.mat; %original data
load aya_ses01_hopf.mat; %power spectrum & FC
f_diff_c=f_diff;
NSUB=24;

%% Group
for nsub=1:NSUB
    ts=subject{nsub}.dbs80ts;  % fMRI CT
    %ts=ts';

    % Preprocessing
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:))); %detrend timeseries by subtracting mean value from detrended timeseries. removes linear trends
    end

    ts2=ts(indexN,10:end-10); %removes first and last ten time points
    Tm=size(ts2,2); % # timepoints
    FCemp=corrcoef(ts2'); % C matrix
    FCCT(nsub,:,:)=FCemp; % store FC matrices for all subjects
    COVemp=cov(ts2'); %covariance matrix of the timeseries
    % COV(tau)
    tst=ts2';
    for i=1:N %for each row
        for j=1:N % for each column
            % in a cov matrix, the square root of an element IS its
            % standard deviation -> diagonal elements represent variances
            % of variables, and those are the square of the standard
            % deviations
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j)); %calculates reciprocal of product of standard deviations of i-th and j-th brain regions, stores in sigratio matrix
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau); %cross-covariance b/w timeseries of i-th and j-th elements with Tau timelag
            indx=find(lags==Tau); %find Tau in the cross-covariance test
            COVtauemp(i,j)=clag(indx)/size(tst,1); % assign clag divided by timepoints to i-th j-th region comparison
        end
    end
    % this is the NR matrix
    COVtauemp=COVtauemp.*sigratio; %scale elements of COVtauemp (the covariance matrix of FC at Tau timelag) by the sigratio -> adjust covariances based on standard deviations of brain regions
    COVtauCT(nsub,:,:)=COVtauemp; % COVtauemp for each subject
end
FCemp=squeeze(mean(FCCT)); %mean FC matrix across all subjects
COVtauemp=squeeze(mean(COVtauCT)); % mean covariance matrix acros subjects (IR?)
Cnew=C; %Initialize effective connectivity matrix with anatomical connectivity
olderror=100000; % convergence, start with high error
for iter=1:5000
    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_c,sigma); % simulates hopf dynamics with existing connectivity, learning rule sigma, and f_diff freq
    COVtausim=expm((Tau*TR)*A)*COVsimtotal; % Calculates simulated cov at time delay Tau by exponentiating matrix A with product of Tau, TR, covimtotal
    COVtausim=COVtausim(1:N,1:N); %extracts submatrix of covtausim corresponding to index regions
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j)); %calculates significance ratio (reciprocal of product of standard deviations i-th and j-th)
        end
    end
    COVtausim=COVtausim.*sigratiosim; %scale simulated time-delayed covariance matrix by significance ratio
    errorFC(iter)=mean(mean((FCemp-FCsim).^2)); %calculate error in effective matrix vs FC
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2)); %calculate error in IR matrix vs FC

    if mod(iter,100)<0.1 %every 100 iterations
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2)); % error = diff b/w empirical and simulated squared (sum both matrices)
        if  (olderror-errornow)/errornow<0.001 % if relative diff b/w old and new error is <0.001, stop optimizing
            break;
        end
        if  olderror<errornow %if old error is less than new, break
            break;
        end
        olderror=errornow;
    end

    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || j==N-i+1) % if C is positive (actually existing connections) or i/j correspond to same element in reverse diagonal
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j)); %update connection at i,j based on FCemp and Covtauemp error, scaled by learning rule
                if Cnew(i,j)<0 %if you somehow end up with negative connection
                    Cnew(i,j)=0; %no connection
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC; %normalize updated connectivity by dividing all elements by maximum value in matrix and scaling to maximum value. keeps all value within some desired range
end
CeffgroupCT=Cnew; %final optimized connectivity for control group (avg'd over all subjects)

%% Individual Control
for nsub=1:NSUB
    nsub
    ts=subject{nsub}.dbs80ts;  % fMRI CT
    %ts=ts';
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCCT(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=CeffgroupCT; %sets eff. conn to group -> closer to indiv. than a totally fresh one
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_c,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.001
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N  %% learning
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    % This part is different for the individuals
    Ceff=Cnew; %Effective = Cnew
    CeffCT(nsub,:,:)=Ceff; %Subject by subject Ceff
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff_c,sigma); %extract FCsim and COVsim
    fittFC_CT(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag)); %Correlation between empirical FC and simulated FC
    COVtausim=expm((Tau*TR)*A)*COVsimtotal; %matrix exponential of Jacobian scaled by covariance matrix (all regions incl.)
    COVtausim=COVtausim(1:N,1:N); %indexed
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim; %scale by sig ratio
    fittCVtau_CT(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag)); %correlation cov matrix
end

%%%%%%%%
%% AYA_SES02 (SZ) AYA
%%%%%%%%
load aya_ses02_hopf.mat;
load aya_ses02.mat;
f_diff_s=f_diff;
NSUB=24;

%% Group
for nsub=1:NSUB
    ts=subject{nsub}.dbs80ts;
    %ts=ts';
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCSZ(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    COVtauSZ(nsub,:,:)=COVtauemp;
end
FCemp=squeeze(mean(FCSZ));
COVtauemp=squeeze(mean(COVtauSZ));
Cnew=C;
olderror=100000;
for iter=1:5000

    % Linear Hopf FC
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_s,sigma);
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    errorFC(iter)=mean(mean((FCemp-FCsim).^2));
    errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

    if mod(iter,100)<0.1
        errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
        if  (olderror-errornow)/errornow<0.001
            break;
        end
        if  olderror<errornow
            break;
        end
        olderror=errornow;
    end

    for i=1:N  %% learning
        for j=1:N
            if (C(i,j)>0 || j==N-i+1)
                Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                    +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                if Cnew(i,j)<0
                    Cnew(i,j)=0;
                end
            end
        end
    end
    Cnew = Cnew/max(max(Cnew))*maxC;
end
CeffgroupSZ=Cnew;

%% Individual
for nsub=1:NSUB
    nsub
    ts=subject{nsub}.dbs80ts;
    %ts=ts';
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
    end
    % FC(0)
    ts2=ts(indexN,10:end-10);
    Tm=size(ts2,2);
    FCemp=corrcoef(ts2');
    FCSZ(nsub,:,:)=FCemp;
    COVemp=cov(ts2');
    % COV(tau)
    tst=ts2';
    for i=1:N
        for j=1:N
            sigratio(i,j)=1/sqrt(COVemp(i,i))/sqrt(COVemp(j,j));
            [clag lags] = xcov(tst(:,i),tst(:,j),Tau);
            indx=find(lags==Tau);
            COVtauemp(i,j)=clag(indx)/size(tst,1);
        end
    end
    COVtauemp=COVtauemp.*sigratio;
    Cnew=CeffgroupSZ;
    olderror=100000;
    for iter=1:5000
        % Linear Hopf FC
        [FCsim,COVsim,COVsimtotal,A]=hopf_int(Cnew,f_diff_s,sigma);
        COVtausim=expm((Tau*TR)*A)*COVsimtotal;
        COVtausim=COVtausim(1:N,1:N);
        for i=1:N
            for j=1:N
                sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
            end
        end
        COVtausim=COVtausim.*sigratiosim;
        errorFC(iter)=mean(mean((FCemp-FCsim).^2));
        errorCOVtau(iter)=mean(mean((COVtauemp-COVtausim).^2));

        if mod(iter,100)<0.1
            errornow=mean(mean((FCemp-FCsim).^2))+mean(mean((COVtauemp-COVtausim).^2));
            if  (olderror-errornow)/errornow<0.001
                break;
            end
            if  olderror<errornow
                break;
            end
            olderror=errornow;
        end

        for i=1:N  %% learning
            for j=1:N
                if (C(i,j)>0 || j==N-i+1)
                    Cnew(i,j)=Cnew(i,j)+epsFC*(FCemp(i,j)-FCsim(i,j)) ...
                        +epsFCtau*(COVtauemp(i,j)-COVtausim(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                end
            end
        end
        Cnew = Cnew/max(max(Cnew))*maxC;
    end
    Ceff=Cnew;
    CeffSZ(nsub,:,:)=Ceff;
    [FCsim,COVsim,COVsimtotal,A]=hopf_int(Ceff,f_diff_s,sigma);
    fittFC_SZ(nsub)=corr2(FCemp(Isubdiag),FCsim(Isubdiag));
    COVtausim=expm((Tau*TR)*A)*COVsimtotal;
    COVtausim=COVtausim(1:N,1:N);
    for i=1:N
        for j=1:N
            sigratiosim(i,j)=1/sqrt(COVsim(i,i))/sqrt(COVsim(j,j));
        end
    end
    COVtausim=COVtausim.*sigratiosim;
    fittCVtau_SZ(nsub)=corr2(COVtauemp(Isubdiag),COVtausim(Isubdiag));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hierarchy Ceff (TROPHIC COHERENCE)

%% Control
for nsub=1:NSUB
    Ceff=squeeze(CeffCT(nsub,:,:));
    Ceff(find(Ceff<0.015))=0;
    A=Ceff';
    d=sum(A)'; %Gin, sum of each column
    delta=sum(A,2); %out, sum of each row
    u=d+delta; % incoming+outgoing
    v=d-delta; % incoming-outgoing
    Lambda=diag(u)-A-A'; %diagonal - outgoing - incoming
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v); % Solve the linear system Lambda * gamma = v for gamma
    gamma=gamma-min(gamma);  % Shift gamma values to make them all non-negative
    hierarchicallevelsCT(nsub,:)=gamma';
    H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2; % Create a matrix H representing squared differences between gamma values
    F0=sum(sum((A.*H)))/sum(sum(A)); % Calculate F0, a measure of trophic coherence, for the current subject
    trophiccoherenceCT(nsub)=1-F0;
end

%% Schizos
for nsub=1:NSUB
    Ceff=squeeze(CeffSZ(nsub,:,:));
    Ceff(find(Ceff<0.015))=0;
    A=Ceff';
    d=sum(A)';
    delta=sum(A,2);
    u=d+delta;
    v=d-delta;

    Lambda=diag(u)-A-A';
    Lambda(1,1)=0;
    gamma=linsolve(Lambda,v);
    gamma=gamma-min(gamma);
    hierarchicallevelsSZ(nsub,:)=gamma';
    H=(meshgrid(gamma)-meshgrid(gamma)'-1).^2;
    F0=sum(sum((A.*H)))/sum(sum(A));
    trophiccoherenceSZ(nsub)=1-F0;
end

% incoming outgoing G vals
for nsub=1:NSUB
    Ceff=squeeze(CeffCT(nsub,:,:));
    Ceff(find(Ceff<0.015))=0;
    A=Ceff';
    G_inCT(nsub,:)=sum(A)';
    G_outCT(nsub,:)=sum(A,2);

    Ceff=squeeze(CeffSZ(nsub,:,:));
    Ceff(find(Ceff<0.015))=0;
    A=Ceff';
    G_inSZ(nsub,:)=sum(A)';
    G_outSZ(nsub,:)=sum(A,2);
end
G_totalCT=G_inCT+G_outCT;
G_totalSZ=G_inSZ+G_outSZ;

ranksum(trophiccoherenceCT, trophiccoherenceSZ)
boxplot([trophiccoherenceCT' trophiccoherenceSZ'])

%{
save results_Ceff_aya.mat CeffCT CeffSZ  ...
    hierarchicallevelsCT hierarchicallevelsSZ trophiccoherenceSZ ...
    trophiccoherenceCT fittFC_CT fittFC_SZ fittCVtau_CT fittCVtau_SZ...
G_inCT G_inSZ G_outCT G_outSZ G_totalCT G_totalSZ;
%}
