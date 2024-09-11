
% this is the model-free script

clear all;
path0=pwd;
path2=[ '/Users/myco/Documents/data/2023-msc_thesis/ayahuasca'];
addpath(path2);
cd(path2);

%indexN=[1:31 50:80];
indexN=[1:80];

N=80;
Tau=3;
NSUB=16;

% Parameters
TR=2.0;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass -> for movies it is 0.08; however for my HADES work I used 0.09 so will continue with that
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Trem = 10;               % number of frames removed from front and back of timeseries


load aya_ses01.mat;
FowRev_sober=zeros(1,NSUB);
FowRev_ayahuasca=zeros(1,NSUB);

for sub=1:NSUB  % over subjects - 24
    ts2=subject{sub}.dbs80ts;
    ts=ts2(indexN,:); % remove regions if e.g. removing subcortical
    clear signal_filt;

    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)),"omitnan");    %standard filtering
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    end
    ts=signal_filt(:,Trem:end-Trem); % removal of the frames at the very beginning and end of each time series

    %ts=ts(:,Trem:end-Trem);
    Tm=size(ts,2);

    FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');       %% Core...FC tau forward
    FCtr=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)'); %% FC tau reversal
    Itauf=-0.5*log(1-FCtf.*FCtf);  %% Mutual information...
    Itaur=-0.5*log(1-FCtr.*FCtr);

    Reference=((Itauf(:)-Itaur(:)).^2)';
    index=find(Reference>quantile(Reference,0.0));
    FowRev_sober(sub)=mean(Reference(index),"omitnan");
    Tenet2_R(sub,:,:)=abs(Itauf-Itaur);
end
for nsub=1:NSUB
    tmp=squeeze(Tenet2_R(nsub,:,:)); %placebo
    m  = triu(true(size(tmp)));
    n = tril(true(size(tmp)));
    tmp=tmp(m)+tmp(n); % exclude diagonal
    hierarchy_sober(nsub)=std(tmp);
end
load aya_ses02.mat;

%NSUB=33;
for sub=1:NSUB
    ts2=subject{sub}.dbs80ts;
    ts=ts2(indexN,:);
    clear signal_filt;
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    end
    ts=signal_filt(:,Trem:end-Trem);
    Tm=size(ts,2);
    FCtf=corr(ts(:,1:Tm-Tau)',ts(:,1+Tau:Tm)');
    FCtr=corr(ts(:,Tm:-1:Tau+1)',ts(:,Tm-Tau:-1:1)');
    Itauf=-0.5*log(1-FCtf.*FCtf);
    Itaur=-0.5*log(1-FCtr.*FCtr);
    Reference=((Itauf(:)-Itaur(:)).^2)';
    index=find(Reference>quantile(Reference,0.0));
    FowRev_ayahuasca(sub)=nanmean(Reference(index));
    Tenet2_A(sub,:,:)=abs(Itauf-Itaur);
end
% Hierarchy
for nsub=1:NSUB
    tmp2=squeeze(Tenet2_A(nsub,:,:)); %drug
    m  = triu(true(size(tmp2)));
    n = tril(true(size(tmp2)));
    tmp2=tmp2(m)+tmp2(n);
    hierarchy_ayahuasca(nsub)=std(tmp2);
end

THRLOW=0;     % was 0
THRHIGH=100;   % was 100
FowRev_sober=rmoutliers(FowRev_sober,'percentiles',[THRLOW THRHIGH]);
FowRev_ayahuasca=rmoutliers(FowRev_ayahuasca,'percentiles',[THRLOW THRHIGH]);



%figure(1);
%boxplot([FowRev_sober' FowRev_ayahuasca']);
%save results_dmt_modelfree.mat FowRev_sober FowRev_ayahuasca Tenet2_R Tenet2_A hierarchy_sober hierarchy_ayahuasca;
