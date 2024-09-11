% iterated tau over model-free results

%% Four Groups
clear all;
path0=pwd; 
path2=[ '/Users/myco/human_brains/dmt_update/'];
addpath(path2);
cd(path2);

%indexN=[1:31 50:80];
indexN=[1:80];

N=80;
NSUB=16;

% Parameters
TR=2.0;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Trem = 10;               % number of frames removed from front and back of timeseries

FowRev_sober=zeros(1,NSUB);
FowRev_ayahuasca=zeros(1,NSUB);
for Tau=1:10
    load dmt_pre.mat;
        for sub=1:NSUB 
            ts2=subject{sub}.dbs80ts;
            ts=ts2(indexN,:); % removing subcortical so would skip for mine
            clear signal_filt;
            
            for seed=1:N
                ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));    %standard filtering
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
            FowRev_sober(Tau,sub)=nanmean(Reference(index));
            Tenet2_R(sub,:,:)=abs(Itauf-Itaur); 
        end

        load dmt_post.mat;
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
            FowRev_ayahuasca(Tau,sub)=nanmean(Reference(index));        
            Tenet2_A(sub,:,:)=abs(Itauf-Itaur);
        end

        load pcb_pre.mat;
        for sub=1:NSUB 
            ts2=subject{sub}.dbs80ts;
            ts=ts2(indexN,:); % removing subcortical so would skip for mine
            clear signal_filt;
            
            for seed=1:N
                ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));    %standard filtering
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
            FowRev_sober2(Tau,sub)=nanmean(Reference(index));
            Tenet2_R(sub,:,:)=abs(Itauf-Itaur); 
        end
    
        load pcb_post.mat;
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
            FowRev_ayahuasca2(Tau,sub)=nanmean(Reference(index));        
            Tenet2_A(sub,:,:)=abs(Itauf-Itaur);
        end

        THRLOW=0;     % was 0 
        THRHIGH=100;   % was 100
        FowRev_sober=rmoutliers(FowRev_sober,'percentiles',[THRLOW THRHIGH]);
        FowRev_ayahuasca=rmoutliers(FowRev_ayahuasca,'percentiles',[THRLOW THRHIGH]);
        FowRev_sober2=rmoutliers(FowRev_sober2,'percentiles',[THRLOW THRHIGH]);
        FowRev_ayahuasca2=rmoutliers(FowRev_ayahuasca2,'percentiles',[THRLOW THRHIGH]);
end

for Tau=1:10
    pre=[FowRev_sober'; FowRev_sober2'];
    post=[FowRev_ayahuasca'; FowRev_ayahuasca2'];
    group=categorical([repelem("DMT",NSUB,1); repelem("Placebo",NSUB,1)]);
    pre=pre(:,Tau);
    post=post(:,Tau);
    subj=repmat(1:NSUB,1,2)';
    X=table(pre,post,group,subj,'VariableNames',{'pre','post','group','subj'});
    X.group=categorical(X.group);
    mdl=fitlme(X,'post~group+pre+(1|subj)');
    p(Tau)=mdl.Coefficients.pValue(3)
end

% make figures
fname="Tau_DMT";
hfig=figure('Name','Autocorrelation time lag');

plot(1:10,p,'Color',[1 0.38 0.53],'LineWidth',1.5,'DisplayName','ACF');
hold on;
xline(find(p==min(p)),'--',{'Ideal','Time lag'});
xlabel('Time lags (\tau)');
ylabel('Significance');

make_pretty_print('Ideal Tau DMT',hfig,'/Users/myco/Desktop/')

%% PLOTS
clear all;
load results_dmt_modelfree.mat
% Sample data for demonstration
numSubjects = 15;

% Bin the pre-injection values into discrete categories (e.g., 5 bins)
numBins = 5;
[~,~,prePlaceboBin] = histcounts(FowRev_pcbpre', numBins);
[~,~,preDMTBin] = histcounts(FowRev_dmtpre', numBins);

% Compute mean post-injection values for each combination of treatment and bin
meanPostPlacebo = accumarray(prePlaceboBin, FowRev_pcbpost', [], @mean);
meanPostDMT = accumarray(preDMTBin, FowRev_dmtpost', [], @mean);

% Plotting
figure; % create a new figure
hold on; % allow multiple plots on the same figure

plot(meanPostPlacebo, '-o', 'DisplayName', 'Placebo');
plot(meanPostDMT, '-x', 'DisplayName', 'DMT');

% Formatting the plot
xticks(1:numBins);
xticklabels(cellstr(num2str((1:numBins)')));
xlabel('Binned Pre-Injection Trophic Coherence');
ylabel('Mean Post-Injection Trophic Coherence');
legend('Location', 'best');
title('Interaction Plot: Effect of Treatment Across Binned Pre-Values');
grid on; % adds a grid for better readability
hold off; % stop plotting on the same figure



%% 2 groups

clear all;
path0=pwd; 
path2=[ '/Users/myco/human_brains/dmt_update/'];
addpath(path2);
cd(path2);

%indexN=[1:31 50:80];
indexN=[1:80];

N=80;
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

FowRev_sober=zeros(1,NSUB);
FowRev_ayahuasca=zeros(1,NSUB);
for Tau=1:10
    load pcb_post.mat;
    %NSUB=28;
    
        for sub=1:NSUB 
            ts2=subject{sub}.dbs80ts;
            ts=ts2(indexN,:); % removing subcortical so would skip for mine
            clear signal_filt;
            
            for seed=1:N
                ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));    %standard filtering
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
            FowRev_sober(Tau,sub)=nanmean(Reference(index));
            Tenet2_R(sub,:,:)=abs(Itauf-Itaur); 
        end
        for nsub=1:NSUB
            tmp=squeeze(Tenet2_R(nsub,:,:)); %placebo
            m  = triu(true(size(tmp)));
            n = tril(true(size(tmp)));
            tmp=tmp(m)+tmp(n); % exclude diagonal
            hierarchy_sober(nsub)=std(tmp);
        end
    
        load dmt_post.mat;
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
            FowRev_ayahuasca(Tau,sub)=nanmean(Reference(index));        
            Tenet2_A(sub,:,:)=abs(Itauf-Itaur);
        end

        THRLOW=0;     % was 0 
        THRHIGH=100;   % was 100
        FowRev_sober=rmoutliers(FowRev_sober,'percentiles',[THRLOW THRHIGH]);
        FowRev_ayahuasca=rmoutliers(FowRev_ayahuasca,'percentiles',[THRLOW THRHIGH]);
end
        %figure(1);
        %violinplot([FowRev_sober' FowRev_ayahuasca']);
        p=zeros(10,1);
        h=zeros(10,1);
        for Tau=1:10
        [p(Tau),h(Tau)] = ranksum(FowRev_sober(Tau,:),FowRev_ayahuasca(Tau,:))
        end

% Calculate hierarchy

for nsub=1:NSUB
    tmp2=squeeze(Tenet2_A(nsub,:,:)); %drug
    m  = triu(true(size(tmp2)));
    n = tril(true(size(tmp2)));
    tmp2=tmp2(m)+tmp2(n);
    hierarchy_ayahuasca(nsub)=std(tmp2);
end

% make figures
fname="Tau_CO";
hfig=figure('Name','Ideal Tau');

plot(1:10,p,'Color',[1 0.38 0.53],'LineWidth',1.5,'DisplayName','ACF');
hold on;
xline(3,'--',{'Ideal','Time lag'});
xlabel('Time lags (\tau)');
ylabel('Significance');

make_pretty_print('Ideal Tau CO',hfig,'/Users/myco/Desktop/')