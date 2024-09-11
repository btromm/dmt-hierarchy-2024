% run this first
% to calculate the optimal tau value
% note: requires econometrics toolbox

clear all;
path0=pwd;
data='/Users/myco/human_brains/ayahuasca';
results_path='/Users/myco/human_brains/INSIDEOUT_THESIS';
addpath(results_path);
addpath(data);
cd(data);

%indexN=[1:31 50:80]; % removes subcortical structures; if you want to include then indexN=[1:80]
indexN=[1:80];
N_Regions=length(indexN); % 80 if subcortical regions included
N=N_Regions;
NSUB=24; % number subjects
numConditions=2;

% Filtering Parameters
TR=1.4;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass -> for movies it is 0.08; however for my HADES work I used 0.09 so will continue with that
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
Trem = 10;               % number of frames removed from front and back of timeseries

tly=1;
for con=1:numConditions
    if con==1
        load aya_ses01.mat
    elseif con==2
        load aya_ses02.mat
    end

    for sub=1:NSUB
        ts2=subject{sub}.dbs80ts;
        ts=ts2(indexN,:); % removing subcortical so would skip for mine
        clear signal_filt;
    
        for seed=1:N
            ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));    %standard filtering
            signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        end
        subj=signal_filt(:,Trem:end-Trem); % removal of the frames at the very beginning and end of each time series
        for n=1:N_Regions
            ts=subj(n,:);
            
            [acf,lags]=autocorr(ts);
            acff(tly,:)=acf;
            lagg(tly,:)=lags;
            tly=tly+1;
            index=find(acf<0,1,'first');
            
            if ~exist('tau','var')
                tau(1,1)=lags(index-1);
                tau(1,2)=acf(index-1);
            else
                tau(end+1,1)=lags(index-1);
                tau(end,2)=acf(index-1);
            end
        end
    end
end

tau_median=median(tau)
tau_mean=mean(tau)

%% make figures
fname="Tau_Aya";
hfig=figure('Name','Autocorrelation time lag');
subplot(2,1,2);
histogram(tau(:,1),'FaceColor',[1 0.38 0.53],'FaceAlpha',0.5);
xlabel('Tau');
ylabel('Frequency');
title('Determination of $\tau$ for time-shifted correlation');
subplot(2,1,1);
plot(nanmean(lagg),nanmean(acff),'Color',[1 0.38 0.53],'LineWidth',1.5,'DisplayName','ACF');
hold on;
xline(tau_mean(1),'--',{'Sufficient','Decay'});
xlabel('Time lags $(\tau)$');
ylabel('Mean autocorrelation');

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',12) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fullfile(results_path,fname),'-dpng','-vector')
