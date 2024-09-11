% This script computes the average z-scored FC matrix (FCemp) and average dominant frequency
% (f_diff) across subjects through power spectrum analysis
% run this first

clear all;
path2data='/Users/myco/human_brains/dmt_update/';
addpath(path2data);
cd(path2data);

load pcb_post.mat; %load your data

N=80;
indexN=1:N; 

NSUB=16;

% Parameters of the data
TR=2.0;  % Repetition Time (seconds)
% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency (max freq that can be accuratley represented in discrete Fourier transform
flp = 0.008;                    % lowpass frequency of filter (Hz)
fhi = 0.08;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

for sub=1:NSUB
    sub    
    clear signal_filt Power_Areas;
    ts=subject{sub}.dbs80ts;
   % ts=ts';
    ts=ts(indexN,:);

    % bandpass filter
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-mean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    end
    signal_filt=signal_filt(:,10:end-10); % Remove first and last 10 frames
   
    [Ns, Tmaxred]=size(signal_filt); %Ns = nodes, Tmaxred = # timepoints after removing first and last 10 frames
    TT=Tmaxred; %timepoints after filtering
    Ts = TT*TR; % length of timeseries in seconds
    freq = (0:TT/2-1)/Ts; % 0Hz to Nyquist frequency
    nfreqs=length(freq);
    for seed=1:N
        pw = abs(fft(zscore(signal_filt(seed,:)))); % "power weight" power decomposition (FFT) of z-score transformed filtered signal
        PowSpect = pw(1:floor(TT/2)).^2/(TT/TR); %take first half of power spectrum, square, normalize by length of timeseries data
        Power_Areas=gaussfilt(freq,PowSpect,0.005); %gaussian filter power spectrum
        [~,index]=max(Power_Areas); % index of timepoint w highest power
        index=squeeze(index);
        f_diff_sub(sub,seed)=freq(index); % takes highest power frequency
    end
    
    ts=zscore(ts,[],2); %z-score normalize timeseries
    FCemp(sub,:,:)=corrcoef(ts'); %Functional connectivity of z-score normalized
end

f_diff = mean(f_diff_sub,1); %mean highest-power frequency across subjects for each region

FCemp=squeeze(nanmean(FCemp)); % mean z-score normalized FC matrix

save pcb_post_hopf.mat FCemp f_diff;

