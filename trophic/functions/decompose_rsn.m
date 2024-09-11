% bobby.tromm@gmail.com
% Bobby Tromm 21/09/23
% Decompose trophic levels from parcellated data into canonical (Yeo)
% resting-state networks

function [yeo] = decompose_rsn(hl)

load('/Users/myco/computer_magic/nifty_utils/dbs2yeo/dk80toyeo7.mat'); %32-50 are subcortical structures, not part of RSNs!
av(32:49)=8;

hlyeobsl=zeros(size(hl,1),7);
hlyeobsl_tally=zeros(size(hl,1),7);
for i=1:80
    if i~=32:49
        if ~isnan(hl(:,i))
            hlyeobsl(:,av(i))=hlyeobsl(:,av(i))+hl(:,i);
            hlyeobsl_tally(:,av(i))=hlyeobsl_tally(:,av(i))+1;
        end
    end
end
yeo=hlyeobsl./hlyeobsl_tally;