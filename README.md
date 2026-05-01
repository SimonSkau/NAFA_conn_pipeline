# NAFA_conn_pipeline
fNIRS pipleine for connectitvty analysis in the NAFA study

% clc
% close all
% clear all

for  sesID = % number of sessions1;
for subID = % study ID % 
  
    snirf = SnirfLoad([ sprintf('NAFA_0%02d',subID) '_rest_' num2str(sesID) '.snirf']);

[ sprintf('NAFA_0%02d',subID) '_rest_' num2str(sesID)]

procResult.mlActMan{1} = ones(size(snirf.GetMeasList,1),1);
procResult.tInc{1}= ones(size(snirf.Get_t,1),1);


mlActAuto = hmrR_PruneChannels(snirf.data, snirf.probe ,procResult.mlActMan , procResult.tInc, [1e-5 1e+3], 1   , [0.0 40.0]);

nSV = 0;
dod = hmrR_Intensity2OD( snirf.data );
[dod, svs, nSV] = hmrR_PCAFilter( dod, mlActAuto, procResult.tInc, nSV );


dod = hmrR_BandpassFilt(dod, 0.03, 0.5);

tMotion = 0.5;% Units of seconds.    Typical value ranges from 0.1 to 0.5.
tMask = 1; %Units of seconds. Typical value ranges from 0.5 to 1
STDEVthresh = 30;
AMPthresh =  1;
procResult.tIncAuto = hmrR_MotionArtifact(snirf.data, snirf.probe, procResult.mlActMan, mlActAuto, procResult.tInc, tMotion, tMask, STDEVthresh, AMPthresh);

[tInc,tIncCh] = hmrR_MotionArtifactByChannel(dod, snirf.probe, procResult.mlActMan, mlActAuto, procResult.tInc, tMotion, tMask, STDEVthresh, AMPthresh);

p = 0.99;
FrameSize_sec = 10;
dod = hmrR_MotionCorrectSpline(dod, procResult.mlActMan, tIncCh, p);

tRange =[-0 0];
[procResult.stim, tRange] = hmrR_StimRejection(dod, snirf.stim, procResult.tIncAuto, procResult.tInc, tRange);

dc = hmrR_OD2Conc( dod, snirf.probe, [6.0 6.0] );


 trange = [-5 300]; 

[data_avg, data_std, nTrials, data_sum2, yTrials] = hmrR_BlockAvg( dc, procResult.stim, trange );


result = data_avg.GetDataTimeSeries('reshape'); %block


coi = sort([1 2 3 4 9 10 11 38 30 31 32 34 13 14 15 24 33 35 46 47 23 ...
25 27 28 5 6 16 17 19 29 12 39 40 41 42 44 18 20 21 45 26 43 50 51]); %channel of interest 
result(result==0)=NaN;

oxy = squeeze(result(:,1,coi));
deoxy = squeeze(result(:,2,coi));
tot = squeeze(result(:,3,coi));
 
all_oxy_data{subID,sesID} = oxy;

[r,p]=corr(oxy);


for n = 1:44;
    for m = 1:44;
        if isnan(r(n,m));
            r(n,m)=0;
        else r(n,m);
        end
    end
end
    
r(r>0.9999)=0;
r(r<0)=0;


N = size(r,1);

% Ta bort diagonal
r_no_diag = r;
r_no_diag(1:N+1:end) = 0;

% Hitta alla unika l�nkar i �vre triangeln
upper_tri_vals = r_no_diag(triu(true(N), 1));

% Sortera i fallande ordning
sorted_vals = sort(upper_tri_vals, 'descend');

% Ber�kna antal l�nkar som motsvarar X%
num_links = length(sorted_vals);
num_keep = round(1 * num_links);

% Tr�skelv�rde = det minsta v�rdet bland de 20% starkaste l�nkarna
threshold = sorted_vals(num_keep);

% Nollst�ll alla v�rden under tr�skelv�rdet
r_thresh = zeros(size(r));
r_thresh(r >= threshold) = r(r >= threshold);

% Beh�ll symmetrin, s�tt diagonal till 0
r_thresh = triu(r_thresh,1);
r_thresh = r_thresh + r_thresh';

% Diagonalen �r noll
r_thresh(1:N+1:end) = 0;

save_r{subID,sesID} = r_thresh;

C_r = clustering_coef_wu(r_thresh);
cavg_r{subID,sesID} = nanmean(C_r);

RC_r = rich_club_wu(r_thresh);
Rich_r{subID,sesID} = RC_r;

Eg_r =efficiency_wei(r_thresh);
Eavg_r{subID,sesID} = nanmean(Eg_r);


[Ci Q] = modularity_und(r_thresh);
Modularity{subID,sesID} = Q;


% figure;
%  imagesc(r_thresh); 
%  axis square


   end
end



