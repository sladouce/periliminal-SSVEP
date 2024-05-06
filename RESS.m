
clear all
close all
eeglab;
cd 'D:\Spatialattention_study2\outicrej\'
directory = 'D:\Spatialattention_study2\outicrej\';
SUBJ = {'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14','P15','P16','P17','P18','P19','P20','P21','P22','P23','P24'};%

COND = {'_subliminal','_periliminal','_control'};
stims = [13,15]; %target frequencies
stimsl = {'_left', '_right'};
zeropad = 4;%factor to zeropad or nfft, 1 corresponds to 1 hz frequency resolution
sr = 500;%sampling rate
compn = 1;%which component to plot (1=highest value)
allsub_RESS_spec = zeros(length(SUBJ),length(COND),length(stims),sr*zeropad);
allsub_RESS_SNR = zeros(length(SUBJ),length(COND),length(stims),sr*zeropad);

for s = 1:length(SUBJ)
    for c = 1:length(COND)
        for stim = 1:length(stims)

        EEG = pop_loadset([SUBJ{s},COND{c},stimsl{stim},'.set']);
            ssvepfrex = stims(stim);
            EEG.data = double(EEG.data);

% to whiten the data (better spatial accuracy)
dozca=1;

% specify RESS parameters
neig      = 1;  % distance to frequency neighbors in Hz 1
fwhm_targ = 1; % FWHM in Hz for target  1
fwhm_neig = 1;  % FWHM in Hz for neighbors  1

% shrinkage proportion for regularization
shr = .01;

% time window for ASSR
tidx = dsearchn(EEG.times',[0 6000]');

% number of time points in filter
pnts = (diff(tidx)+1)*EEG.trials;

% FFT param
nfft = EEG.srate*zeropad;
hz   = linspace(0,EEG.srate,nfft);

%%
% %%zca whitening
% 
tmpdat = reshape( EEG.data,EEG.nbchan,[] );
tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));

covmat = tmpdat*tmpdat' / (length(tmpdat)-1);

%covmat =cov(EEG.data')
[V,D] = eig(covmat);

% ZCA (note that yz*yz' = diag and almost lI)
yz = V*D^(-1/2)*V'*reshape(EEG.data,EEG.nbchan,[]);

% replace data with whitened data
if dozca
    EEG.data = reshape(yz,size(EEG.data));
end

%%
% now for RESS

for fi=1:length(ssvepfrex)
    clear comptsalltrials
   % S covariance matrix
    % filter
    data = filterFGx(EEG.data,EEG.srate,ssvepfrex(fi),fwhm_targ);
    % extract and mean-center data
    data = data(:,tidx(1):tidx(2),:);
    data = reshape( data ,EEG.nbchan,[] );
    data = bsxfun(@minus,data,mean(data,2));
    % covariance
    covS = (data*data') / (pnts-1);
 
%%%%%%%%%%%%%%%%%%%%%%%    
    
    % R covariance matrices    
    % lower R/2
   data = filterFGx(EEG.data,EEG.srate,ssvepfrex(fi)-neig,fwhm_neig);    
    % extract and mean-center data
    data = data(:,tidx(1):tidx(2),:);
    data = reshape( data ,EEG.nbchan,[] );
    data = bsxfun(@minus,data,mean(data,2));   
    % covariance
    covRl = (data*data') / (pnts-1);
  
    % filter
    data = filterFGx(EEG.data,EEG.srate,ssvepfrex(fi)+neig,fwhm_neig);
    % extract and mean-center data
    data = data(:,tidx(1):tidx(2),:);
    data = reshape( data ,EEG.nbchan,[] );
    data = bsxfun(@minus,data,mean(data,2));
    %covariance
    covRu = (data*data') / (pnts-1);
    
    %full R matrix is average of lower/upper
    covR = (covRl+covRu)/2;
    
    %GED with optional shrinkage
    % apply shrinkage to covR
    covR = (1-shr)*covR + shr*mean(eig(covR))*eye(size(covR));
    
    %GED and sort components
    [evecs,evals] = eig(covS,covR);
    [evals,sidx]  = sort(diag(evals),'descend');
    evecs = evecs(:,sidx);
    
    % compute filter forward model and flip sign
    map = covS*evecs(:,compn);
    [~,maxchan] = max(abs(map));
    map = map*sign(map(maxchan));

  % back projecting EEG data to component space
    compts_filtered = evecs(:,compn)'*reshape(EEG.data,EEG.nbchan,[]);
   
    compts_filtered = filterFGx(compts_filtered,EEG.srate,ssvepfrex(fi),1);
    compts_filtered=reshape(compts_filtered,EEG.pnts,EEG.trials);
   
    % psd over component with hilbert transform (better for band passed signal)
    compts1 = mean(abs(hilbert(compts_filtered(:,1:EEG.trials))).^2,2);
    compts2 = [];
    for ntrial = 1:EEG.trials
      compts2(:, ntrial) = abs(hilbert(compts_filtered(:,ntrial))).^2;
    end
        

    % compute component time series
    compts = evecs(:,compn)'*reshape(EEG.data,EEG.nbchan,[]);
    compts = reshape(compts,EEG.pnts,EEG.trials);
    % power spectrum averaged over trials
    powr = mean(abs(fft(compts(tidx(1):tidx(2),:),nfft,1)/EEG.pnts).^2 ,2);

    % SNR spectrum
    skipbins =  5; % .5 Hz 5
    numbins  = 2+skipbins; % 2 Hz 20
    snr = zeros(size(powr));
    
    for hzi=numbins+1:length(hz)-numbins-1
        % SNR over all time points and conditions
        numer = powr(hzi);
        denom = mean(powr([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snr(hzi) = numer./denom;
    end   
    
%saving in variable for averaging
allsub_RESS_spec(s,c,stim,1:zeropad*EEG.srate) = powr;
allsub_RESS_SNR(s,c,stim,1:zeropad*EEG.srate) = snr;
allsub_RESS_time(s,c,stim,1:size(EEG.data,2)) = compts1;
allsub_RESS_time_st(s,c,stim,1:size(EEG.data,2),1:EEG.trials) = compts2;
allsub_RESS_time_st_cong(s,c,stim,1:size(EEG.data,2),1:EEG.trials) = compts2;
allsub_RESS_map(s,c,stim,1:size(EEG.data,1)) = map;
allsub_RESS_eigval(s,c,stim,1:size(EEG.data,1)) = evals;  
    
   end

        end

    end
end

avg_RESS_SNR = permute(mean(allsub_RESS_SNR,1),[3 4 1 2]);
avg_RESS_SPEC = permute(mean(allsub_RESS_spec,1),[3 4 1 2]);
avg_RESS_map = permute(mean(allsub_RESS_map,1),[3 4 1 2]);
avg_RESS_time = permute(mean(allsub_RESS_time,1),[3 4 1 2]);
avg_RESS_eigval = permute(mean(allsub_RESS_eigval,1),[3 4 1 2]);


%% computing RESS main component for the non-target stimulation

stims3 = [15,13]; %non target frequencies
stims4 = {'_left', '_right'};


allsub_RESS_spec_b = zeros(length(SUBJ),length(COND),length(stims3),sr*zeropad);
allsub_RESS_SNR_b = zeros(length(SUBJ),length(COND),length(stims3),sr*zeropad);


for s = 1:length(SUBJ)
    for c = 1:length(COND)
        for stim = 1:length(stims3)

        EEG = pop_loadset([SUBJ{s},COND{c},stims4{stim},'.set']);

            ssvepfrex = stims3(stim);
            EEG.data = double(EEG.data);

% to whiten the data (better spatial accuracy)
dozca=1;

% specify RESS parameters
neig      = 1;  % distance to frequency neighbors in Hz
fwhm_targ = 1; % FWHM in Hz for target  
fwhm_neig = 1;  % FWHM in Hz for neighbors  

% shrinkage proportion for regularization
shr = .01;

% time window for ASSR
tidx = dsearchn(EEG.times',[0 6000]');

% number of time points in filter
pnts = (diff(tidx)+1)*EEG.trials;

% FFT param
nfft = EEG.srate*zeropad;
hz   = linspace(0,EEG.srate,nfft);

%%
% %%zca whitening
% 
tmpdat = reshape( EEG.data,EEG.nbchan,[] );
tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));

covmat = tmpdat*tmpdat' / (length(tmpdat)-1);

%covmat =cov(EEG.data')
[V,D] = eig(covmat);

% ZCA (note that yz*yz' = diag and almost lI)
yz = V*D^(-1/2)*V'*reshape(EEG.data,EEG.nbchan,[]);

% replace data with whitened data
if dozca
    EEG.data = reshape(yz,size(EEG.data));
end

%%
% now for RESS

for fi=1:length(ssvepfrex)
    clear comptsalltrials
   % S covariance matrix
    % filter
    data = filterFGx(EEG.data,EEG.srate,ssvepfrex(fi),fwhm_targ);
    % extract and mean-center data
    data = data(:,tidx(1):tidx(2),:);
    data = reshape( data ,EEG.nbchan,[] );
    data = bsxfun(@minus,data,mean(data,2));
    % covariance
    covS = (data*data') / (pnts-1);
 
%%%%%%%%%%%%%%%%%%%%%%%    
    
    % R covariance matrices    
    % lower R/2
   data = filterFGx(EEG.data,EEG.srate,ssvepfrex(fi)-neig,fwhm_neig);    
    % extract and mean-center data
    data = data(:,tidx(1):tidx(2),:);
    data = reshape( data ,EEG.nbchan,[] );
    data = bsxfun(@minus,data,mean(data,2));   
    % covariance
    covRl = (data*data') / (pnts-1);
  
    % filter
    data = filterFGx(EEG.data,EEG.srate,ssvepfrex(fi)+neig,fwhm_neig);
    % extract and mean-center data
    data = data(:,tidx(1):tidx(2),:);
    data = reshape( data ,EEG.nbchan,[] );
    data = bsxfun(@minus,data,mean(data,2));
    %covariance
    covRu = (data*data') / (pnts-1);
    
    %full R matrix is average of lower/upper
    covR = (covRl+covRu)/2;
    
    %GED with optional shrinkage
    % apply shrinkage to covR
    covR = (1-shr)*covR + shr*mean(eig(covR))*eye(size(covR));
    
    %GED and sort components
    [evecs,evals] = eig(covS,covR);
    [evals,sidx]  = sort(diag(evals),'descend');
    evecs = evecs(:,sidx);
    
    % compute filter forward model and flip sign
    map = covS*evecs(:,compn);
    [~,maxchan] = max(abs(map));
    map = map*sign(map(maxchan));

  % back projecting EEG data to component space
    compts_filtered = evecs(:,compn)'*reshape(EEG.data,EEG.nbchan,[]);
   
    compts_filtered = filterFGx(compts_filtered,EEG.srate,ssvepfrex(fi),1);
    compts_filtered=reshape(compts_filtered,EEG.pnts,EEG.trials);
   
    % psd over component with hilbert transform (better for band passed signal)
    compts1 = mean(abs(hilbert(compts_filtered(:,1:EEG.trials))).^2,2);
    
    %retrieve values for single trials
    compts2 = [];
    for ntrial = 1:EEG.trials
      compts2(:, ntrial) = abs(hilbert(compts_filtered(:,ntrial))).^2;
    end
    % compute component time series
    compts = evecs(:,compn)'*reshape(EEG.data,EEG.nbchan,[]);
    compts = reshape(compts,EEG.pnts,EEG.trials);
    % power spectrum averaged over trials
    powr = mean(abs(fft(compts(tidx(1):tidx(2),:),nfft,1)/EEG.pnts).^2 ,2);

    % SNR spectrum
    skipbins =  5; % .5 Hz 5
    numbins  = 2+skipbins; % 2 Hz 20
    snr = zeros(size(powr));
    
    for hzi=numbins+1:length(hz)-numbins-1
        % SNR over all time points and conditions
        numer = powr(hzi);
        denom = mean(powr([hzi-numbins:hzi-skipbins hzi+skipbins:hzi+numbins]) );
        snr(hzi) = numer./denom;
    end   
    
%saving in variable for averaging
allsub_RESS_SNR_b(s,c,stim,1:zeropad*EEG.srate) = snr;
allsub_RESS_time_b(s,c,stim,1:size(EEG.data,2)) = compts1;
allsub_RESS_time_st_b(s,c,stim,1:size(EEG.data,2),1:EEG.trials) = compts2;
allsub_RESS_time_st_incong(s,c,stim,1:size(EEG.data,2),1:EEG.trials) = compts2;

allsub_RESS_map_b(s,c,stim,1:size(EEG.data,1)) = map;
    
   end

        end

    end
end

avg_RESS_SNR_b = permute(mean(allsub_RESS_SNR_b,1),[3 4 1 2]);
avg_RESS_map_b = permute(mean(allsub_RESS_map_b,1),[3 4 1 2]);
avg_RESS_time_b = permute(mean(allsub_RESS_time_b,1),[3 4 1 2]);


%% Normalization of the attentional effect by subtracting non-target RESS signal from the target RESS signal
allsub_RESS_SNR_norm = allsub_RESS_SNR - allsub_RESS_SNR_b;
allsub_RESS_time_norm = allsub_RESS_time - allsub_RESS_time_b;
allsub_RESS_map_norm = allsub_RESS_map - allsub_RESS_map_b;
avg_RESS_map_norm = avg_RESS_map - avg_RESS_map_b;
allsub_RESS_time_st_norm = allsub_RESS_time_st - allsub_RESS_time_st_b;



%% plots
figure;
    for c = 1:length(COND)
    subplot(1,length(COND),c)
        temp = permute(mean(allsub_RESS_SNR_norm(:,c,1,:),1),[4 1 2 3]);
        temp2 = permute(std(allsub_RESS_SNR_norm(:,c,1,:),1)/sqrt(size(allsub_RESS_SNR_norm(:,c,1,:),1)),[4 1 2 3]);

        temp3 = permute(mean(allsub_RESS_SNR_norm(:,c,2,:),1),[4 1 2 3]);
        temp4 = permute(std(allsub_RESS_SNR_norm(:,c,2,:),1)/sqrt(size(allsub_RESS_SNR_norm(:,c,1,:),1)),[4 1 2 3]);
        [h, p] = boundedline(hz, temp, temp2, '-b', 'alpha', hz,temp3, temp4,'-r','alpha');
        axis([5 35 -5 40]); ax = gca;
        ax.XLim = [5 35];
        ax.YLim = [-5 40];
        ax.XTick = [5:5:35];
        ax.YTick = [-5:5:40];
        box off; 
        %PlotAxisAtOrigin;PlotAxisAtTarget;PlotAxisAtStop;
        set(h, 'LineWidth', 2);
    end

thresholds = [3.29, 2.14, 2.41];    
    
figure;
    for c = 1:length(COND)
    subplot(1,length(COND),c)
        temp = permute(mean(allsub_RESS_time_norm(:,c,1,:),1),[4 1 2 3]);
        temp2 = permute(std(allsub_RESS_time_norm(:,c,1,:),1)/sqrt(size(allsub_RESS_time_norm(:,c,1,:),1)),[4 1 2 3]);

        temp3 = permute(mean(allsub_RESS_time_norm(:,c,2,:),1),[4 1 2 3]);
        temp4 = permute(std(allsub_RESS_time_norm(:,c,2,:),1)/sqrt(size(allsub_RESS_time_norm(:,c,1,:),1)),[4 1 2 3]);
        [h, p] = boundedline(EEG.times, temp, temp2, '-b', 'alpha', EEG.times,temp3, temp4,'-r','alpha');
        axis([-3000 7000 -10 30]); ax = gca;
        ax.XLim = [-3000 7000];
        ax.YLim = [-10 30];
        ax.XTick = [-3000:1000:7000];
        ax.YTick = [-10:10:30];
        box off; PlotAxisAtOrigin;PlotAxisAtTarget;PlotAxisAtStop;
        yl = yline(thresholds(c),'-','Threshold')
        yl.LabelHorizontalAlignment = 'left';
        set(h, 'LineWidth', 2);
    end


%plotting one subject
% subjn = 24;
%     
% figure;
%     for c = 1:length(COND)
%     subplot(1,length(COND),c)
%         temp = permute(mean(allsub_RESS_time_norm(subjn,c,1,:),1),[4 1 2 3]);
%         temp2 = permute(std(allsub_RESS_time_norm(subjn,c,1,:),1)/sqrt(size(allsub_RESS_time_norm(subjn,c,1,:),1)),[4 1 2 3]);
% 
%         temp3 = permute(mean(allsub_RESS_time_norm(subjn,c,2,:),1),[4 1 2 3]);
%         temp4 = permute(std(allsub_RESS_time_norm(subjn,c,2,:),1)/sqrt(size(allsub_RESS_time_norm(subjn,c,1,:),1)),[4 1 2 3]);
%         [h, p] = boundedline(EEG.times, temp, temp2, '-b', 'alpha', EEG.times,temp3, temp4,'-r','alpha');
%         axis([-2000 6900 -10 30]); ax = gca;
%         ax.XLim = [-2000 6900];
%         ax.YLim = [-10 30];
%         ax.XTick = [-2000:1000:6900];
%         ax.YTick = [-10:10:30];
%         box off; PlotAxisAtOrigin;PlotAxisAtTarget;PlotAxisAtStop;
%         set(h, 'LineWidth', 2);
%     end    
    
    
figure;
cm_viridis = viridis(200);
for c = 1:length(COND)
for stim = 1:length(stims)
subplot(length(COND),length(stims),stim+(length(stims)*c)-length(stims))
temp = avg_RESS_map_norm(stim,:,1,c)';
EEGplot = EEG;
EEGplot.icawinv = temp(:,1);%only first component

% some arguments for headplot function: 'maplimits'='maxmin' or 'absmax' / 'colormap'
% 'electrode3d','view'
sm_view = [32 29];

%jet colormap
pop_headplot(EEGplot, 0, 1,[stimsl{stim}, ' Hz'],[1  1], 'setup',{'D:\\ISAE-SUPAERO\\Teaching\\demoBCI\\Dependencies\\spline.spl' 'meshfile' 'mheadnew.mat' 'transform' [-0.35579 -6.3369 12.3705 0.053324 0.018746 -1.5526 1.0637 0.98772 0.93269] },'colorbar','off','maplimits',[0 0.1],'electrodes','off','view',[sm_view(1) sm_view(2)]);
end
end
title('Grand average RESS Scalp maps')


%%
% Extracting single trial values
% Subj x amplitude depth x direction x freqs x trials
test = mean(allsub_RESS_time_st_norm(1,1,1,501:2001,:),4);

fixations = squeeze(mean(allsub_RESS_time_st_norm(:,:,:,501:2001,:),4));
cues = squeeze(mean(allsub_RESS_time_st_norm(:,:,:,2279:3501,:),4));
targets = squeeze(mean(allsub_RESS_time_st_norm(:,:,:,3501:5001,:),4));

subfixleft = [];
subfixright = [];
subcueleft = [];
subcueright = [];
subtargleft = [];
subtargright = [];

threshfixleft = [];
threshfixright = [];
threshcueleft = [];
threshcueright = [];
threshtargleft = [];
threshtargright = [];

controlfixleft = [];
controlfixright = [];
controlcueleft = [];
controlcueright = [];
controltargleft = [];
controltargright = [];

for s = 1:length(SUBJ)
    
subfixleft = [subfixleft; squeeze(fixations(s,1,1,:))];
subfixright = [subfixright; squeeze(fixations(s,1,2,:))];
subcueleft = [subcueleft; squeeze(cues(s,1,1,:))];
subcueright = [subcueright; squeeze(cues(s,1,2,:))];
subtargleft = [subtargleft; squeeze(targets(s,1,1,:))];
subtargright = [subtargright; squeeze(targets(s,1,2,:))];

threshfixleft = [threshfixleft; squeeze(fixations(s,2,1,:))];
threshfixright = [threshfixright; squeeze(fixations(s,2,2,:))];
threshcueleft = [threshcueleft; squeeze(cues(s,2,1,:))];
threshcueright = [threshcueright; squeeze(cues(s,2,2,:))];
threshtargleft = [threshtargleft; squeeze(targets(s,2,1,:))];
threshtargright = [threshtargright; squeeze(targets(s,2,2,:))];

controlfixleft = [controlfixleft; squeeze(fixations(s,3,1,:))];
controlfixright = [controlfixright; squeeze(fixations(s,3,2,:))];
controlcueleft = [controlcueleft; squeeze(cues(s,3,1,:))];
controlcueright = [controlcueright; squeeze(cues(s,3,2,:))];
controltargleft = [controltargleft; squeeze(targets(s,3,1,:))];
controltargright = [controltargright; squeeze(targets(s,3,2,:))];
end







%% plot time course of the RESS SSVEP responses to target and non-target (non normalized)

    
figure;
    for c = 1:length(COND)
    subplot(1,length(COND),c)
        temp = permute(mean(allsub_RESS_time(:,c,1,:),1),[4 1 2 3]);
        temp2 = permute(std(allsub_RESS_time(:,c,1,:),1)/sqrt(size(allsub_RESS_time_norm(:,c,1,:),1)),[4 1 2 3]);

        temp3 = permute(mean(allsub_RESS_time_b(:,c,1,:),1),[4 1 2 3]);
        temp4 = permute(std(allsub_RESS_time_b(:,c,1,:),1)/sqrt(size(allsub_RESS_time_norm(:,c,1,:),1)),[4 1 2 3]);
        [h, p] = boundedline(EEG.times, temp, temp2, '-b', 'alpha', EEG.times,temp3, temp4,':b','alpha');
        axis([-3000 7000 -10 30]); ax = gca;
        ax.XLim = [-3000 7000];
        ax.YLim = [-10 30];
        ax.XTick = [-3000:1000:7000];
        ax.YTick = [-10:10:30];
        box off; PlotAxisAtOrigin;PlotAxisAtTarget;PlotAxisAtStop;
        set(h, 'LineWidth', 2);

        temp5 = permute(mean(allsub_RESS_time(:,c,2,:),1),[4 1 2 3]);
        temp6 = permute(std(allsub_RESS_time(:,c,2,:),1)/sqrt(size(allsub_RESS_time_norm(:,c,1,:),1)),[4 1 2 3]);
        temp7 = permute(mean(allsub_RESS_time_b(:,c,2,:),1),[4 1 2 3]);
        temp8 = permute(std(allsub_RESS_time_b(:,c,2,:),1)/sqrt(size(allsub_RESS_time_norm(:,c,1,:),1)),[4 1 2 3]);


        [h, p] = boundedline(EEG.times, temp5, temp6, '-r', 'alpha', EEG.times,temp7, temp8,':r','alpha');
        axis([-3000 7000 -10 30]); ax = gca;
        ax.XLim = [-3000 7000];
        ax.YLim = [-10 30];
        ax.XTick = [-3000:1000:7000];
        ax.YTick = [-10:10:30];
        box off; PlotAxisAtOrigin;PlotAxisAtTarget;PlotAxisAtStop;
        set(h, 'LineWidth', 2);
    end

