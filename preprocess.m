%preprocessing
%% loading eeglab
clear all;
close all;
eeglab; 

PATHIN = 'D:\Spatialattention_study\';
SUBJ = {'P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14','P15','P16','P17','P18','P19','P20','P21','P22','P23','P24'};%
COND = {'_subliminal','_periliminal','_control'};

%% loading data set

for s = 1:length(SUBJ)
    for c = 1:length(COND)
cd 'D:\Spatialattention_study\'
EEG = pop_loadset([SUBJ{s},COND{c},'.set']);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
%% loading channel loc

%EEG = pop_select( EEG, 'nochannel',{'ACC_X' 'ACC_Y' 'ACC_Z'});
EEG=pop_chanedit(EEG, 'lookup','D:\\Softwares\\eeglab_latest\\plugins\\dipfit\\standard_BESA\\standard-10-5-cap385.elp','load',{'D:\\ISAE-SUPAERO\\Dependencies\\live_amp.ced' 'filetype' 'autodetect'});
%% Band pass filtering 

%EEG.data = detrend(EEG.data);
EEG = pop_eegfiltnew(EEG,1,40,[]);
[EEGout,indelec]= pop_rejchan(EEG, 'elec',[1:32] ,'threshold',3,'norm','on','measure','spec','freqrange',[1 250] );
badelecs(s,c) = length(indelec);
EEG = pop_interp(EEG, indelec, 'spherical');

originalEEG = EEG;
EEG.nbinterp = length(indelec);
EEG.interp = indelec;

%average referencing
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:) = zeros(1, EEG.pnts);
EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, []);
EEG = pop_select( EEG,'nochannel',{'initialReference'});

% ICA
EEG = pop_runica(EEG, 'icatype', 'acsobiro', 'pca', EEG.nbchan-length(indelec));
EEG = eeg_checkset(EEG);

%ICLabel automatically tag components 
EEG = pop_iclabel(EEG, 'default');%%
BrainComp= find(EEG.etc.ic_classification.ICLabel.classifications(:,1) > 0.8);
EyeComp= find(EEG.etc.ic_classification.ICLabel.classifications(:,3) > 0.7);
OtherComp= find(EEG.etc.ic_classification.ICLabel.classifications(:,7) > 0.7);
MuscleComp= find(EEG.etc.ic_classification.ICLabel.classifications(:,2) > 0.7);
HeartComp= find(EEG.etc.ic_classification.ICLabel.classifications(:,4) > 0.7);
ElectrodeComp= find(EEG.etc.ic_classification.ICLabel.classifications(:,6) > 0.7);
LineNoise= find(EEG.etc.ic_classification.ICLabel.classifications(:,5) > 0.4);
badcomps = [EyeComp',MuscleComp',LineNoise',OtherComp',HeartComp',ElectrodeComp'];

EEG = pop_subcomp(EEG,[badcomps],0,0);%discard components classified as artifacts
EEG.badcomps = badcomps;

%% Epoching

    
    EEG1 = pop_epoch( EEG, {'cueleft'}, [-4 7], 'newname', strcat('wholeepochleft'), 'epochinfo', 'yes');
    EEG2 = pop_epoch( EEG, {'cueright'}, [-4 7], 'newname', strcat('wholeepochright'), 'epochinfo', 'yes');

cd 'D:\Spatialattention_study\outicrej\'

    pop_saveset(EEG1,strcat([SUBJ{s},COND{c},'_left']));
    pop_saveset(EEG2,strcat([SUBJ{s},COND{c},'_right']));

    end
end

