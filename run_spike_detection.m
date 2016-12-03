
clear all;
close all;

% ADD IEEG PORTAL AND UTILITIES TO PATH
addpath(genpath('portal-matlab-tools/Libraries/ieeg-matlab-1.13.2'))
addpath(genpath('portal-matlab-tools/Utilities'))
addpath('portal-matlab-tools/Analysis')

userID= 'hoameng';
pwdFile = 'hoa_ieeglogin.bin';

session = IEEGSession('I022_P001_D01',userID,pwdFile);
subj = 1:19;
subjsess = 1:6;
for i = subj
    for j = subjsess
        dataName = sprintf('I022_P0%.2d_D%.2d',i,j);
        try
            session.openDataSet(dataName);
        catch
            fprintf('Unable to load %s\n',dataName)
        end
    end
end
subj = 1:54;
subjsess = 1:15;
for i = subj
    for j = subjsess
        dataName = sprintf('I027_P0%.2d_D%.2d',i,j);
        try
            session.openDataSet(dataName);
        catch
            fprintf('Unable to load %s\n',dataName)
        end
    end
end

datasetNames = cell(numel(session.data),1);
for i = 1:numel(session.data)
    datasetNames{i} = session.data(i).snapName;
end

%% remove channels by name
channelsToSkip={'EKG','ECG','RESP','EMG','OSAT','EVENT','X','CHIN','CHEST','SNORE','FLOW','FPZ','ABD'}; % exclude these channels
channelIdxs = cell(numel(session.data),1);
toRemove = cell(numel(session.data),1);
for i = 1:numel(session.data)
    chanLabels = {session.data(i).rawChannels.label}; % get channel labels for dataset
    channelIdxs{i} = 1:numel(chanLabels);
    removeChans = cellfun(@(x)cell2mat(regexpi(x,channelsToSkip)),chanLabels,'UniformOutput',0); %find matches
    toRemove{i} = find(cellfun(@(x)(~isempty(x)),removeChans)~=0);
    channelIdxs{i}(toRemove{i}) = []; %remove channels
    chanLabels(toRemove{i}) = [];
end

params.datasetID = datasetNames;
params.IEEGid = userID;
params.IEEGpwd = pwdFile;
%% remove artifact channels
thres = 3
[ch_keep, ch_remove, slopes] = identifyArtifactChannels_par(datasetNames,params,thres);

%% removing manual channels
a = cell2mat(cellfun(@(x)strcmp(x,'I027_P004_D06'),datasetNames,'UniformOutput',0));
ch_keep{a==1}(ismember(ch_keep{a==1},[48 49])) = [];

%merge channels to keep 
ch_keep_final = cellfun(@(x,y)(intersect(x,y)),ch_keep,channelIdxs,'UniformOutput',0);


%% DETECT SPIKES 
warning('off','all')
%find all spikes
parfor i =1:numel(datasetNames)
    try
        fprintf('Detecting spikes in : %s -\n',datasetNames{i});
        if (exist(fullfile(cd, sprintf('%s_spiketimes_ja.mat',datasetNames{i})), 'file') == 0)
            tmpsession = IEEGSession(datasetNames{i},userID,pwdFile);
            [spikeTimes, spikeChannels,DE] = spike_ja_wrapper(tmpsession.data,ch_keep_final{i});
            parsave(sprintf('%s_spiketimes_ja.mat',tmpsession.data.snapName),spikeTimes,spikeChannels,DE);
        else
            fprintf(' Mat found.\n')
        end
    catch ME
        fprintf('Error in %s, dataset %d, %s\n, ',datasetNames{i},i,ME.identifier);
    end
end

%Clean and upload
overwriteFlag = 0;
for i =1:numel(datasetNames)
    load(sprintf('%s_spiketimes_ja.mat',datasetNames{i}))
    spikeTimes2 = spikeTimes(DE.pdf<(0.05),:);
    spikeChannels2 = spikeChannels(:,DE.pdf<(0.05));
    
    %spatial integration
    [eventTimes, eventChannels] = spatialIntegration(spikeTimes2,num2cell(spikeChannels2),[50 150]);
    
    %remove spikes outside of task
    eventAnnotationName = 'pyFR_Events';
    [annots, timesUSec, ch] = getAnnotations(session.data(i),eventAnnotationName);
    eventChannels(eventTimes(:,2) < timesUSec(1,1) | eventTimes(:,1)>timesUSec(end,2)) = [];
    eventTimes(eventTimes(:,2)<timesUSec(1,1) | eventTimes(:,1)>timesUSec(end,2),:) = [];
    
    [a] = cellfun(@(x)numel(x),eventChannels);
    uploadAnnotations(session.data(i),'spike_ja_spatial',eventTimes(a>1,:),eventChannels(a>1),'spike',overwriteFlag);
end

