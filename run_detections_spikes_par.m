
%% Establish IEEG Sessions
% Establish IEEG Sessions through the IEEGPortal. This will allow on demand
% data access
addpath(genpath('Z:\public\USERS\hoameng\Libraries\ieeg-matlab-1.13.2'));
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox\portalGit\Analysis'));
addpath(genpath('Z:\public\USERS\hoameng\Projects\p05-IEEGPortalToolbox\portalGit\Utilities'));

clear all;
userID= 'hoameng';
pwdFile = 'hoa_ieeglogin.bin';
session = IEEGSession('I022_P001_D01','hoameng','hoa_ieeglogin.bin');
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

[subjid durations] = getAllDurations(session.data);

%% remove channels
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
[ch_keep, ch_remove, slopes] = identifyArtifactChannels_par(datasetNames,params,3);

%% removing manual channels
a = cell2mat(cellfun(@(x)strcmp(x,'I027_P004_D06'),datasetNames,'UniformOutput',0));
ch_keep{a==1}(ismember(ch_keep{a==1},[48 49])) = [];

%merge channels to keep 
ch_keep_final = cellfun(@(x,y)(intersect(x,y)),ch_keep,channelIdxs,'UniformOutput',0);


%% DETECT SPIKES 
warning('off','all')
%find all spikes
for i =1:numel(datasetNames)
    try
        fprintf('Detecting spikes in : %s -\n',datasetNames{i});
        if (exist(fullfile(cd, sprintf('%s_spiketimes_ja2.mat',datasetNames{i})), 'file') == 0)
            tmpsession = IEEGSession(datasetNames{i},userID,pwdFile);
            [spikeTimes, spikeChannels,DE] = spike_ja_wrapper(tmpsession.data,ch_keep_final{i});
            parsave(sprintf('%s_spiketimes_ja2.mat',tmpsession.data.snapName),spikeTimes,spikeChannels,DE);
        else
            fprintf(' Mat found.\n')
        end
    catch ME
        fprintf('Error in %s, dataset %d, %s\n, ',datasetNames{i},i,ME.identifier);
    end
end

%UPLOAD all and bonf corrected
overwriteFlag = 0;
for i =1:numel(datasetNames)
    load(sprintf('%s_spiketimes_ja.mat',datasetNames{i}))
    spikeTimes2 = spikeTimes(DE.pdf<(0.05),:);
    spikeChannels2 = spikeChannels(:,DE.pdf<(0.05));
    %uploadAnnotations(session.data(i),'spike_ja',spikeTimes,num2cell(spikeChannels),'spike',overwriteFlag);
    [eventTimes, eventChannels] = spatialIntegration(spikeTimes2,num2cell(spikeChannels2),[50 150]);
    [a] = cellfun(@(x)numel(x),eventChannels);
    uploadAnnotations(session.data(i),'spike_ja_spatial',eventTimes(a>1,:),eventChannels(a>1),'spike',overwriteFlag);
    spikeTimes = spikeTimes(DE.pdf<(0.05/numel(DE.pdf)),:);
    spikeChannels = spikeChannels(:,DE.pdf<(0.05/numel(DE.pdf)));
    [eventTimes, eventChannels] = spatialIntegration(spikeTimes,num2cell(spikeChannels),[50 150]);
    [a] = cellfun(@(x)numel(x),eventChannels);
    uploadAnnotations(session.data(i),'spike_ja_bonf_spatial',eventTimes(a>1,:),eventChannels(a>1),'spike',overwriteFlag);
    uploadAnnotations(session.data(i),'spike_ja_bonf',eventTimes,eventChannels,'spike',overwriteFlag)
end

overwriteFlag = 0;
for i =59:numel(datasetNames)
    load(sprintf('%s_spiketimes_ja2.mat',datasetNames{i}))
    spikeTimes2 = spikeTimes(DE.pdf<(0.05),:);
    spikeChannels2 = spikeChannels(:,DE.pdf<(0.05));
    %uploadAnnotations(session.data(i),'spike_ja',spikeTimes,num2cell(spikeChannels),'spike',overwriteFlag);
    [eventTimes, eventChannels] = spatialIntegration(spikeTimes2,num2cell(spikeChannels2),[50 150]);
    
    %remove spikes outside of task
    eventAnnotationName = 'pyFR_Events';
    try
        [~, timesUSec, ch] = getAnnotations(session.data(i),eventAnnotationName);
        eventChannels(eventTimes(:,2) < timesUSec(1,1) | eventTimes(:,1)>timesUSec(end,2)) = [];
        eventTimes(eventTimes(:,2)<timesUSec(1,1) | eventTimes(:,1)>timesUSec(end,2),:) = [];
        [a] = cellfun(@(x)numel(x),eventChannels);

        uploadAnnotations(session.data(i),'spike_ja_trim',eventTimes(a>1,:),eventChannels(a>1),'spike',overwriteFlag);
    catch
    end
end


%% spatialIntegration

%% DETECT SPIKES WITH SPIKE KU

params.spkku.SPKDURATION = 200;
params.spkku.ABSTHRESH = 250;
params.spkku.AFTDUR = 100;
params.spkku.LLTHRESH = 9;
warning('off','all')
%find all spikes
parfor i =1:numel(datasetNames)
    try
        fprintf('Detecting spikes in : %s -\n',datasetNames{i});
        if (exist(fullfile(cd, sprintf('%s_spiketimes_LL.mat',datasetNames{i})), 'file') == 0)
            session = IEEGSession(datasetNames{i},userID,pwdFile);
            [spikeTimes, spikeChannels] = spike_ku_wrapper(session.data,ch_keep_final{i},params);
            parsave(sprintf('%s_spiketimes_kuv4.mat',session.data.snapName),spikeTimes,spikeChannels);
        else
            fprintf(' Mat found.\n')
        end
    catch
        fprintf('...error!\n');
    end
end

%% upload
uploadAnnotations(session.data,'spike_candidates',spikeTimes,num2cell(spikeChannels),'spikes_ja')
%% train classifier
fn = 28
% get true layer
[~, trueSpikesTimes tchannels] = getAnnotations(session.data(fn),'brian_oommen''s I022_P009_D01 annotations');
% get spike layer
[~, candSpikesTimes cchannels] = getAnnotations(session.data(fn),'spike_candidates');
% find all detected spikes within true layer, use as true set


idx = arrayfun(@(x)sum(x>trueSpikesTimes(:,1) &x<trueSpikesTimes(:,2)),candSpikesTimes(:,1));
idx = logical(idx);
trueSpikeSet = candSpikesTimes(idx,1);
trueSpikeCh = cchannels(idx);

falseSpikeSet = candSpikesTimes(~idx,1);
falseSpikeCh = cchannels(~idx);


%% for each true and false spike, get waveform
fs = session.data(fn).sampleRate;
timeBefore = 0.5;
timeAfter = 0.5 ;
scales = 1:60;

allSpikeEEG = zeros(size(candSpikesTimes,1),(timeBefore+timeAfter)*fs + 1);
for i = 1:size(candSpikesTimes,1)
    tmpdata = session.data(fn).getvalues((candSpikesTimes(i,1)/1e6 - timeBefore)*fs : (candSpikesTimes(i,1)/1e6 + timeAfter)*fs,cchannels{i});
%     c = cwt(tmpdata(:,1),scales,'mexh',1,'scal');
%     if idx(i)==1
%         fprintf('True Spike\n');
%     else
%         fprintf('False spike\n');
%     end
%     set(gca,'FontSize',20);
%     a = input('Any string to continue: ','s');
%     clf
    allSpikeEEG(i,:) = tmpdata;
end

%% feature extraction
coeffs = zeros(size(allSpikeEEG,1),numel(scales)*size(allSpikeEEG,2));
for i = 1:size(allSpikeEEG,1)
    coeffs(i,:) = reshape(cwt(allSpikeEEG(i,:),scales,'mexh'),1,numel(scales)*size(allSpikeEEG,2));
end
coeffs = abs(coeffs);
norm_coeffs = bsxfun (@minus, coeffs, min(coeffs,[],2));
divnorm = (max(coeffs,[],2)-min(coeffs,[],2));
norm_coeffs = bsxfun (@rdivide, norm_coeffs, divnorm);
spike_norm_coeffs = norm_coeffs(:,(timeBefore-0.05)*fs:end-(timeAfter+0.2)*fs);
colMeans = mean(spike_norm_coeffs);
[evectors, score, evalues] = pca(spike_norm_coeffs);

feats = score;

%% set division
trueSpikeSet = feats(idx,1:60);
trueSpikeCh = cchannels(idx);

falseSpikeSet = feats(~idx,1:60);
falseSpikeCh = cchannels(~idx);

%% plot pc1 vs 2
scatter(trueSpikeSet(:,1),trueSpikeSet(:,2),'o');
hold on;
scatter(falseSpikeSet(:,1),falseSpikeSet(:,2),'x');

%% cv classify
labels = [ones(size(trueSpikeSet,1),1); ones(size(falseSpikeSet,1),1)*0];
trainingset = [trueSpikeSet;falseSpikeSet];

mod = TreeBagger(1000,trainingset,labels,'method','Classification','OOBPredictorImportance','on','Cost',[0 1; 25 0]);
oobErrorBaggedEnsemble = oobError(mod);
plot(oobErrorBaggedEnsemble)
xlabel 'Number of grown trees';
ylabel 'Out-of-bag classification error';


[yhat,scores] = oobPredict(mod);
[conf, classorder] = confusionmat(categorical(labels), categorical(yhat))
disp(dataset({conf,classorder{:}}, 'obsnames', classorder));

imp = mod.OOBPermutedPredictorDeltaError;
predictorNames = {};
for i = 1:60
    predictorNames{i} = sprintf('%d',i');
end
figure;
bar(imp);
ylabel('Predictor importance estimates');
xlabel('PC');
h = gca;
h.XTick = 1:2:60
h.XTickLabel = predictorNames
h.XTickLabelRotation = 45;
h.TickLabelInterpreter = 'none';


addpath(genpath('../../Libraries/libsvm/'))
cv_acc = svmtrain(labels,trainingset,'-s 0 -t 0 -v 10');
mod = svmtrain(labels,trainingset,'-s 0 -t 0');

%% test on new dataset
params.spkku.SPKDURATION = 200;
params.spkku.ABSTHRESH = 250;
params.spkku.AFTDUR = 100;
params.spkku.LLTHRESH = 9;
params.spkku.pca.colmeans = colMeans;
params.spkku.scales = 1:60

i = 3
[spikeTimes spikeChannels] = detect_spikes(session.data(i),ch_keep_final{i},params,mod)
spikeTimes = [spikeTimes(:,1)-(0.05*1e6) spikeTimes+(0.05*1e6)];
uploadAnnotations(session.data(i),'detected_spikes',spikeTimes,num2cell(spikeChannels),'spike')

i = 1
[spikeTimes spikeChannels] = detect_spikes(session.data(i),ch_keep_final{i},params,mod)
uploadAnnotations(session.data(i),'detected_spikes',spikeTimes,num2cell(spikeChannels),'spike')

i = 4
[spikeTimes spikeChannels] = detect_spikes(session.data(i),ch_keep_final{i},params,mod)
uploadAnnotations(session.data(i),'detected_spikes',spikeTimes,num2cell(spikeChannels),'spike')

%% CLUSTER SPIKES: TO BE DEVELOPED
[spikeWV, spikeInfo] = loadSpikesMulti_par(params,'spikekuv4_WV','spiketimes_kuv4');


%view spikes
for i = 1:size(spikeWV,1)
    plot(spikeWV(i,:));
    hold on;
    pause(.1)
end


[idx b] = kmeans(spikeWV,50,'MaxIter',300);
[idx2 b2] = kmeans(spikeWV,50,'MaxIter',300,'Distance','cityblock');

res = plotExFromCluster(session.data,spikeInfo,idx,1);


[~, score] = pca(spikeWV);
[idx b] = kmeans(score(:,1:4),25,'MaxIter',300);

for i = 1:25
   scatter(score(idx==i,1),score(idx==i,2),1)
   hold on;
end
session = IEEGSession(datasetNames{1},userID,pwdFile);
for i = 2:numel(datasetNames)
    session.openDataSet(datasetNames{i})
end

tmpidx = randperm(size(score,1));
tmpidx = tmpidx(1:round(numel(tmpidx)*.1))
res= reviewClusters(session.data,spikeInfo,tmpidx);


[x y] = knnsearch(spikeWV(res(:,1),:),spikeWV);
resSpkIdx = find(res(:,2)==1);
spikeIdx = ismember(x,resSpkIdx)
res = plotExFromCluster(session.data,spikeInfo,spikeIdx,1);



[x y] = knnsearch(score(1:163,:),score(164:end,:));
res=  [];
for i = 1:max(idx)
    res = [res; plotExFromCluster(session.data,spikeInfo,idx,i)];
end

%save('initialCluster.mat','idx');
load('initialCluster.mat');
%cluster PCS
[a b c] = pca(spikes);
E = evalclusters(b(:,1:2),'kmeans','GAP','klist',[1:30]);
[idx b] = kmeans(b(:,1:2),E.OptimalK);

%% Plot and Remove first
t = 1/400:1/400:.2;
for i = 1:max(idx)
    tmp = spikeWV(idx==i,:);
    subaxis(floor(sqrt(max(idx))),ceil(sqrt(max(idx))),i,'SpacingVert',0.04,'SpacingHoriz',0);
    errorbar(t,mean(tmp),std(tmp),'Color',[.7 .7 .7])
    hold on;
    plot(t,mean(tmp),'r');
    xlim([0 .2])
    
%     for j = 1:size(tmp,1)
%         plot(tmp(j,:));
%         hold on;
%     end
    title(['Cluster ' num2str(i) ': ' num2str(size(tmp,1))])
end

%remove 
removeClusters = [15 30 41 50 52 56 64 68 82 84 90 99];
removeIdx = ismember(idx,removeClusters);
spikeWV(removeIdx,:) = [];
spikeInfo(removeIdx,:) = [];
idx(removeIdx,:) = [];

%% recluster
load('secondroundcluster.mat');
% [idx b] = kmeans(spikeWV,100);
% for i = 1:max(idx)
%     tmp = spikeWV(idx==i,:);
%     subaxis(floor(sqrt(max(idx))),ceil(sqrt(max(idx))),i,'SpacingVert',0.04,'SpacingHoriz',0);
%     errorbar(t,mean(tmp),std(tmp),'Color',[.7 .7 .7])
%     hold on;
%     plot(t,mean(tmp),'r');
%     xlim([0 .2])
%     
% %     for j = 1:size(tmp,1)
% %         plot(tmp(j,:));
% %         hold on;
% %     end
%     title(['Cluster ' num2str(i) ': ' num2str(size(tmp,1))])
% end
removeClusters = [3 39 83 94];
removeIdx = ismember(idx,removeClusters);
spikeWV(removeIdx,:) = [];
spikeInfo(removeIdx,:) = [];
idx(removeIdx,:) = [];

%% recheck
% [idx b] = kmeans(spikeWV,100);
% for i = 1:max(idx)
%     tmp = spikeWV(idx==i,:);
%     subaxis(floor(sqrt(max(idx))),ceil(sqrt(max(idx))),i,'SpacingVert',0.04,'SpacingHoriz',0);
%     errorbar(t,mean(tmp),std(tmp),'Color',[.7 .7 .7])
%     hold on;
%     plot(t,mean(tmp),'r');
%     xlim([0 .2])
%     
% %     for j = 1:size(tmp,1)
% %         plot(tmp(j,:));
% %         hold on;
% %     end
%     title(['Cluster ' num2str(i) ': ' num2str(size(tmp,1))])
% end

%%
%[idx1] = kmeans(spikeWV,36,'MaxIter',400);
%save('kmeans36maxiter400wv.mat','idx1','spikeWV','spikeInfo');
% [p1, pc, p2] = pca(spikeWV);
% [idxp] = kmeans(pc(:,1:20),49,'MaxIter',400);
% save('kmeans49maxiter400pca20.mat','idxp','spikeWV','spikeInfo');

%% choose clusters
load('kmeans36maxiter400wv.mat');
figure;
t = 1/400:1/400:.2;
for i = 1:max(idx1)
    tmp = spikeWV(idx1==i,:);
    subaxis(floor(sqrt(max(idx1))),ceil(sqrt(max(idx1))),i,'SpacingVert',0.04,'SpacingHoriz',0);
    errorbar(t,mean(tmp),std(tmp),'Color',[.7 .7 .7])
    hold on;
    plot(t,mean(tmp),'r');
    xlim([0 .2])
    
%     for j = 1:size(tmp,1)
%         plot(tmp(j,:));
%         hold on;
%     end
    title(['Cluster ' num2str(i) ': ' num2str(size(tmp,1))])
end
figure;

for i = 1:max(idx)
    res = [res plotExFromCluster(datasetNames,spikeInfo,idx1,i)];
end

load('kmeans49maxiter400pca20.mat');
figure;
t = 1/400:1/400:.2;
for i = 1:max(idx1)
    tmp = spikeWV(idx1==i,:);
    subaxis(floor(sqrt(max(idx1))),ceil(sqrt(max(idx1))),i,'SpacingVert',0.04,'SpacingHoriz',0);
    errorbar(t,mean(tmp),std(tmp),'Color',[.7 .7 .7])
    hold on;
    plot(t,mean(tmp),'r');
    xlim([0 .2])
    
%     for j = 1:size(tmp,1)
%         plot(tmp(j,:));
%         hold on;
%     end
    title(['Cluster ' num2str(i) ': ' num2str(size(tmp,1))])
end
figure;
resp= [];
for i = 1:max(idxp)
resp = [resp plotExFromCluster(datasetNames,spikeInfo,idxp,i)];
end

keep = str2num(input('Input clusters to keep: ','s'));
idx = ismember(idx,keep);

%plot all spikes
% for i = 1:size(spikes2,1)
%     subplot(2,1,1);
%     plot(spikes(i,:));
%     hold on;
%     subplot(2,1,2);
%     plot(spikes2(i,:));
%     hold all;
% end


% %plot all spikes
% 
% %gaussian mixture model
% 
% %find optimum number of clusters
% % eva = evalclusters(spikes,'gmdistribution','gap','Klist',[3 4 5 6]);
% % plot(eva)
options = statset('Display','final');
gm = gmdistribution.fit(spikes2,4,'Options',options);
% 
% idx = cluster(gm,spikes2);
% for i = 1:max(idx)
%     figure(i);
%     clf;
%     tmp = spikes2(idx==i,:);
%     for j = 1:size(tmp,1)
%         plot(tmp(j,:));
%         hold on
%     end
%     plot(mean(tmp),'r','LineWidth',2)
% end




%        uploadAnnotations(session.data(i), 'spike_LL',spikeTimes(idx,:),spikeChannels(idx),'spike');       
%         fprintf('...done!\n');
%         spikeChannels = num2cell(spikeChannels);
%         save(sprintf('%s_spikes.mat',session.data(i).snapName),'spikeTimes','spikeChannels');
%     catch
%         disp('error');
%     end
% end
    %match to template spikes

%     fs = session.data(i).sampleRate;
%     template = load('template_spikes.mat');
%     
%     %validate template spikes
%     figure;
%     keepIdx = [];
%     for j = 1:size(template.spikes2,1);
%         clf;
%         x = size(template.spikes2,2);
%         x = (1:x)/fs;
%         plot(x,template.spikes2(j,:));
%         title(num2str(j));
%         keep = input('Keep?: ','s');
%         if strcmp(keep,'y')
%             keepIdx = [keepIdx j];
%         end
%     end
%     template.spikes2 = template.spikes2(keepIdx,:);
%        
%     
%     idx = zeros(size(spikeTimes,1),size(template.spikes2,1));
%     for j = i:size(spikeTimes,1)
%         startPt = max((spikeTimes(j,1)/1e6-0.04)*fs,1);
%         endPt = startPt+0.2*fs;
%         dat = session.data(i).getvalues(startPt:endPt,spikeChannels{j});
%         for k = 1:size(template.spikes2,1);
%             [C, Lags] = xcorr(dat,template.spikes2(k,:),'coeff');
%             if max(abs(C)) > 0.9
%                 %subplot(3,1,1);
%                 %plot(template.spikes2(k,:));
%                 %xlim([0 81]);
%                 %subplot(3,1,2);
%                	%shift = Lags(C==max(C));
%                 %x = [1:81]+shift;
%                 %plot(x,dat);
%                 %xlim([0 81]);
%                 %subplot(3,1,3);
%                 %plot(C);
%                 idx(j,k) = max(abs(C));
%             end
%         end
%     end
%     imagesc(corrmat);
%     colormap('hot');
%     colorbar;

     
     
    % idx = clusterSpikes(session.data(i),spikeTimes, spikeChannels);
       
   
     %standalone_spike_keating_std(session.data(i),'spike_keating_auto',1200,channelIdxs{i},7,180,80); 
    % standalone_spike_LL(session.data(i),'spike_LL_filtered',channelIdxs{i},.1,4,1);
%          [a, b, TP, FP] = evalDetections(session.data(i),'dd_MARKS','spike_keating_auto')
%          [a, b, TP, FP] = evalDetections(session.data(i),'dd_MARKS','spike_LL_filtered')
%          [a, b, TP, FP] = evalDetections(session.data(i),'dd_MARKS','spike_AR')


%% cluster spikes
[annots, timesUSec, eventChannels] = getAllAnnots(session.data,'brian_oommen''s I022_P009_D01 annotations');
peakTimes = zeros(timesUSec,1);

 

% generate stat table
%Patient | Trial | Word order | Word | TimeToRecall | Recalled | NumSpikes 
eventAnnotationName = 'pyFR_Events';
spikeLayerName = 'spikes';
parfor i =1:numel(datasetNames)
    session = IEEGSession(datasetNames{i},userID,pwdFile);
    [annots, timesUSec, ch] = getAllAnnots(session.data,eventAnnotationName);
    desc = {annots.description};
    desc2 = {annots.type};
    words = cellfun(@(x)x(find(x==':',1)+2:find(x==';',1)-2),desc2,'UniformOutput',0);
    tRecalled = cellfun(@(x)str2num(x(find(x==':',1,'last')+2:end)),desc2,'UniformOutput',0); %extract recalled: 1 = true, 0 = false, -999 = NA
    %count spikes
    fs = session.data.sampleRate;
    startIdxs = startT/1e6*fs;
    tPatient = cell(size(startT,1),1);
    tTrial = cell(size(startT,1),1);
    tWordOrder = cell(size(startT,1),1);
    tWord = cell(size(startT,1),1);
    tIsRecWord = cell(size(startT,1),1);
    tTimeToRecall = cell(size(startT,1),1);
    tNumSpikes = cell(size(startT,1),1);
    
    currentTrial = 0
    for w = 1:numel(startIdxs) % for each annotation
        recFlag = 0;
        switch desc{w}
            case 'Trial'
                currentTrial = currentTrial+1;
                wordOrder = 0
            case 'REC_START'
                rec_start = timesUsec(w,1);
                wordOrder = 0
            case 'REC_WORD'
                tTimeToRecall{w} = timesUSec(w,1) - rec_start
        end
        wordOrder = wordOrder + 1;
        tPatient{w} = w;
        tTrial{w} = currentTrial;
        tWordOrder{w} = wordOrder;
        
        %find total spikes
        
    end
        
end