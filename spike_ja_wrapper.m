function [spikeTimes, spikeChannels, DE] = spike_ja_wrapper(IEEGDataset,channels,params)
    fs = IEEGDataset.sampleRate;
    eventIdxs = cell(numel(channels),1);
    eventChannels = cell(numel(channels),1);
    dataLim = 500*130*2000;
    durInPts = IEEGDataset.rawChannels(1).get_tsdetails.getDuration/1e6*fs;
    data = getAllData(IEEGDataset,channels,3600);
    [DE]=spike_detector_hilbert_v16_byISARG(data,fs, '-h 60');
    spikeTimes = DE.pos*1e6;
    spikeChannels = channels(DE.chan);
    fprintf('%d spikes found\n',size(spikeTimes,1));
end



