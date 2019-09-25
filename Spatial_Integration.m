function [EventTimeUSec, EventChannels] = Spatial_Integration(spikeTimesUSec, spikeChannels,halfSearchRangeInMiliSec)
% Hanlin,June 17 2015
% Hanlin,June 18 2015
% spikeTimesUSec:  a 1D vector recording the timesInUSec where spikes are detected.
% spikeChannels:    a 1D vector having the same size as spikeTimesUSec,
%                               telling us which corresponding channels the spikes are. 
% halfSearchRange:   unit: Miliseconds.  this variable specify the time range we search for other spikes in other channels whenever a spike is found in one of the channels at a specific time: 
halfSearchRange=halfSearchRangeInMiliSec*1000; % convert to microSeconds; 
[SortedTime,SortIdx]=sort(spikeTimesUSec);  % sort the 2 input variable in ascending order
SortedChannels = spikeChannels(SortIdx);
EventTimeUSec=0;   % output initialized
EventChannels=cell(1,1); 

DiffList=diff(SortedTime)<=halfSearchRange; % apply searching criteria. u
if numel(DiffList)==0
    return
end
%DiffList(2:end+1)=DiffList;%DiffList(1)=DiffList(2);  % in the next step we lost the first element, this step compensates for that.
%MultiList=DiffList(1:end-1).*DiffList(2:end); % each element here dictates whether to start or new record or not

% the first input record is also the first output record. 
EventTimeUSec(1)=SortedTime(1);
EventChannels{1} = SortedChannels(1);
% now determine if other input records can be grouped or should be kept as individual
for i=1:numel(DiffList)
     if logical(DiffList(i))  % group the channel record to previous ones, ignore event timing record. 
        EventChannels{end,1}=[EventChannels{end} SortedChannels(i+1)];
         continue;
     end
     %start a new record.
EventTimeUSec(end+1,1)=SortedTime(i+1);
EventChannels{end+1,1} = SortedChannels(i+1);    
end
end