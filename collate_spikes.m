function collate_spikes(userID, pwdFile)
%Usage: run_spike_detection(userID,pwdFile)
% userID - IEEG.org username
% pwdFile - password file for associated username
% Program will collate spikes all datasets included in study and write
% results to table

session = IEEGSession('I022_P009_D01', userID, pwdFile);
subj = [1:2 4:17];
subjsess = 1:6;
for i = subj
    for j = subjsess
        dataName = sprintf('I022_P0%0.2d_D%0.2d',i,j);
        try
            session.openDataSet(dataName);
        catch
            fprintf('Unable to load %s\n',dataName);
        end
    end
end
subj =  1:54;
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

%%
% Create table consisting of Free Recall Task and specified annotation data
allFn = [{session.data.snapName}];
ptFn = cellfun(@(x)x(1:9),allFn(2:size(allFn,2)),'UniformOutput',0); %size(allFn,2)
uniquePtFn = unique(ptFn);
for i = 1:numel(uniquePtFn) %for each patient
    try
        dIdx= regexp(allFn,uniquePtFn(i));
        dIdx = cellfun(@(x)~isempty(x),dIdx);
        datasets = session.data(dIdx);
        
        % Load Seizure onset zone channels from premade .csv 
        fid = fopen('Szchannels.csv');
        a = textscan(fid,'%s%s','Delimiter',',');
        patientID = uniquePtFn{i};
        b = strfind(a{1,1}, patientID);
        patientIndx =find(~cellfun(@isempty,b));
        Szchannels = strsplit(a{1,2}{patientIndx,1},';');
        fclose(fid);
        
        % Load patient demographic info from premade .csv
        fid = fopen('pt_dem.csv');
        a = textscan(fid,'%s%s%s','Delimiter',',');
        patientID = uniquePtFn{i};
        b = strfind(a{1,1}, patientID);
        patientIndx =find(~cellfun(@isempty,b));
        age = a{1,2}{patientIndx,1};
        if isempty(age)
            age = NaN;
        end
        sex = a{1,3}{patientIndx,1};
        if isempty(sex)
            sex = NaN;
        end
        fclose(fid);
        
        
        Tab = cell2table({}); % initialize empty table for patient's data
        for j = 1:numel(datasets)%1:numel(datasets)
            dataset = datasets(j);
            sessnumber = dataset.snapName(12:13);
            fprintf('Extracting %s \n',dataset.snapName);
            % Extract spike and Free Recall (FR) task annotations
            try
                [FR_allEvents, FR_timesUSec, FR_channels] = getAllAnnots(dataset,'pyFR_Events');
                [spike_spat_allEvents, spike_spat_timesUSec, spike_spat_channels] = getAllAnnots(dataset,'spike_ja_spatial');
            catch
                fprintf('Missing annotation layer in session %s\n',r_sessindx(i))
            end
            
            % Create array of start times for SzZone spike annotations only
            SzchanIndxs = [];
            for r = 1:size(Szchannels, 2)
                c = strfind(dataset.channelLabels(:,1), Szchannels{1,r});
                chanIndx = find(~cellfun(@isempty,c));
                SzchanIndxs = [SzchanIndxs chanIndx];
            end
                       
            SzZonecheck = zeros(size(spike_spat_channels,1),1);
            if strcmp(dataset.snapName, 'I027_P012_D01') || (strcmp(uniquePtFn{i}, 'I027_P004') && ismember(j,1:7))
                A = find(SzZonecheck);
            else
                for g = 1:size(spike_spat_channels,1)
                    SzZonecheck(g) = any(ismember(SzchanIndxs, spike_spat_channels{g,1}));
                end
            end
            A = find(SzZonecheck);
            Szspike_timesUSec =  spike_spat_timesUSec(A,1);
            
            % Create cell array of spikes between each Free Recall task annotation
            [trash, edges_s, ~] = histcounts(spike_spat_timesUSec(:,1), FR_timesUSec(:,1));
            [trash, edges_sz, ~] = histcounts(Szspike_timesUSec, FR_timesUSec(:,1));

            spikes_from_rec = 0; % start summation of spike_ja_spatial from recall task start
            spikes_from_rec_Sz = 0; % start summation of spike_Sz from recall task start
            spikes_from_rec_nSz = 0; % start summation of spike_nSz from recall task start
            
            % Initialize cell matrix of FR annotation information (i.e. TASK, WORD, ONSET,
            % RECALLED, etc.)
            FR_annots = num2cell(zeros(size(FR_channels,1),10)); % zero cell matrix of FR task annotations
            x = 0; % TRIAL #
            order = 0; % WORD ORDER within each TRIAL
            
            for p = 1:size(FR_annots,1)
                % Extract TASK from FR_allEvents
                FR_annots{p,1} = FR_allEvents(1,p).description; % B
                % Extract WORD from FR_allEvents
                [wordstart, wordend] = regexp(FR_allEvents(1,p).type,'\<[A-Z]{2,8}\>');
                if isempty(wordstart)
                    FR_annots{p,2} = 'X';
                else
                    FR_annots{p,2} = FR_allEvents(1,p).type(wordstart:wordend);
                end
                % Extract stimulus ONSET from FR_timesUSec
                FR_annots{p,3} = FR_timesUSec(p,1);
                % Identify whether there was RECALL '1' or not '0', '2'
                recall_value = regexp(FR_allEvents(1,p).type,'\<[0-9]\>');
                if isempty(recall_value)
                    FR_annots{p,4} = 2;
                else
                    FR_annots{p,4} = str2double(FR_allEvents(1,p).type(recall_value));
                end
                % Identify TRIAL #
                if strcmp(FR_annots{p,1},'TRIAL') && p == size(FR_annots,1)
                    FR_annots{p,5} = 0;
                elseif strcmp(FR_annots{p,1},'TRIAL')
                    x = x + 1;
                    FR_annots{p,5} = x;
                else
                    FR_annots{p,5} = x;
                end
                % Identify presentation ORDER of WORD stimuli within each block
                if strcmp(FR_annots{p,1},'REC_WORD_VV') && (strcmp(FR_annots{p,2},'X') || strcmp(FR_annots{p,2},'VV'))
                    FR_annots{p,6} = 0;
                elseif strcmp(FR_annots{p,1},'REC_WORD') && strcmp(FR_annots{p,2},'X')
                    FR_annots{p,6} = 0;
                elseif strcmp(FR_annots{p,2},'X')
                    order = 0;
                    FR_annots{p,6} = order;
                else
                    order = order + 1;
                    FR_annots{p,6} = order;
                end
                % Identify TIME FROM RECALL START (from REC_START to each REC_WORD)
                if (strcmp(FR_annots{p,1},'REC_WORD') || strcmp(FR_annots{p,1},'REC_WORD_VV')) && strcmp(FR_annots{p-1,1},'WORD')
                    rec_start_time = 0;
                elseif strcmp(FR_annots{p,1},'REC_START')
                    rec_start_time = FR_timesUSec(p);
                elseif strcmp(FR_annots{p,1},'REC_WORD') || strcmp(FR_annots{p,1},'REC_WORD_VV')
                    FR_annots{p,7} = FR_timesUSec(p) - rec_start_time;
                else
                    FR_annots{p,7} = 0;
                end
            end
            
            % Modify SPIKE_COUNT to exclude spikes after a WORD's ENCODING window
            edges_mod = edges_s;                
            edges_modsz = edges_sz;
            try
                for t = 2:size(FR_annots,1)
                    if (strcmp(FR_annots{t,1},'REC_START') || strcmp(FR_annots{t,1},'SESS_START') ||...
                            ((strcmp(FR_annots{t,1},'REC_WORD_VV') || strcmp(FR_annots{t,1},'REC_WORD')) && strcmp(FR_annots{t-1,1},'WORD')))...
                            && ((FR_annots{t,3} - FR_annots{t-1,3}) > 2.8e6)
                        edges_mod(t) = FR_annots{t-1,3} + 2.8e6;
                        edges_modsz(t) = FR_annots{t-1,3} + 2.8e6;                    
                    end
                end
            catch
                fprintf('Cannot exclude spike in session %s\n',t)
            end
            [spike_s_num, ~, ~] = histcounts(spike_spat_timesUSec(:,1), edges_mod(:,1));
            [spike_Sz_num, ~, ~] = histcounts(Szspike_timesUSec, edges_modsz);
            
            % Modify SPIKE_COUNT to match RECALL indices
            spike_num2 = zeros(size(spike_s_num,2),1);
            spike_num2_Sz = zeros(size(spike_Sz_num,2),1);            
            for y = 1:size(spike_s_num,2)
                if (strcmp(FR_annots{y,1},'REC_WORD') || strcmp(FR_annots{y,1},'REC_WORD_VV')) && strcmp(FR_annots{y-1,1},'REC_START')
                    spike_num2(y) = spike_s_num(y);
                    spike_num2_Sz(y) = spike_Sz_num(y);                   
                elseif strcmp(FR_annots{y,1},'REC_WORD') || strcmp(FR_annots{y,1},'REC_WORD_VV')
                    spike_num2(y) = spike_s_num(y-1);
                    spike_num2_Sz(y) = spike_Sz_num(y-1);                 
                else
                    spike_num2(y) = spike_s_num(y);
                    spike_num2_Sz(y) = spike_Sz_num(y);                    
                end
            end
            
            % Create array of spikes from non-onset seizure channels
            spike_num2_nSz = spike_num2 - spike_num2_Sz;           
            %Pad 0 to match FR_annots size
            spike_count = num2cell([spike_num2; 0]);
            spike_count_Sz = num2cell([spike_num2_Sz; 0]);
            spike_count_nSz = num2cell([spike_num2_nSz; 0]);

            %Identify TOTAL SPIKES FROM RECALL START for each recall trial
            for p = 1:size(FR_annots,1)-1
                if FR_annots{p,7} > 0
                    spikes_from_rec = spikes_from_rec + spike_num2(p);
                    FR_annots{p,8} = spikes_from_rec;
                    spikes_from_rec_Sz = spikes_from_rec_Sz + spike_num2_Sz(p);
                    FR_annots{p,9} = spikes_from_rec_Sz;
                    spikes_from_rec_nSz = spikes_from_rec_nSz + spike_num2_nSz(p);
                    FR_annots{p,10} = spikes_from_rec_nSz;
                else
                    FR_annots{p,8} = 0;
                    spikes_from_rec = 0;
                    FR_annots{p,9} = 0;
                    spikes_from_rec_Sz = 0;
                    FR_annots{p,10} = 0;
                    spikes_from_rec_nSz = 0;
                end
            end

            %Create structure array of session's annotation and spike data
            FR_trials = size(FR_allEvents,2);
            patient_struct = struct('patient', repmat({dataset.snapName(1:9)}, FR_trials, 1),...
                'age', repmat({age}, FR_trials,1),...
                'sex', repmat({sex}, FR_trials,1),...
                'session', repmat({dataset.snapName(11:13)}, FR_trials, 1),...
                'task', FR_annots(:,1), 'trial', FR_annots(:,5),'word', FR_annots(:,2),...
                'word_order', FR_annots(:,6),...
                'onset', FR_annots(:,3),'recalled', FR_annots(:,4), 'time_from_recall_start', FR_annots(:,7),...
                'spikes_ja_spatial', spike_count,...
                'spikes_ja_spatial_from_recall_start', FR_annots(:,8),...
                'SzOnset_spikes', spike_count_Sz,...
                'SzOnset_spikes_from_recall_start', FR_annots(:,9),...
                'non_SzOnset_spikes', spike_count_nSz,...
                'non_SzOnset_spikes_from_recall_start', FR_annots(:,10));
            
            %Convert to table
            try
                D = struct2table(patient_struct);
                % Add Brodmann areas/hemisphere side
                T = areaChannels6(patientID,spike_spat_channels,edges_mod, FR_annots, spike_spat_timesUSec, spike_s_num,dataset);
                C = horzcat(D, T);
                % Align columns for datasets where new electrodes were
                % implanted after first session
                if (strcmp(uniquePtFn{i}, 'I027_P012') && strcmp(dataset.snapName, 'I027_P012_D01')) ||...
                        (strcmp(uniquePtFn{i}, 'I027_P004') && strcmp(dataset.snapName, 'I027_P004_D01'))
                    headings = C.Properties.VariableNames;
                end
                try
                    Tab = vertcat(Tab, C);
                catch
                    headings2 = C.Properties.VariableNames;
                    [~,headind] = setdiff(headings2,headings);
                    C2 = C;
                    C2(:,headind) = [];
                    Tab = vertcat(Tab, C2);
                end
            catch
                fprintf('Unable to tabulate %s\n',sessnumber)
            end
        end
    catch
        fprintf('Unable to analyze session %s\n',sessnumber)
    end
    filename = [uniquePtFn{i} '.csv'];
     % Save merged session patient file
    writetable(Tab, filename, 'Delimiter', ',');
end
