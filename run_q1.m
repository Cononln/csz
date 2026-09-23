function run_q1(sourceDir, outputDir)
%RUN_Q1 Reproducible analysis for C-question 1 using raw Fz/F3/F4 channels.
%   run_q1(sourceDir, outputDir) reads the four VisualCog*.mat files,
%   extracts cue-locked epochs, flags suspect epochs, estimates robust ERPs,
%   and evaluates left/right cue information across recordings.

    arguments
        sourceDir (1,1) string
        outputDir (1,1) string
    end
    if ~isfolder(sourceDir)
        error('Source directory does not exist: %s', sourceDir);
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    figureDir = fullfile(outputDir, "figures");
    if ~isfolder(figureDir)
        mkdir(figureDir);
    end

    cfg = default_config();
    rng(cfg.seed, 'twister');
    files = ["VisualCogA_Task-1.mat", "VisualCogA_Task-2.mat", ...
             "VisualCogB_Task-1.mat", "VisualCogB_Task-2.mat"];
    allEvents = table();
    qualitySummary = table();
    erpCurves = table();
    windowStatistics = table();
    reliability = table();
    contrastReliability = table();
    filterSensitivity = table();
    lateWave = table();
    lateWaveCommon = table();
    spatialCurves = table();
    records = struct('name', {}, 'task', {}, 'features', {}, ...
                     'spatialFeatures', {}, ...
                     'cue', {}, 'accepted', {}, 'epochs', {}, 'time', {});

    for fileIdx = 1:numel(files)
        sourceFile = fullfile(sourceDir, files(fileIdx));
        if ~isfile(sourceFile)
            error('Missing source file: %s', sourceFile);
        end
        s = load(sourceFile, 'data', 'SampleRate', 'DataLabel');
        validate_source(s, files(fileIdx));
        fs = s.SampleRate;
        [name, task] = parse_record_name(files(fileIdx));
        time = (-round(cfg.preSec*fs):round(cfg.postSec*fs))/fs;
        eventTable = parse_events(s.data, fs, name, task);
        [epochs, isClipped, peakToPeak, maxJump, inBounds] = ...
            make_epochs(s.data, eventTable.cueIndex, fs, cfg, time);
        [accepted, highPeak, highJump, peakThreshold, jumpThreshold] = ...
            mark_artifacts(isClipped, peakToPeak, maxJump, inBounds, cfg);

        eventTable.inBounds = inBounds;
        eventTable.rawClipped = isClipped;
        eventTable.peakToPeak = peakToPeak;
        eventTable.maxJump = maxJump;
        eventTable.highPeak = highPeak;
        eventTable.highJump = highJump;
        eventTable.accepted = accepted;
        allEvents = [allEvents; eventTable]; %#ok<AGROW>

        summaryRow = table(name, task, height(eventTable), ...
            sum(eventTable.cueDir == -1), sum(eventTable.cueDir == 1), ...
            sum(isClipped), sum(highPeak), sum(highJump), sum(accepted), ...
            sum(accepted & eventTable.cueDir == -1), ...
            sum(accepted & eventTable.cueDir == 1), ...
            peakThreshold, jumpThreshold, ...
            'VariableNames', {'record','task','trials','leftCues','rightCues', ...
            'rawClipped','highPeak','highJump','accepted','acceptedLeft', ...
            'acceptedRight','peakThreshold','jumpThreshold'});
        qualitySummary = [qualitySummary; summaryRow]; %#ok<AGROW>

        if any(sum(accepted & eventTable.cueDir == [-1, 1], 1) < cfg.minTrialsPerClass)
            error('Too few retained trials in %s. Inspect quality thresholds.', name);
        end
        [curves, fitted] = estimate_curves(epochs, eventTable.cueDir, ...
            accepted, time, name, task, cfg);
        erpCurves = [erpCurves; curves]; %#ok<AGROW>
        spatialCurves = [spatialCurves; spatial_mode_curves( ...
            fitted,time,name,task)]; %#ok<AGROW>
        windowStatistics = [windowStatistics; window_stats(epochs, ...
            eventTable.cueDir, accepted, time, name, task, cfg)]; %#ok<AGROW>
        reliability = [reliability; split_half(epochs, eventTable.cueDir, ...
            accepted, time, name, task, cfg)]; %#ok<AGROW>
        contrastReliability = [contrastReliability; contrast_split_half( ...
            epochs, eventTable.cueDir, accepted, time, name, task, cfg)]; %#ok<AGROW>
        filterSensitivity = [filterSensitivity; filter_sensitivity( ...
            s.data, eventTable.cueIndex, eventTable.cueDir, accepted, ...
            fs, name, task, cfg)]; %#ok<AGROW>
        [lateRows, commonRows] = late_wave_diagnostic(epochs, ...
            eventTable.cueDir, accepted, time, name, task, cfg);
        lateWave = [lateWave; lateRows]; %#ok<AGROW>
        lateWaveCommon = [lateWaveCommon; commonRows]; %#ok<AGROW>
        draw_erp(time, fitted, name, fullfile(figureDir, name + "_erp.png"));

        records(fileIdx).name = name;
        records(fileIdx).task = task;
        records(fileIdx).features = extract_features(epochs, time, cfg);
        records(fileIdx).spatialFeatures = ...
            extract_spatial_features(epochs, time, cfg);
        records(fileIdx).cue = eventTable.cueDir;
        records(fileIdx).accepted = accepted;
        records(fileIdx).epochs = epochs;
        records(fileIdx).time = time;
        fprintf('%s: %d/%d retained, %d clipped, left/right %d/%d.\n', ...
            name, sum(accepted), height(eventTable), sum(isClipped), ...
            sum(accepted & eventTable.cueDir == -1), ...
            sum(accepted & eventTable.cueDir == 1));
    end

    crossRecord = [evaluate_cross_record(records,cfg,'features'); ...
        evaluate_cross_record(records,cfg,'spatialFeatures')];
    crossRecord.holmP = holm_adjust(crossRecord.permutationP);
    [withinOriginal,foldsOriginal] = evaluate_within_record( ...
        records,cfg,'features');
    [withinSpatial,foldsSpatial] = evaluate_within_record( ...
        records,cfg,'spatialFeatures');
    withinRecord = [withinOriginal;withinSpatial];
    withinRecord.holmP = holm_adjust(withinRecord.circularShiftP);
    withinFolds = [foldsOriginal;foldsSpatial];
    writetable(allEvents, fullfile(outputDir, 'events_and_quality.csv'));
    writetable(qualitySummary, fullfile(outputDir, 'quality_summary.csv'));
    writetable(erpCurves, fullfile(outputDir, 'erp_curves.csv'));
    writetable(spatialCurves, fullfile(outputDir, 'spatial_mode_curves.csv'));
    writetable(windowStatistics, fullfile(outputDir, 'window_statistics.csv'));
    writetable(reliability, fullfile(outputDir, 'split_half_reliability.csv'));
    writetable(contrastReliability, fullfile(outputDir, ...
        'contrast_split_half.csv'));
    writetable(filterSensitivity, fullfile(outputDir, ...
        'filter_sensitivity.csv'));
    writetable(crossRecord, fullfile(outputDir, 'cross_record_decoding.csv'));
    writetable(withinRecord, fullfile(outputDir, 'within_record_decoding.csv'));
    writetable(withinFolds, fullfile(outputDir, 'within_record_folds.csv'));
    writetable(lateWave, fullfile(outputDir, 'late_wave_summary.csv'));
    writetable(lateWaveCommon, fullfile(outputDir, ...
        'late_wave_common_mode.csv'));
    save(fullfile(outputDir, 'q1_analysis.mat'), 'records', 'cfg', ...
         'allEvents', 'qualitySummary', 'windowStatistics', ...
         'contrastReliability', 'filterSensitivity', ...
         'withinRecord', 'withinFolds', 'spatialCurves', 'lateWave', ...
         'lateWaveCommon', '-v7.3');
    fprintf('Question 1 outputs saved to %s\n', outputDir);
end

function cfg = default_config()
    cfg.preSec = 0.2;
    cfg.postSec = 0.8;
    cfg.highPassHz = 0.5;
    cfg.lowPassHz = 30;
    cfg.clipLevel = 999;
    cfg.artifactMadMultiplier = 6;
    cfg.minTrialsPerClass = 15;
    cfg.huberConstant = 1.345;
    cfg.huberIterations = 15;
    cfg.smoothLambda = 12;
    cfg.bootstrapIterations = 500;
    cfg.splitIterations = 100;
    cfg.permutationIterations = 1000;
    cfg.withinRecordFolds = 5;
    cfg.seed = 20260922;
    cfg.windows = [0.10 0.25; 0.25 0.50];
    cfg.windowNames = ["early_100_250ms", "late_250_500ms"];
    cfg.electrodeNames = ["Fz", "F3", "F4"];
    cfg.diagnosticHighPassHz = [0.1, 0.5, 1.0];
    cfg.diagnosticPreSec = 0.8;
end

function validate_source(s, fileName)
    required = {'data','SampleRate','DataLabel'};
    if ~all(isfield(s, required))
        error('%s lacks required MAT variables.', fileName);
    end
    if ~isnumeric(s.data) || size(s.data,1) ~= 10 || ...
            ~isscalar(s.SampleRate) || s.SampleRate ~= 256
        error('%s has unexpected data dimensions or sampling rate.', fileName);
    end
    if any(~isfinite(s.data(:)))
        error('%s contains nonfinite source values.', fileName);
    end
    timestamp = s.data(10,:);
    if any(diff(timestamp) <= 0)
        error('%s has a nonmonotonic timestamp channel.', fileName);
    end
end

function [name, task] = parse_record_name(fileName)
    match = regexp(char(fileName), 'VisualCog([AB])_Task-([12])\.mat$', 'tokens', 'once');
    if isempty(match)
        error('Unexpected source filename: %s', fileName);
    end
    name = "VisualCog" + string(match{1}) + "_Task-" + string(match{2});
    task = str2double(match{2});
end

function events = parse_events(data, fs, name, task)
    cue = data(8,:);
    targetAction = data(9,:);
    cueStart = find(cue ~= 0 & [true, cue(1:end-1) == 0]);
    targetStart = find(abs(targetAction) == 1 & ...
        [true, abs(targetAction(1:end-1)) ~= 1]);
    clickStart = find(abs(targetAction) == 2 & ...
        [true, abs(targetAction(1:end-1)) ~= 2]);
    n = numel(cueStart);
    if n ~= 100
        error('%s has %d cue events, expected 100.', name, n);
    end
    trial = (1:n)';
    cueIndex = cueStart(:);
    cueDir = cue(cueIndex)';
    cueTime = data(10,cueIndex)';
    cueOffsetIndex = nan(n,1);
    cueDurationSec = nan(n,1);
    targetIndex = nan(n,1); targetDir = nan(n,1); targetTime = nan(n,1);
    clickIndex = nan(n,1); clickDir = nan(n,1); clickTime = nan(n,1);
    for i = 1:n
        nextCue = size(data,2) + 1;
        if i < n
            nextCue = cueStart(i+1);
        end
        offset = find(cue(cueStart(i):nextCue-1) == 0,1,'first');
        if ~isempty(offset)
            cueOffsetIndex(i) = cueStart(i)+offset-1;
            cueDurationSec(i) = data(10,cueOffsetIndex(i))-cueTime(i);
        end
        targetCandidate = targetStart(targetStart > cueStart(i) & targetStart < nextCue);
        if ~isempty(targetCandidate)
            p = targetCandidate(1);
            targetIndex(i) = p;
            targetDir(i) = sign(targetAction(p));
            targetTime(i) = data(10,p);
            clickCandidate = clickStart(clickStart > p & clickStart < nextCue);
            if ~isempty(clickCandidate)
                q = clickCandidate(1);
                clickIndex(i) = q;
                clickDir(i) = sign(targetAction(q));
                clickTime(i) = data(10,q);
            end
        end
    end
    if any(isnan(targetIndex))
        warning('%s has %d cues without a matched target.', name, sum(isnan(targetIndex)));
    end
    record = repmat(name,n,1);
    task = repmat(task,n,1);
    cueToTargetSec = targetTime - cueTime;
    targetToClickSec = clickTime - targetTime;
    events = table(record,task,trial,cueIndex,cueTime,cueDir, ...
        cueOffsetIndex,cueDurationSec, ...
        targetIndex,targetTime,targetDir,clickIndex,clickTime,clickDir, ...
        cueToTargetSec,targetToClickSec);
    if abs(median(cueToTargetSec,'omitnan') - 2.21) > 0.15
        warning('%s has unexpected cue-to-target timing.', name);
    end
    if task(1) == 1 && any(~isnan(clickIndex))
        warning('%s contains explicit click pulses; inspect Task-1 labeling.', name);
    end
    if task(1) == 2 && sum(~isnan(clickIndex)) ~= n
        warning('%s has missing explicit click pulses.', name);
    end
    if abs(median(diff(data(10,:))) - 1/fs) > 0.001
        warning('%s timestamp spacing differs from nominal sample rate.', name);
    end
end

function [epochs, clipped, peakToPeak, maxJump, inBounds] = ...
        make_epochs(data, cueIndex, fs, cfg, time)
    raw = data(1:3,:);
    [bh,ah] = butter(2, cfg.highPassHz/(fs/2), 'high');
    [bl,al] = butter(4, cfg.lowPassHz/(fs/2), 'low');
    filtered = filtfilt(bl,al,filtfilt(bh,ah,raw') )';
    n = numel(cueIndex);
    nTime = numel(time);
    epochs = nan(n,nTime,3);
    clipped = false(n,1); peakToPeak = nan(n,1);
    maxJump = nan(n,1); inBounds = false(n,1);
    leftCount = round(cfg.preSec*fs);
    rightCount = round(cfg.postSec*fs);
    baseline = time < 0;
    for i = 1:n
        lo = cueIndex(i)-leftCount;
        hi = cueIndex(i)+rightCount;
        if lo < 1 || hi > size(data,2)
            continue
        end
        inBounds(i) = true;
        rawEpoch = raw(:,lo:hi);
        clipped(i) = any(abs(rawEpoch(:)) >= cfg.clipLevel);
        peakToPeak(i) = max(max(rawEpoch,[],2)-min(rawEpoch,[],2));
        maxJump(i) = max(abs(diff(rawEpoch,1,2)),[],'all');
        segment = filtered(:,lo:hi);
        segment = segment - mean(segment(:,baseline),2);
        epochs(i,:,:) = reshape(segment',1,nTime,3);
    end
end

function [accepted, highPeak, highJump, peakThreshold, jumpThreshold] = ...
        mark_artifacts(clipped, peakToPeak, maxJump, inBounds, cfg)
    reference = inBounds & ~clipped;
    peakThreshold = robust_upper(peakToPeak(reference), cfg.artifactMadMultiplier);
    jumpThreshold = robust_upper(maxJump(reference), cfg.artifactMadMultiplier);
    highPeak = peakToPeak > peakThreshold;
    highJump = maxJump > jumpThreshold;
    accepted = inBounds & ~clipped & ~highPeak & ~highJump;
end

function threshold = robust_upper(values, multiplier)
    m = median(values);
    scale = 1.4826 * median(abs(values-m));
    threshold = m + multiplier * max(scale,eps(max(abs(m),1)));
end

function [curveTable, fitted] = estimate_curves(epochs, cueDir, accepted, ...
        time, name, task, cfg)
    curveTable = table();
    fitted = nan(numel(time),3,2);
    for sideIdx = 1:2
        direction = 2*sideIdx-3;
        subset = epochs(accepted & cueDir == direction,:,:);
        for e = 1:3
            x = reshape(subset(:,:,e),size(subset,1),numel(time));
            ordinary = mean(x,1)';
            robust = huber_location(x,cfg)';
            smooth = smooth_curve(robust,cfg.smoothLambda);
            fitted(:,e,sideIdx) = smooth;
            n = numel(time);
            row = table(repmat(name,n,1),repmat(task,n,1), ...
                repmat(direction,n,1),repmat(cfg.electrodeNames(e),n,1), ...
                time(:),ordinary,robust,smooth, ...
                'VariableNames', {'record','task','cueDir','electrode', ...
                'timeSec','ordinaryMean','huberMean','fittedCurve'});
            curveTable = [curveTable;row]; %#ok<AGROW>
        end
    end
end

function location = huber_location(x,cfg)
    location = median(x,1);
    scale = 1.4826*median(abs(x-location),1);
    scale = max(scale,1e-9);
    for iteration = 1:cfg.huberIterations
        residual = abs(x-location);
        weight = min(1,cfg.huberConstant*scale./max(residual,1e-12));
        updated = sum(weight.*x,1)./sum(weight,1);
        if max(abs(updated-location),[],'all') < 1e-8
            location = updated;
            break
        end
        location = updated;
    end
end

function smooth = smooth_curve(values, lambda)
    n = numel(values);
    secondDifference = spdiags([ones(n-2,1),-2*ones(n-2,1), ...
        ones(n-2,1)],[0,1,2],n-2,n);
    smooth = (speye(n)+lambda*(secondDifference'*secondDifference))\values(:);
end

function result = window_stats(epochs,cueDir,accepted,time,name,task,cfg)
    result = table();
    for windowIdx = 1:size(cfg.windows,1)
        mask = time >= cfg.windows(windowIdx,1) & ...
               time < cfg.windows(windowIdx,2);
        perTrial = reshape(mean(epochs(:,mask,:),2),size(epochs,1),3);
        perTrial(:,4) = perTrial(:,2)-perTrial(:,3);
        electrodes = [cfg.electrodeNames,"F3_minus_F4"];
        for e = 1:4
            left = perTrial(accepted & cueDir == -1,e);
            right = perTrial(accepted & cueDir == 1,e);
            for condition = 1:3
                if condition == 1
                    estimate = huber_location(left,cfg);
                    boot = bootstrap_location(left,cfg);
                    label = "left";
                elseif condition == 2
                    estimate = huber_location(right,cfg);
                    boot = bootstrap_location(right,cfg);
                    label = "right";
                else
                    estimate = huber_location(right,cfg)-huber_location(left,cfg);
                    boot = bootstrap_location(right,cfg)-bootstrap_location(left,cfg);
                    label = "right_minus_left";
                end
                ci = prctile(boot,[2.5,97.5]);
                row = table(name,task,cfg.windowNames(windowIdx), ...
                    electrodes(e),label,estimate,ci(1),ci(2), ...
                    numel(left),numel(right), ...
                    'VariableNames', {'record','task','window','electrode', ...
                    'condition','estimate','ciLow','ciHigh','nLeft','nRight'});
                result = [result;row]; %#ok<AGROW>
            end
        end
    end
end

function boot = bootstrap_location(x,cfg)
    n = numel(x);
    boot = zeros(cfg.bootstrapIterations,1);
    for iteration = 1:cfg.bootstrapIterations
        sample = x(randi(n,n,1));
        boot(iteration) = huber_location(sample,cfg);
    end
end

function result = split_half(epochs,cueDir,accepted,time,name,task,cfg)
    result = table();
    mask = time >= 0 & time <= 0.8;
    for sideIdx = 1:2
        direction = 2*sideIdx-3;
        ids = find(accepted & cueDir == direction);
        values = zeros(cfg.splitIterations,1);
        half = floor(numel(ids)/2);
        for iteration = 1:cfg.splitIterations
            ordered = ids(randperm(numel(ids)));
            a = squeeze(mean(epochs(ordered(1:half),mask,:),1));
            b = squeeze(mean(epochs(ordered(half+1:end),mask,:),1));
            c = corrcoef(a(:),b(:));
            values(iteration) = c(1,2);
        end
        ci = prctile(values,[2.5,97.5]);
        row = table(name,task,direction,numel(ids),mean(values), ...
            ci(1),ci(2),'VariableNames',{'record','task','cueDir', ...
            'nTrials','meanCorrelation','ciLow','ciHigh'});
        result = [result;row]; %#ok<AGROW>
    end
end

function result = contrast_split_half(epochs,cueDir,accepted,time,name,task,cfg)
% Test the repeatability of the RIGHT-LEFT waveform, not the common ERP.
    mask = time >= 0 & time <= cfg.postSec;
    leftIds = find(accepted & cueDir == -1);
    rightIds = find(accepted & cueDir == 1);
    leftHalf = floor(numel(leftIds)/2);
    rightHalf = floor(numel(rightIds)/2);
    values = nan(cfg.splitIterations,1);
    for iteration = 1:cfg.splitIterations
        leftOrder = leftIds(randperm(numel(leftIds)));
        rightOrder = rightIds(randperm(numel(rightIds)));
        leftA = squeeze(mean(epochs(leftOrder(1:leftHalf),mask,:),1));
        leftB = squeeze(mean(epochs(leftOrder(leftHalf+1:end),mask,:),1));
        rightA = squeeze(mean(epochs(rightOrder(1:rightHalf),mask,:),1));
        rightB = squeeze(mean(epochs(rightOrder(rightHalf+1:end),mask,:),1));
        a = rightA-leftA;
        b = rightB-leftB;
        correlation = corrcoef(a(:),b(:));
        values(iteration) = correlation(1,2);
    end
    ci = prctile(values,[2.5,97.5]);
    result = table(name,task,numel(leftIds),numel(rightIds), ...
        mean(values,'omitnan'),ci(1),ci(2), ...
        'VariableNames',{'record','task','nLeft','nRight', ...
        'meanCorrelation','ciLow','ciHigh'});
end

function result = filter_sensitivity(data,cueIndex,cueDir,accepted,fs,name,task,cfg)
% Fixed accepted trials; inspect whether slow wave and shape contrast
% survive plausible high-pass choices. This does not select a filter by p.
    result = table();
    preCount = round(cfg.diagnosticPreSec*fs);
    postCount = round(cfg.postSec*fs);
    time = (-preCount:postCount)/fs;
    baseline = time >= -cfg.preSec & time < 0;
    preEarly = time >= -0.7 & time < -0.5;
    preLate = time >= -0.2 & time < -0.05;
    mid = time >= 0.25 & time < 0.5;
    late = time >= 0.6 & time <= 0.8;
    [bl,al] = butter(4,cfg.lowPassHz/(fs/2),'low');
    for highPass = cfg.diagnosticHighPassHz
        [bh,ah] = butter(2,highPass/(fs/2),'high');
        filtered = filtfilt(bl,al,filtfilt(bh,ah,data(1:3,:)'))';
        wave = nan(numel(cueIndex),numel(time),3);
        valid = accepted;
        for i = 1:numel(cueIndex)
            lo = cueIndex(i)-preCount;
            hi = cueIndex(i)+postCount;
            if ~valid(i)
                continue
            end
            if lo < 1 || hi > size(data,2) || ...
                    any(abs(data(1:3,lo:hi)) >= cfg.clipLevel,'all')
                valid(i) = false;
                continue
            end
            segment = filtered(:,lo:hi);
            segment = segment-mean(segment(:,baseline),2);
            wave(i,:,:) = reshape(segment',1,numel(time),3);
        end
        for e = 1:3
            curves = nan(3,numel(time));
            counts = zeros(3,1);
            for side = 1:2
                direction = 2*side-3;
                subset = wave(valid & cueDir == direction,:,e);
                subset = reshape(subset,size(subset,1),numel(time));
                counts(side) = size(subset,1);
                curves(side,:) = huber_location(subset,cfg);
            end
            curves(3,:) = curves(2,:)-curves(1,:);
            counts(3) = min(counts(1:2));
            labels = ["left","right","right_minus_left"];
            for condition = 1:3
                curve = curves(condition,:);
                row = table(name,task,highPass,cfg.electrodeNames(e), ...
                    labels(condition),counts(condition), ...
                    mean(curve(preEarly)),mean(curve(preLate)), ...
                    mean(curve(mid)),mean(curve(late)), ...
                    mean(curve(late))-mean(curve(mid)), ...
                    'VariableNames',{'record','task','highPassHz', ...
                    'electrode','condition','nTrials','preEarlyMean', ...
                    'preLateMean','midMean','lateMean','lateMinusMid'});
                result = [result;row]; %#ok<AGROW>
            end
        end
    end
end

function [summary,commonSummary] = late_wave_diagnostic( ...
        epochs,cueDir,accepted,time,name,task,cfg)
% Exploratory characterization of the large wave after 0.55 s. A common
% scalp waveform alone cannot distinguish neural activity from eye motion.
    summary = table();
    commonSummary = table();
    mid = time >= 0.25 & time < 0.50;
    late = time >= 0.55 & time <= 0.80;
    correlationWindow = time >= 0.45 & time <= 0.80;
    transformed = cat(3,epochs,mean(epochs,3), ...
        epochs(:,:,2)-epochs(:,:,3));
    labels = [cfg.electrodeNames,"common_FzF3F4","F3_minus_F4"];
    for side = 1:2
        direction = 2*side-3;
        selected = accepted & cueDir == direction;
        n = sum(selected);
        curves = nan(3,numel(time));
        lateEstimates = nan(5,1);
        for e = 1:5
            x = reshape(transformed(selected,:,e),n,numel(time));
            middleTrial = mean(x(:,mid),2);
            lateTrial = mean(x(:,late),2);
            changeTrial = lateTrial-middleTrial;
            estimate = huber_location(changeTrial,cfg);
            ci = prctile(bootstrap_location(changeTrial,cfg),[2.5,97.5]);
            curve = huber_location(x,cfg);
            lateTime = time(late);
            [peak,peakIndex] = max(curve(late));
            lateEstimates(e) = huber_location(lateTrial,cfg);
            row = table(name,task,direction,labels(e),n, ...
                huber_location(middleTrial,cfg),lateEstimates(e), ...
                estimate,ci(1),ci(2),peak,lateTime(peakIndex), ...
                'VariableNames',{'record','task','cueDir','channel', ...
                'nTrials','midMean','lateMean','lateMinusMid', ...
                'ciLow','ciHigh','latePeak','peakTimeSec'});
            summary = [summary;row]; %#ok<AGROW>
            if e <= 3
                curves(e,:) = curve;
            end
        end
        r = corrcoef(curves(:,correlationWindow)');
        commonRow = table(name,task,direction,n,r(1,2),r(1,3), ...
            r(2,3),lateEstimates(4),lateEstimates(5), ...
            'VariableNames',{'record','task','cueDir','nTrials', ...
            'corrFzF3','corrFzF4','corrF3F4','commonLateMean', ...
            'lateralLateMean'});
        commonSummary = [commonSummary;commonRow]; %#ok<AGROW>
    end
end

function features = extract_features(epochs,time,cfg)
    n = size(epochs,1);
    features = nan(n,8);
    for windowIdx = 1:2
        mask = time >= cfg.windows(windowIdx,1) & ...
               time < cfg.windows(windowIdx,2);
        amplitudes = reshape(mean(epochs(:,mask,:),2),n,3);
        columns = (windowIdx-1)*4+(1:4);
        features(:,columns) = [amplitudes, ...
            amplitudes(:,2)-amplitudes(:,3)];
    end
end

function features = extract_spatial_features(epochs,time,cfg)
% Two independent frontal contrasts remove the three-channel common mode.
% This is an analysis montage, not proof that common activity is artifact.
    n = size(epochs,1);
    features = nan(n,4);
    for windowIdx = 1:2
        mask = time >= cfg.windows(windowIdx,1) & ...
               time < cfg.windows(windowIdx,2);
        amplitude = reshape(mean(epochs(:,mask,:),2),n,3);
        columns = (windowIdx-1)*2+(1:2);
        features(:,columns) = [amplitude(:,1)- ...
            0.5*(amplitude(:,2)+amplitude(:,3)), ...
            amplitude(:,2)-amplitude(:,3)];
    end
end

function result = spatial_mode_curves(fitted,time,name,task)
    result = table();
    for side = 1:2
        direction = 2*side-3;
        fz = fitted(:,1,side);
        f3 = fitted(:,2,side);
        f4 = fitted(:,3,side);
        n = numel(time);
        row = table(repmat(name,n,1),repmat(task,n,1), ...
            repmat(direction,n,1),time(:), ...
            (fz+f3+f4)/3,fz-(f3+f4)/2,f3-f4, ...
            'VariableNames',{'record','task','cueDir','timeSec', ...
            'commonMode','midlineContrast','lateralContrast'});
        result = [result;row]; %#ok<AGROW>
    end
end

function result = evaluate_cross_record(records,cfg,featureField)
    result = table();
    for task = 1:2
        ids = find([records.task] == task);
        if numel(ids) ~= 2
            error('Expected exactly two recordings for Task %d.', task);
        end
        for direction = 1:2
            train = records(ids(direction));
            test = records(ids(3-direction));
            trainX = train.(featureField)(train.accepted,:);
            trainY = train.cue(train.accepted);
            testX = test.(featureField)(test.accepted,:);
            testY = test.cue(test.accepted);
            [scores,prediction] = shrinkage_lda(trainX,trainY,testX);
            observed = balanced_accuracy(testY,prediction);
            auc = roc_auc(scores,testY);
            permuted = zeros(cfg.permutationIterations,1);
            for iteration = 1:cfg.permutationIterations
                labels = testY(randperm(numel(testY)));
                permuted(iteration) = balanced_accuracy(labels,prediction);
            end
            p = (1+sum(permuted >= observed))/(1+cfg.permutationIterations);
            row = table(string(featureField),task,train.name,test.name, ...
                numel(trainY),numel(testY), ...
                observed,auc,p, ...
                'VariableNames',{'featureSet','task','trainRecord', ...
                'testRecord','nTrain', ...
                'nTest','balancedAccuracy','auc','permutationP'});
            result = [result;row]; %#ok<AGROW>
        end
    end
    result.holmP = holm_adjust(result.permutationP);
end

function [summary, folds] = evaluate_within_record(records,cfg,featureField)
% Chronological held-out blocks. Each test trial receives one prediction.
% Circularly shifting the original cue sequence supplies a null that keeps
% its run structure and temporal autocorrelation while breaking alignment.
    summary = table();
    folds = table();
    for recordIdx = 1:numel(records)
        rec = records(recordIdx);
        trialIds = find(rec.accepted);
        X = rec.(featureField)(trialIds,:);
        y = rec.cue(trialIds);
        n = numel(y);
        foldId = min(cfg.withinRecordFolds, ...
            ceil((1:n)'*cfg.withinRecordFolds/n));
        [scores,prediction] = crossvalidated_lda(X,y,foldId);
        observed = balanced_accuracy(y,prediction);
        auc = roc_auc(scores,y);
        shifted = nan(numel(rec.cue)-1,1);
        for shift = 1:numel(shifted)
            shiftedCue = circshift(rec.cue,shift);
            shiftedY = shiftedCue(trialIds);
            [~,shiftedPrediction] = crossvalidated_lda(X,shiftedY,foldId);
            shifted(shift) = balanced_accuracy(shiftedY,shiftedPrediction);
        end
        p = (1+sum(shifted >= observed))/(1+numel(shifted));
        row = table(string(featureField),rec.name,rec.task,n, ...
            observed,auc,p, ...
            mean(shifted),prctile(shifted,95), ...
            'VariableNames',{'featureSet','record','task', ...
            'nTrials','balancedAccuracy', ...
            'auc','circularShiftP','nullMeanAccuracy','null95Accuracy'});
        summary = [summary;row]; %#ok<AGROW>
        for k = 1:cfg.withinRecordFolds
            test = foldId == k;
            if any(y(test) == -1) && any(y(test) == 1)
                foldAccuracy = balanced_accuracy(y(test),prediction(test));
                foldAuc = roc_auc(scores(test),y(test));
            else
                foldAccuracy = NaN;
                foldAuc = NaN;
            end
            foldRow = table(string(featureField),rec.name,rec.task,k, ...
                min(trialIds(test)), ...
                max(trialIds(test)),sum(test),sum(y(test)==-1), ...
                sum(y(test)==1),foldAccuracy,foldAuc, ...
                'VariableNames',{'featureSet','record','task', ...
                'fold','firstTrial', ...
                'lastTrial','nTest','nLeft','nRight', ...
                'balancedAccuracy','auc'});
            folds = [folds;foldRow]; %#ok<AGROW>
        end
    end
    summary.holmP = holm_adjust(summary.circularShiftP);
end

function [scores,prediction] = crossvalidated_lda(X,y,foldId)
    scores = nan(numel(y),1);
    prediction = nan(numel(y),1);
    for fold = unique(foldId)'
        test = foldId == fold;
        train = ~test;
        if sum(y(train)==-1) < 2 || sum(y(train)==1) < 2
            error('A training fold lacks one cue class.');
        end
        [scores(test),prediction(test)] = ...
            shrinkage_lda(X(train,:),y(train),X(test,:));
    end
end

function [scores,prediction] = shrinkage_lda(trainX,trainY,testX)
    center = mean(trainX,1);
    spread = std(trainX,0,1);
    spread(spread < 1e-10) = 1;
    trainX = (trainX-center)./spread;
    testX = (testX-center)./spread;
    left = trainX(trainY == -1,:);
    right = trainX(trainY == 1,:);
    muLeft = mean(left,1)';
    muRight = mean(right,1)';
    pooled = ((size(left,1)-1)*cov(left)+(size(right,1)-1)*cov(right)) / ...
             (size(left,1)+size(right,1)-2);
    diagonal = diag(diag(pooled));
    ridge = max(trace(pooled)/size(pooled,1),1e-6)*1e-3;
    regularized = 0.7*pooled+0.3*diagonal+ridge*eye(size(pooled));
    w = regularized\(muRight-muLeft);
    intercept = -0.5*(muRight+muLeft)'*w;
    scores = testX*w+intercept;
    prediction = ones(numel(scores),1);
    prediction(scores < 0) = -1;
end

function value = balanced_accuracy(actual,prediction)
    value = 0.5*(mean(prediction(actual == -1) == -1) + ...
                 mean(prediction(actual == 1) == 1));
end

function auc = roc_auc(scores,actual)
    [~,order] = sort(scores);
    ranks = zeros(numel(scores),1);
    ranks(order) = 1:numel(scores);
    [uniqueScores,~,group] = unique(scores);
    if numel(uniqueScores) < numel(scores)
        for i = 1:numel(uniqueScores)
            tied = group == i;
            ranks(tied) = mean(ranks(tied));
        end
    end
    nPos = sum(actual == 1);
    nNeg = sum(actual == -1);
    auc = (sum(ranks(actual == 1))-nPos*(nPos+1)/2)/(nPos*nNeg);
end

function adjusted = holm_adjust(p)
    [sorted,order] = sort(p);
    m = numel(p);
    corrected = min(1,(m:-1:1)'.*sorted);
    corrected = cummax(corrected);
    adjusted = zeros(size(p));
    adjusted(order) = corrected;
end

function draw_erp(time,fitted,name,outputFile)
    figureHandle = figure('Visible','off','Color','w', ...
        'Position',[100,100,1000,760]);
    electrodeNames = ["Fz","F3","F4"];
    for e = 1:3
        subplot(3,1,e);
        plot(time,fitted(:,e,1),'LineWidth',1.6,'Color',[0.15 0.38 0.78]);
        hold on;
        plot(time,fitted(:,e,2),'LineWidth',1.6,'Color',[0.80 0.26 0.20]);
        xline(0,':k');
        yline(0,':','Color',[0.5 0.5 0.5]);
        xlim([min(time),max(time)]);
        grid on;
        ylabel('Amplitude (source units)');
        title(electrodeNames(e));
        if e == 1
            legend('Left cue','Right cue','Location','best');
        end
        if e == 3
            xlabel('Time from cue (s)');
        end
    end
    sgtitle(strrep(name,'_','\_') + " raw-channel robust ERP");
    exportgraphics(figureHandle,outputFile,'Resolution',180);
    close(figureHandle);
end
