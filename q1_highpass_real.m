function q1_highpass_real(sourceDir, outputDir, fixedEventsCsv)
%Q1_HIGHPASS_REAL Fixed-trial, four-cutoff EEG high-pass sensitivity.
%   q1_highpass_real(sourceDir, outputDir, fixedEventsCsv)
%   uses the *existing* accepted flags and cue labels from
%   results/events_and_quality.csv. It never repeats artifact rejection.
%   Units are the source MAT units throughout. Requires Signal Processing
%   Toolbox (butter/filtfilt) and Statistics and Machine Learning Toolbox
%   (prctile). This is a standalone diagnostic, not a change to run_q1.m.

    arguments
        sourceDir (1,1) string
        outputDir (1,1) string
        fixedEventsCsv (1,1) string
    end
    if ~isfolder(sourceDir)
        error('Source directory does not exist: %s', sourceDir);
    end
    if ~isfile(fixedEventsCsv)
        error('Fixed event table does not exist: %s', fixedEventsCsv);
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    cfg = local_config();
    rng(cfg.seed, 'twister');
    events = readtable(fixedEventsCsv);
    required = ["record","task","trial","cueIndex","cueDir","accepted"];
    if ~all(ismember(required, string(events.Properties.VariableNames)))
        error('The fixed event table is missing a required column.');
    end
    if height(events) ~= 400 || sum(logical(events.accepted)) ~= 350
        error('Expected the canonical 400 events and exactly 350 accepted trials.');
    end

    files = ["VisualCogA_Task-1.mat","VisualCogA_Task-2.mat", ...
             "VisualCogB_Task-1.mat","VisualCogB_Task-2.mat"];
    trialBalance = table();
    waveformRows = table();
    metricRows = struct([]);
    for fileNo = 1:numel(files)
        fileName = files(fileNo);
        token = regexp(char(fileName),'VisualCog([AB])_Task-([12])\.mat$', ...
                       'tokens','once');
        record = string(token{1});
        task = str2double(token{2});
        canonicalName = "VisualCog" + record + "_Task-" + string(task);
        take = string(events.record) == canonicalName & events.task == task;
        ev = events(take,:);
        ev = sortrows(ev,'trial');
        if height(ev) ~= 100 || ~isequal(ev.trial(:),(1:100)')
            error('%s does not have precisely 100 canonical ordered events.',fileName);
        end
        accepted = logical(ev.accepted);
        direction = double(ev.cueDir);
        if any(~ismember(direction,[-1,1]))
            error('%s has unexpected cue directions.',fileName);
        end
        for condition = [-1,1]
            original = sum(direction == condition);
            kept = sum(accepted & direction == condition);
            trialBalance = [trialBalance; table(record,task,condition, ...
                original,kept,kept/original,'VariableNames', ...
                {'record','task','condition','original_trials', ...
                 'retained_trials','retention_rate'})]; %#ok<AGROW>
        end

        inputFile = fullfile(sourceDir,fileName);
        if ~isfile(inputFile)
            error('Missing source MAT: %s',inputFile);
        end
        source = load(inputFile,'data','SampleRate');
        if ~isfield(source,'data') || ~isfield(source,'SampleRate') || ...
                size(source.data,1) < 8 || ...
                ~isscalar(source.SampleRate)
            error('Unexpected MAT layout: %s',inputFile);
        end
        fs = double(source.SampleRate);
        cueIndex = double(ev.cueIndex);
        if any(cueIndex < 1 | cueIndex > size(source.data,2) | ...
                cueIndex ~= round(cueIndex)) || ...
                any(source.data(8,cueIndex) ~= direction')
            error('Canonical cue index/direction differs from %s.',inputFile);
        end
        if fileNo > 1 && fs ~= referenceFs
            error('Sampling rates differ between the four source MAT files.');
        end
        referenceFs = fs;
        preCount = round(cfg.diagnosticPreSec*fs);
        formalPreCount = round(cfg.preSec*fs);
        postCount = round(cfg.postSec*fs);
        time = (-preCount:postCount)/fs;
        formalMask = time >= -formalPreCount/fs;
        formalTime = time(formalMask);
        baselineMask = time >= -formalPreCount/fs & time < 0;
        diagnosticMask = time >= -cfg.diagnosticPreSec & ...
                         time < -cfg.preSec;
        earlyMask = time >= -0.7 & time < -0.5;
        latePreMask = time >= -0.4 & time < -0.2;
        postMask = formalTime >= 0 & formalTime <= cfg.postSec;
        windowMask = formalTime >= 0.25 & formalTime < 0.5;
        if ~any(windowMask) || ~any(diagnosticMask)
            error('Sampling rate is incompatible with diagnostic windows.');
        end
        acceptedIndex = find(accepted);
        if any(cueIndex(acceptedIndex)-preCount < 1 | ...
               cueIndex(acceptedIndex)+postCount > size(source.data,2))
            error(['An accepted canonical trial lacks the -0.8 to 0.8 s ' ...
                   'diagnostic window. No extra trial rejection is allowed.']);
        end
        acceptedDir = direction(acceptedIndex);
        leftIndex = find(acceptedDir == -1);
        rightIndex = find(acceptedDir == 1);
        [bootLeft,bootRight,splitLeft,splitRight] = ...
            fixed_resamples(numel(leftIndex),numel(rightIndex),cfg);
        smoother = smoothing_matrix(numel(formalTime),cfg.smoothLambda);
        [bl,al] = butter(4,cfg.lowPassHz/(fs/2),'low');
        referenceLeft = nan(numel(formalTime),3);
        referenceRight = nan(numel(formalTime),3);
        referenceDelta = nan(numel(formalTime),3);

        for cutoff = cfg.highPassHz
            [bh,ah] = butter(2,cutoff/(fs/2),'high');
            filtered = filtfilt(bl,al, ...
                filtfilt(bh,ah,double(source.data(1:3,:))'))';
            epoch = nan(numel(acceptedIndex),numel(time),3);
            for k = 1:numel(acceptedIndex)
                at = cueIndex(acceptedIndex(k));
                segment = filtered(:,at-preCount:at+postCount);
                segment = segment - mean(segment(:,baselineMask),2);
                epoch(k,:,:) = reshape(segment',1,numel(time),3);
            end
            if any(~isfinite(epoch),'all')
                error('Nonfinite filtered samples in %s at %.2f Hz.',fileName,cutoff);
            end

            for channelNo = 1:3
                channel = cfg.channelNames(channelNo);
                leftTrials = reshape(epoch(leftIndex,:,channelNo), ...
                    numel(leftIndex),numel(time));
                rightTrials = reshape(epoch(rightIndex,:,channelNo), ...
                    numel(rightIndex),numel(time));
                leftHuber = huber_location(leftTrials,cfg);
                rightHuber = huber_location(rightTrials,cfg);
                deltaHuber = rightHuber-leftHuber;
                leftFit = smoother\leftHuber(formalMask)';
                rightFit = smoother\rightHuber(formalMask)';
                deltaFit = rightFit-leftFit;
                if cutoff == cfg.highPassHz(1)
                    referenceLeft(:,channelNo) = leftFit;
                    referenceRight(:,channelNo) = rightFit;
                    referenceDelta(:,channelNo) = deltaFit;
                end
                leftAll = nan(numel(time),1);
                rightAll = leftAll;
                deltaAll = leftAll;
                leftAll(formalMask) = leftFit;
                rightAll(formalMask) = rightFit;
                deltaAll(formalMask) = deltaFit;
                nTime = numel(time);
                wave = table(repmat(record,nTime,1),repmat(task,nTime,1), ...
                    repmat(channel,nTime,1),repmat(cutoff,nTime,1), ...
                    time(:),leftAll,rightAll,deltaAll,leftHuber(:), ...
                    rightHuber(:),deltaHuber(:), ...
                    'VariableNames',{'record','task','channel','highpass', ...
                    'timeSec','left','right','delta','left_huber', ...
                    'right_huber','delta_huber'});
                waveformRows = [waveformRows;wave]; %#ok<AGROW>

                baselineSamples = [leftTrials(:,baselineMask); ...
                                   rightTrials(:,baselineMask)];
                baselineSamples = baselineSamples(:);
                baselineRms = sqrt(mean(baselineSamples.^2));
                baselineMad = 1.4826*median(abs( ...
                    baselineSamples-median(baselineSamples)));
                baselineVariance = var(baselineSamples,0);
                leftWindowTrials = mean(leftTrials(:, ...
                    time >= 0.25 & time < 0.5),2);
                rightWindowTrials = mean(rightTrials(:, ...
                    time >= 0.25 & time < 0.5),2);
                trialVariance = ((numel(leftIndex)-1)*var(leftWindowTrials,0) ...
                    +(numel(rightIndex)-1)*var(rightWindowTrials,0)) / ...
                    (numel(leftIndex)+numel(rightIndex)-2);

                [rLeft,rRight,rDelta,splitMeanAbsDiff] = split_reliability( ...
                    leftTrials(:,formalMask),rightTrials(:,formalMask), ...
                    splitLeft,splitRight,smoother,postMask,windowMask,cfg);
                [bootLow,bootHigh,bootWaveSd] = bootstrap_waveform( ...
                    leftTrials(:,formalMask),rightTrials(:,formalMask), ...
                    bootLeft,bootRight,smoother,postMask,windowMask,cfg);

                leftPreMean = huber_location(mean(leftTrials(:,diagnosticMask),2),cfg);
                rightPreMean = huber_location(mean(rightTrials(:,diagnosticMask),2),cfg);
                leftSlope = robust_slope(leftTrials(:,diagnosticMask), ...
                    time(diagnosticMask),cfg);
                rightSlope = robust_slope(rightTrials(:,diagnosticMask), ...
                    time(diagnosticMask),cfg);
                leftEarly = huber_location(mean(leftTrials(:,earlyMask),2),cfg);
                rightEarly = huber_location(mean(rightTrials(:,earlyMask),2),cfg);
                leftLate = huber_location(mean(leftTrials(:,latePreMask),2),cfg);
                rightLate = huber_location(mean(rightTrials(:,latePreMask),2),cfg);
                leftStats = curve_measure(leftFit,formalTime,windowMask);
                rightStats = curve_measure(rightFit,formalTime,windowMask);
                deltaStats = curve_measure(deltaFit,formalTime,windowMask);

                row = struct();
                row.record = record;
                row.task = task;
                row.channel = channel;
                row.highpass = cutoff;
                row.n_left = numel(leftIndex);
                row.n_right = numel(rightIndex);
                row.baseline_rms = baselineRms;
                row.baseline_mad = baselineMad;
                row.baseline_variance = baselineVariance;
                row.trial_to_trial_variance = trialVariance;
                row.split_half_r = mean(rDelta,'omitnan');
                row.split_half_left_r = mean(rLeft,'omitnan');
                row.split_half_right_r = mean(rRight,'omitnan');
                row.split_half_mean_abs_diff = median(splitMeanAbsDiff,'omitnan');
                row.bootstrap_ci_low = bootLow;
                row.bootstrap_ci_high = bootHigh;
                row.bootstrap_width = bootHigh-bootLow;
                row.bootstrap_waveform_sd = bootWaveSd;
                row.left_mean_250_500 = leftStats.mean;
                row.right_mean_250_500 = rightStats.mean;
                row.delta_mean_250_500 = deltaStats.mean;
                row.left_peak = leftStats.peak;
                row.right_peak = rightStats.peak;
                row.delta_peak = deltaStats.peak;
                row.left_latency = leftStats.latency;
                row.right_latency = rightStats.latency;
                row.delta_latency = deltaStats.latency;
                row.left_auc = leftStats.auc;
                row.right_auc = rightStats.auc;
                row.delta_auc = deltaStats.auc;
                row.delta_min = min(deltaFit(windowMask));
                row.delta_positive_fraction = mean(deltaFit(windowMask)>0);
                row.baseline_mean_left = leftPreMean;
                row.baseline_mean_right = rightPreMean;
                row.baseline_difference = rightPreMean-leftPreMean;
                row.baseline_slope_left = leftSlope;
                row.baseline_slope_right = rightSlope;
                row.baseline_slope_difference = rightSlope-leftSlope;
                row.baseline_early_difference = rightEarly-leftEarly;
                row.baseline_late_difference = rightLate-leftLate;
                row.left_waveform_r_vs_0p1 = safe_corr(leftFit(postMask), ...
                    referenceLeft(postMask,channelNo));
                row.right_waveform_r_vs_0p1 = safe_corr(rightFit(postMask), ...
                    referenceRight(postMask,channelNo));
                row.delta_waveform_r_vs_0p1 = safe_corr(deltaFit(postMask), ...
                    referenceDelta(postMask,channelNo));
                row.delta_waveform_rmse_vs_0p1 = sqrt(mean( ...
                    (deltaFit(postMask)-referenceDelta(postMask,channelNo)).^2));
                row.left_peak_change_vs_0p1 = leftStats.peak - ...
                    max(referenceLeft(windowMask,channelNo));
                row.right_peak_change_vs_0p1 = rightStats.peak - ...
                    max(referenceRight(windowMask,channelNo));
                row.delta_peak_change_vs_0p1 = deltaStats.peak - ...
                    max(referenceDelta(windowMask,channelNo));
                referenceLeftStats = curve_measure( ...
                    referenceLeft(:,channelNo),formalTime,windowMask);
                referenceRightStats = curve_measure( ...
                    referenceRight(:,channelNo),formalTime,windowMask);
                referenceDeltaStats = curve_measure( ...
                    referenceDelta(:,channelNo),formalTime,windowMask);
                row.left_latency_shift_ms_vs_0p1 = 1000* ...
                    (leftStats.latency-referenceLeftStats.latency);
                row.right_latency_shift_ms_vs_0p1 = 1000* ...
                    (rightStats.latency-referenceRightStats.latency);
                row.delta_latency_shift_ms_vs_0p1 = 1000* ...
                    (deltaStats.latency-referenceDeltaStats.latency);
                row.pre_stimulus_delta_rms = sqrt(mean( ...
                    deltaHuber(time < 0).^2));
                row.pre_stimulus_delta_zero_crossings = sum( ...
                    diff(sign(deltaHuber(time >= -0.2 & time < 0))) ~= 0);
                row.waveform_stability = NaN;
                metricRows = [metricRows;row]; %#ok<AGROW>
            end
            fprintf('%s %.2f Hz complete, fixed trials %d (%d Left, %d Right).\n', ...
                fileName,cutoff,numel(acceptedIndex),numel(leftIndex), ...
                numel(rightIndex));
        end
    end

    metrics = struct2table(metricRows);
    if height(metrics) ~= 48 || sum(trialBalance.retained_trials) ~= 350
        error('Output dimensionality does not match the fixed-trial design.');
    end
    metrics = fill_cross_cutoff_stability(metrics,waveformRows,cfg);
    summary = build_summary(metrics,cfg);
    writetable(trialBalance,fullfile(outputDir,'q1_trial_balance.csv'));
    writetable(metrics,fullfile(outputDir,'q1_highpass_sensitivity_real.csv'));
    writetable(waveformRows,fullfile(outputDir,'q1_highpass_waveforms.csv'));
    writetable(summary,fullfile(outputDir,'q1_highpass_summary.csv'));
    fprintf('Real-data high-pass sensitivity outputs saved to %s\n',outputDir);
end

function cfg = local_config()
    cfg.highPassHz = [0.10,0.20,0.25,0.50];
    cfg.lowPassHz = 30;
    cfg.preSec = 0.2;
    cfg.diagnosticPreSec = 0.8;
    cfg.postSec = 0.8;
    cfg.huberConstant = 1.345;
    cfg.huberIterations = 15;
    cfg.smoothLambda = 12;
    cfg.bootstrapIterations = 500;
    cfg.splitIterations = 100;
    cfg.seed = 20260922;
    cfg.channelNames = ["Fz","F3","F4"];
end

function [bootLeft,bootRight,splitLeft,splitRight] = ...
        fixed_resamples(nLeft,nRight,cfg)
% Common random numbers across all channels and HPs of one recording.
    bootLeft = randi(nLeft,nLeft,cfg.bootstrapIterations);
    bootRight = randi(nRight,nRight,cfg.bootstrapIterations);
    splitLeft = zeros(nLeft,cfg.splitIterations);
    splitRight = zeros(nRight,cfg.splitIterations);
    for iteration = 1:cfg.splitIterations
        splitLeft(:,iteration) = randperm(nLeft)';
        splitRight(:,iteration) = randperm(nRight)';
    end
end

function smoother = smoothing_matrix(n,lambda)
    difference = spdiags([ones(n-2,1),-2*ones(n-2,1), ...
        ones(n-2,1)],[0,1,2],n-2,n);
    smoother = decomposition(speye(n)+lambda*(difference'*difference),'chol');
end

function location = huber_location(x,cfg)
% Identical pointwise Huber implementation and parameters to run_q1.m.
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

function fitted = fit_trials(x,smoother,cfg)
    fitted = smoother\huber_location(x,cfg)';
end

function [rLeft,rRight,rDelta,meanAbsDiff] = split_reliability( ...
        left,right,splitLeft,splitRight,smoother,postMask,windowMask,cfg)
% Trial-independent halves; each half is estimated with the fixed Method B.
    iterations = cfg.splitIterations;
    rLeft = nan(iterations,1);
    rRight = nan(iterations,1);
    rDelta = nan(iterations,1);
    meanAbsDiff = nan(iterations,1);
    nLeftA = floor(size(left,1)/2);
    nRightA = floor(size(right,1)/2);
    for iteration = 1:iterations
        pL = splitLeft(:,iteration);
        pR = splitRight(:,iteration);
        lA = fit_trials(left(pL(1:nLeftA),:),smoother,cfg);
        lB = fit_trials(left(pL(nLeftA+1:end),:),smoother,cfg);
        rA = fit_trials(right(pR(1:nRightA),:),smoother,cfg);
        rB = fit_trials(right(pR(nRightA+1:end),:),smoother,cfg);
        dA = rA-lA;
        dB = rB-lB;
        rLeft(iteration) = safe_corr(lA(postMask),lB(postMask));
        rRight(iteration) = safe_corr(rA(postMask),rB(postMask));
        rDelta(iteration) = safe_corr(dA(postMask),dB(postMask));
        meanAbsDiff(iteration) = abs(mean(dA(windowMask))- ...
                                     mean(dB(windowMask)));
    end
end

function [ciLow,ciHigh,waveSd] = bootstrap_waveform( ...
        left,right,bootLeft,bootRight,smoother,postMask,windowMask,cfg)
% Stratified bootstrap of the *same* Huber-plus-smoothing contrast.
    n = cfg.bootstrapIterations;
    deltaBoot = zeros(n,size(left,2));
    windowBoot = zeros(n,1);
    for iteration = 1:n
        leftFit = fit_trials(left(bootLeft(:,iteration),:),smoother,cfg);
        rightFit = fit_trials(right(bootRight(:,iteration),:),smoother,cfg);
        delta = rightFit-leftFit;
        deltaBoot(iteration,:) = delta';
        windowBoot(iteration) = mean(delta(windowMask));
    end
    ci = prctile(windowBoot,[2.5,97.5]);
    ciLow = ci(1);
    ciHigh = ci(2);
    waveSd = mean(std(deltaBoot(:,postMask),0,1));
end

function slope = robust_slope(trials,t,cfg)
% Regression slope of each trial's extended pre-cue segment, then Huber.
    centered = t(:)-mean(t);
    slopes = (trials*centered)/(centered'*centered);
    slope = huber_location(slopes,cfg);
end

function result = curve_measure(curve,time,mask)
    x = curve(mask);
    t = time(mask);
    [peak,idx] = max(x);
    result.mean = mean(x);
    result.peak = peak;
    result.latency = t(idx);
    result.auc = trapz(t,x);
end

function r = safe_corr(a,b)
    a = a(:);
    b = b(:);
    if any(~isfinite(a)) || any(~isfinite(b)) || ...
            std(a) < 1e-12 || std(b) < 1e-12
        r = NaN;
    else
        c = corrcoef(a,b);
        r = c(1,2);
    end
end

function metrics = fill_cross_cutoff_stability(metrics,waves,cfg)
% Symmetric median shape correlation with the other three cutoffs.
% This metric does not presume that the 0.1-Hz waveform is ground truth.
    fields = ["left","right","delta"];
    for rowNo = 1:height(metrics)
        row = metrics(rowNo,:);
        own = waves(waves.record == row.record & ...
            waves.task == row.task & waves.channel == row.channel & ...
            waves.highpass == row.highpass & waves.timeSec >= 0,:);
        similarities = nan(1,3*(numel(cfg.highPassHz)-1));
        pos = 0;
        for otherHp = cfg.highPassHz
            if otherHp == row.highpass
                continue
            end
            other = waves(waves.record == row.record & ...
                waves.task == row.task & waves.channel == row.channel & ...
                waves.highpass == otherHp & waves.timeSec >= 0,:);
            if height(own) ~= height(other) || ...
                    any(own.timeSec ~= other.timeSec)
                error('Cross-cutoff waveform time grids differ.');
            end
            for field = fields
                pos = pos+1;
                similarities(pos) = safe_corr(own.(field),other.(field));
            end
        end
        metrics.waveform_stability(rowNo) = median(similarities,'omitnan');
    end
end

function result = build_summary(metrics,cfg)
% Descriptive trade-offs only. Change relative to 0.1 Hz is not ground truth.
    result = table();
    for hp = cfg.highPassHz
        x = metrics(metrics.highpass == hp,:);
        Noise = median(x.baseline_rms,'omitnan');
        ERP_Reliability = median((x.split_half_left_r+ ...
                                  x.split_half_right_r)/2,'omitnan');
        Waveform_Change_vs_0p1 = median( ...
            1-x.delta_waveform_r_vs_0p1,'omitnan');
        LR_Repeatability = median(x.split_half_r,'omitnan');
        Baseline_Contamination = median(abs(x.baseline_difference),'omitnan');
        Bootstrap_Width = median(x.bootstrap_width,'omitnan');
        Overall_Observation = "Descriptive trade-off; no weighted score or selected cutoff";
        row = table(hp,Noise,ERP_Reliability,Waveform_Change_vs_0p1, ...
            LR_Repeatability,Baseline_Contamination,Bootstrap_Width, ...
            Overall_Observation,'VariableNames',{'HP','Noise', ...
            'ERP_Reliability','Waveform_Change_vs_0p1', ...
            'LR_Repeatability','Baseline_Contamination', ...
            'Bootstrap_Width','Overall_Observation'});
        result = [result;row]; %#ok<AGROW>
    end
end
