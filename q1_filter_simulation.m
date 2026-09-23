function q1_filter_simulation(sourceDir, outputDir, fixedEventsCsv, nReplicates)
%Q1_FILTER_SIMULATION Semi-synthetic high-pass recovery check for Q1.
%   q1_filter_simulation(sourceDir, outputDir) uses the canonical 350-trial
%   mask in results/events_and_quality.csv. Optional fixedEventsCsv may
%   provide the same frozen event table. Optional nReplicates defaults to 12.
%
%   The 0.1- and 0.5-Hz training-derived templates are deliberately called
%   semi-synthetic references, not clean biological ground truth. Their
%   differing low-frequency content exposes reference-choice dependence.
%   Within each record, trials 1:60 provide template/noise calibration,
%   trials 61:64 form a temporal guard, and trials 65:100 provide held-out
%   recovery locations. All candidate cutoffs receive identical synthetic
%   signals and identical held-out events; filtering is continuous before
%   epoching. Source amplitudes are never relabeled as microvolts.

    arguments
        sourceDir (1,1) string
        outputDir (1,1) string
        fixedEventsCsv (1,1) string = ""
        nReplicates (1,1) double {mustBeInteger,mustBePositive} = 12
    end

    repoDir = string(fileparts(mfilename('fullpath')));
    if fixedEventsCsv == ""
        fixedEventsCsv = fullfile(repoDir, "results", ...
            "events_and_quality.csv");
    end
    if ~isfolder(sourceDir) || ~isfile(fixedEventsCsv)
        error('Missing source directory or frozen event CSV.');
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    cfg.fs = 256;
    cfg.lowPassHz = 30;
    cfg.candidateHp = [0.1 0.2 0.25 0.5];
    cfg.templateHp = [0.1 0.5];
    cfg.seed = 20260923;
    cfg.preSec = 0.2;
    cfg.postSec = 0.8;
    cfg.extPreSec = 1.5;
    cfg.extPostSec = 1.8;
    cfg.huberConstant = 1.345;
    cfg.huberIterations = 15;
    cfg.smoothLambda = 12;
    cfg.highQualityFraction = 0.8;
    cfg.lastTrainingTrial = 60;
    cfg.firstValidationTrial = 65;
    cfg.spikeProbability = 0.12;
    cfg.burstProbability = 0.15;
    rng(cfg.seed,'twister');

    events = readtable(fixedEventsCsv,'TextType','string');
    required = ["record","task","trial","cueIndex","cueDir", ...
        "accepted","peakToPeak","maxJump"];
    if ~all(ismember(required,string(events.Properties.VariableNames)))
        error('Frozen event CSV lacks required columns.');
    end
    if height(events) ~= 400 || sum(events.accepted == 1) ~= 350
        error('Frozen event CSV does not contain 400 trials / 350 accepted.');
    end
    if any(~ismember(events.cueDir,[-1,1])) || ...
            any(~ismember(events.accepted,[0,1]))
        error('Frozen event CSV has unexpected labels or acceptance flags.');
    end

    files = ["VisualCogA_Task-1.mat","VisualCogA_Task-2.mat", ...
             "VisualCogB_Task-1.mat","VisualCogB_Task-2.mat"];
    splitRows = table();
    calibrationRows = table();
    recoveryRows = table();
    exampleRows = table();
    analysisOffsets = -round(cfg.preSec*cfg.fs):round(cfg.postSec*cfg.fs);
    analysisTime = analysisOffsets(:)/cfg.fs;
    extendedOffsets = -round(cfg.extPreSec*cfg.fs): ...
        round(cfg.extPostSec*cfg.fs);
    extendedTime = extendedOffsets/cfg.fs;
    analysisBaseline = analysisTime < 0;
    extendedBaseline = extendedTime >= -cfg.preSec & extendedTime < 0;
    [bl,al] = butter(4,cfg.lowPassHz/(cfg.fs/2),'low');

    for fileIdx = 1:numel(files)
        fileName = files(fileIdx);
        sourceFile = fullfile(sourceDir,fileName);
        if ~isfile(sourceFile)
            error('Missing source file: %s',sourceFile);
        end
        source = load(sourceFile,'data','SampleRate');
        if ~isfield(source,'data') || ~isfield(source,'SampleRate') || ...
                source.SampleRate ~= cfg.fs || size(source.data,1) ~= 10
            error('Unexpected MAT structure or sampling rate: %s',fileName);
        end
        name = erase(fileName,".mat");
        task = str2double(extractAfter(name,"Task-"));
        rows = events(events.record == name,:);
        if height(rows) ~= 100 || any(rows.trial(:) ~= (1:100)')
            error('Frozen event rows are missing or unordered for %s.',name);
        end
        cue = source.data(8,:);
        cueStart = find(cue ~= 0 & [true,cue(1:end-1) == 0]);
        if numel(cueStart) ~= 100 || ...
                any(rows.cueIndex(:) ~= cueStart(:)) || ...
                any(rows.cueDir(:) ~= reshape( ...
                cue(cueStart(:)),[],1))
            error('Frozen event CSV does not match cue markers in %s.',fileName);
        end
        raw = double(source.data(1:3,:));
        nSamples = size(raw,2);
        if any(rows.cueIndex(rows.accepted == 1) + ...
                extendedOffsets(1) < 1) || ...
                any(rows.cueIndex(rows.accepted == 1) + ...
                extendedOffsets(end) > nSamples)
            error('An accepted trial lacks the extended reference window.');
        end

        [trainMask,validationMask,referenceMask,qualityScore] = ...
            fixed_partition(rows,cfg);
        if any(trainMask & validationMask) || ...
                any(referenceMask & validationMask) || ...
                any((trainMask | validationMask) & ~logical(rows.accepted))
            error('Simulation split violates train/validation separation.');
        end
        guardMask = logical(rows.accepted) & ...
            ~trainMask & ~validationMask;
        split = table(repmat(name,100,1),repmat(task,100,1), ...
            rows.trial,rows.cueIndex,rows.cueDir,rows.accepted, ...
            trainMask,guardMask,validationMask,referenceMask, ...
            qualityScore, ...
            'VariableNames',{'record','task','trial','cueIndex','cueDir', ...
            'accepted','train','guard','validation','reference', ...
            'qualityScore'});
        splitRows = [splitRows;split]; %#ok<AGROW>

        calibration = calibrate_noise(raw,rows,trainMask,cfg);
        for channelIdx = 1:3
            label = ["Fz","F3","F4"];
            c = calibration(channelIdx);
            calibrationRows = [calibrationRows;table(name,task, ...
                label(channelIdx),sum(trainMask), ...
                c.baselineLevelMad,c.baselineStdMedian, ...
                c.peakToPeakMedian,c.peakToPeakQ75,c.maxJumpQ75, ...
                c.driftAmplitude,c.shiftSd,c.coloredNoiseSd, ...
                c.spikeAmplitude,c.burstAmplitude, ...
                c.driftFrequencyLow,c.driftFrequencyHigh, ...
                c.coloredPhi,c.burstFrequency,c.burstDurationSec, ...
                cfg.spikeProbability,cfg.burstProbability, ...
                'VariableNames',{'record','task','channel', ...
                'nCalibration','baselineLevelMad','baselineStdMedian', ...
                'peakToPeakMedian','peakToPeakQ75','maxJumpQ75', ...
                'driftAmplitude','shiftSd','coloredNoiseSd', ...
                'spikeAmplitude','burstAmplitude', ...
                'driftFrequencyLow','driftFrequencyHigh', ...
                'coloredPhi','burstFrequency','burstDurationSec', ...
                'spikeProbability','burstProbability'})]; %#ok<AGROW>
        end

        for templateIdx = 1:numel(cfg.templateHp)
            templateHp = cfg.templateHp(templateIdx);
            [bhTemplate,ahTemplate] = butter(2, ...
                templateHp/(cfg.fs/2),'high');
            trainingFiltered = filtfilt(bl,al, ...
                filtfilt(bhTemplate,ahTemplate,raw'))';
            template = make_template(trainingFiltered,rows, ...
                referenceMask,extendedOffsets,extendedTime, ...
                extendedBaseline,cfg);
            clear trainingFiltered

            clean = zeros(3,nSamples);
            validationIds = find(validationMask);
            for j = 1:numel(validationIds)
                i = validationIds(j);
                sideIdx = 1 + (rows.cueDir(i) == 1);
                indices = rows.cueIndex(i) + extendedOffsets;
                clean(:,indices) = clean(:,indices) + ...
                    template(:,:,sideIdx);
            end
            oracleSignal = filtfilt(bl,al,clean')';
            oracle = estimate_heldout(oracleSignal,rows, ...
                validationMask,analysisOffsets,analysisBaseline,cfg);

            for replicate = 1:nReplicates
                % Reuse the exact contamination draw for both template
                % families, independently of loop order and cutoff.
                drawSeed = cfg.seed + 1000*fileIdx + replicate;
                rng(drawSeed,'twister');
                noise = synthesize_noise(nSamples,rows,validationMask, ...
                    calibration,cfg);
                synthetic = clean + noise;
                noisyLowPass = filtfilt(bl,al,synthetic')';
                noisyCurves = estimate_heldout(noisyLowPass,rows, ...
                    validationMask,analysisOffsets,analysisBaseline,cfg);
                for hpIdx = 1:numel(cfg.candidateHp)
                    highpass = cfg.candidateHp(hpIdx);
                    [bh,ah] = butter(2,highpass/(cfg.fs/2),'high');
                    filtered = filtfilt(bl,al, ...
                        filtfilt(bh,ah,synthetic'))';
                    recovered = estimate_heldout(filtered,rows, ...
                        validationMask,analysisOffsets,analysisBaseline,cfg);
                    for channelIdx = 1:3
                        label = ["Fz","F3","F4"];
                        nLeftReference = sum(referenceMask & ...
                            rows.cueDir == -1);
                        nRightReference = sum(referenceMask & ...
                            rows.cueDir == 1);
                        nLeftValidation = sum(validationMask & ...
                            rows.cueDir == -1);
                        nRightValidation = sum(validationMask & ...
                            rows.cueDir == 1);
                        for conditionIdx = 1:3
                            if conditionIdx < 3
                                trueWave = oracle(:,channelIdx, ...
                                    conditionIdx);
                                fitWave = recovered(:,channelIdx, ...
                                    conditionIdx);
                                noisyWave = noisyCurves(:,channelIdx, ...
                                    conditionIdx);
                                if conditionIdx == 1
                                    condition = "Left";
                                    nReference = nLeftReference;
                                    nValidation = nLeftValidation;
                                else
                                    condition = "Right";
                                    nReference = nRightReference;
                                    nValidation = nRightValidation;
                                end
                            else
                                trueWave = oracle(:,channelIdx,2) - ...
                                    oracle(:,channelIdx,1);
                                fitWave = recovered(:,channelIdx,2) - ...
                                    recovered(:,channelIdx,1);
                                noisyWave = noisyCurves(:,channelIdx,2) - ...
                                    noisyCurves(:,channelIdx,1);
                                condition = "Delta";
                                nReference = nLeftReference + ...
                                    nRightReference;
                                nValidation = nLeftValidation + ...
                                    nRightValidation;
                            end
                            m = recovery_metrics(trueWave,fitWave, ...
                                analysisTime);
                            deltaMeanError250_500 = NaN;
                            if conditionIdx == 3
                                deltaMeanError250_500 = ...
                                    m.meanError250_500;
                            end
                            c = calibration(channelIdx);
                            recoveryRows = [recoveryRows;table( ...
                                name,task,label(channelIdx),templateHp, ...
                                replicate,highpass,condition,nReference, ...
                                nValidation,drawSeed,m.waveformRmse, ...
                                m.waveformR,m.referenceMean250_500, ...
                                m.recoveredMean250_500, ...
                                m.meanError250_500, ...
                                m.referencePeak250_500, ...
                                m.recoveredPeak250_500,m.peakError, ...
                                m.referenceLatencySec,m.recoveredLatencySec, ...
                                m.latencyErrorSec,m.preStimRmse, ...
                                m.postStimRmse,deltaMeanError250_500, ...
                                c.driftAmplitude,c.shiftSd, ...
                                c.coloredNoiseSd,c.spikeAmplitude, ...
                                c.burstAmplitude, ...
                                'VariableNames',{'record','task','channel', ...
                                'templateHp','replicate','highpass', ...
                                'condition','nReference','nValidation', ...
                                'seed','waveformRmse','waveformR', ...
                                'referenceMean250_500', ...
                                'recoveredMean250_500', ...
                                'meanError250_500', ...
                                'referencePeak250_500', ...
                                'recoveredPeak250_500','peakError', ...
                                'referenceLatencySec', ...
                                'recoveredLatencySec','latencyErrorSec', ...
                                'preStimRmse','postStimRmse', ...
                                'deltaMeanError250_500','driftAmplitude', ...
                                'shiftSd','coloredNoiseSd', ...
                                'spikeAmplitude','burstAmplitude'})]; %#ok<AGROW>

                            if name == "VisualCogA_Task-2" && ...
                                    channelIdx == 2 && replicate == 1
                                nTime = numel(analysisTime);
                                exampleRows = [exampleRows;table( ...
                                    repmat(name,nTime,1), ...
                                    repmat(task,nTime,1), ...
                                    repmat(label(channelIdx),nTime,1), ...
                                    repmat(templateHp,nTime,1), ...
                                    repmat(replicate,nTime,1), ...
                                    repmat(highpass,nTime,1), ...
                                    repmat(condition,nTime,1), ...
                                    analysisTime,trueWave,noisyWave,fitWave, ...
                                    'VariableNames',{'record','task', ...
                                    'channel','templateHp','replicate', ...
                                    'highpass','condition','timeSec', ...
                                    'reference','noisyBeforeHighpass', ...
                                    'recovered'})]; %#ok<AGROW>
                            end
                        end
                    end
                end
            end
        end
        fprintf('%s: %d fixed accepted, %d train, %d held out.\n', ...
            name,sum(rows.accepted),sum(trainMask),sum(validationMask));
    end

    if sum(splitRows.accepted) ~= 350 || ...
            sum(splitRows.train) + sum(splitRows.validation) + ...
            sum(splitRows.guard) ~= 350
        error('Simulation trial split does not reconcile to 350, including guard trials.');
    end
    writetable(splitRows,fullfile(outputDir, ...
        'q1_filter_simulation_split.csv'));
    writetable(calibrationRows,fullfile(outputDir, ...
        'q1_filter_simulation_calibration.csv'));
    writetable(recoveryRows,fullfile(outputDir, ...
        'q1_filter_simulation_recovery.csv'));
    writetable(exampleRows,fullfile(outputDir, ...
        'q1_filter_simulation_example_waveform.csv'));
    fprintf('Semi-synthetic recovery outputs saved to %s\n',outputDir);
end

function [trainMask,validationMask,referenceMask,score] = ...
        fixed_partition(rows,cfg)
    n = height(rows);
    trainMask = false(n,1);
    validationMask = false(n,1);
    referenceMask = false(n,1);
    score = nan(n,1);
    for direction = [-1,1]
        trainIds = find(rows.accepted == 1 & ...
            rows.cueDir == direction & ...
            rows.trial <= cfg.lastTrainingTrial);
        validationIds = find(rows.accepted == 1 & ...
            rows.cueDir == direction & ...
            rows.trial >= cfg.firstValidationTrial);
        if numel(trainIds) < 15 || numel(validationIds) < 10
            error('Too few accepted training or validation trials per condition.');
        end
        trainMask(trainIds) = true;
        validationMask(validationIds) = true;
        peakScale = max(median(rows.peakToPeak(trainIds)),eps);
        jumpScale = max(median(rows.maxJump(trainIds)),eps);
        score(trainIds) = rows.peakToPeak(trainIds)/peakScale + ...
            rows.maxJump(trainIds)/jumpScale;
        [~,rank] = sortrows([score(trainIds),rows.trial(trainIds)], ...
            [1,2]);
        nReference = max(10,round( ...
            cfg.highQualityFraction*numel(trainIds)));
        referenceMask(trainIds(rank(1:nReference))) = true;
    end
end

function calibration = calibrate_noise(raw,rows,trainMask,cfg)
    ids = find(trainMask);
    preOffsets = -round(cfg.preSec*cfg.fs):-1;
    epochOffsets = -round(cfg.preSec*cfg.fs): ...
        round(cfg.postSec*cfg.fs);
    calibration = repmat(struct('baselineLevelMad',0, ...
        'baselineStdMedian',0,'peakToPeakMedian',0, ...
        'peakToPeakQ75',0,'maxJumpQ75',0, ...
        'driftAmplitude',0,'shiftSd',0,'coloredNoiseSd',0, ...
        'spikeAmplitude',0,'burstAmplitude',0, ...
        'driftFrequencyLow',0,'driftFrequencyHigh',0, ...
        'coloredPhi',0,'burstFrequency',0, ...
        'burstDurationSec',0),3,1);
    for e = 1:3
        levels = zeros(numel(ids),1);
        spread = zeros(numel(ids),1);
        peakToPeak = zeros(numel(ids),1);
        maxJump = zeros(numel(ids),1);
        lagNumerator = 0;
        lagDenominator = 0;
        burstFrequencies = zeros(numel(ids),1);
        excursionLengths = [];
        for j = 1:numel(ids)
            x = raw(e,rows.cueIndex(ids(j)) + preOffsets);
            levels(j) = median(x);
            spread(j) = std(x);
            centered = x-mean(x);
            lagNumerator = lagNumerator + sum( ...
                centered(1:end-1).*centered(2:end));
            lagDenominator = lagDenominator + sum( ...
                centered(1:end-1).^2);
            spectrum = abs(fft(centered.*hamming(numel(x))',256)).^2;
            f = (0:128)*cfg.fs/256;
            band = f >= 4 & f <= 20;
            candidates = f(band);
            bandPower = spectrum(1:129);
            [~,ix] = max(bandPower(band));
            burstFrequencies(j) = candidates(ix);
            robustSpread = 1.4826*median(abs(x-median(x)));
            excursions = abs(centered) > 1.5*robustSpread;
            boundaries = diff([false,excursions,false]);
            starts = find(boundaries == 1);
            ends = find(boundaries == -1);
            excursionLengths = [excursionLengths,ends-starts]; %#ok<AGROW>
            segment = raw(e,rows.cueIndex(ids(j)) + epochOffsets);
            peakToPeak(j) = max(segment)-min(segment);
            maxJump(j) = max(abs(diff(segment)));
        end
        levelMad = 1.4826*median(abs(levels-median(levels)));
        baselineStd = median(spread);
        p2pMedian = median(peakToPeak);
        p2pQ75 = prctile(peakToPeak,75);
        jumpQ75 = prctile(maxJump,75);
        c = calibration(e);
        c.baselineLevelMad = levelMad;
        c.baselineStdMedian = baselineStd;
        c.peakToPeakMedian = p2pMedian;
        c.peakToPeakQ75 = p2pQ75;
        c.maxJumpQ75 = jumpQ75;
        c.driftAmplitude = min(0.5*levelMad,p2pQ75/4);
        c.shiftSd = min(levelMad,p2pQ75/4);
        c.coloredNoiseSd = min(0.5*baselineStd,p2pQ75/20);
        c.spikeAmplitude = min(jumpQ75,p2pQ75/5);
        c.burstAmplitude = min(1.5*baselineStd,p2pQ75/6);
        trainingEnd = rows.cueIndex(rows.trial == ...
            cfg.lastTrainingTrial) + round(cfg.postSec*cfg.fs);
        continuousTrain = detrend(raw(e,1:trainingEnd));
        welchLength = min(32768,2^floor(log2( ...
            numel(continuousTrain)/2)));
        [power,f] = pwelch(continuousTrain, ...
            hamming(welchLength),floor(welchLength/2), ...
            welchLength,cfg.fs);
        c.driftFrequencyLow = dominant_frequency(power,f,0.02,0.08);
        c.driftFrequencyHigh = dominant_frequency(power,f,0.08,0.20);
        c.coloredPhi = min(0.995,max(0.4, ...
            lagNumerator/max(lagDenominator,eps)));
        c.burstFrequency = median(burstFrequencies);
        if isempty(excursionLengths)
            c.burstDurationSec = 0.08;
        else
            c.burstDurationSec = min(0.25,max(0.04, ...
                3*median(excursionLengths)/cfg.fs));
        end
        calibration(e) = c;
    end
end

function template = make_template(filtered,rows,referenceMask,offsets, ...
        time,baselineMask,cfg)
    nTime = numel(offsets);
    template = zeros(3,nTime,2);
    taper = ones(1,nTime);
    early = time < -0.7;
    taper(early) = max(0,min(1,(time(early) + ...
        cfg.extPreSec)/(cfg.extPreSec-0.7)));
    late = time > 1.0;
    taper(late) = max(0,min(1,(cfg.extPostSec - ...
        time(late))/(cfg.extPostSec-1.0)));
    for sideIdx = 1:2
        direction = 2*sideIdx-3;
        ids = find(referenceMask & rows.cueDir == direction);
        for e = 1:3
            x = zeros(numel(ids),nTime);
            for j = 1:numel(ids)
                x(j,:) = filtered(e,rows.cueIndex(ids(j)) + offsets);
                x(j,:) = x(j,:) - mean(x(j,baselineMask));
            end
            fitted = smooth_curve(huber_location(x,cfg), ...
                cfg.smoothLambda)';
            fitted = fitted - mean(fitted(baselineMask));
            template(e,:,sideIdx) = fitted.*taper;
        end
    end
end

function noise = synthesize_noise(nSamples,rows,validationMask,cal,cfg)
    noise = zeros(3,nSamples);
    t = (0:nSamples-1)/cfg.fs;
    ids = find(validationMask);
    for e = 1:3
        c = cal(e);
        f1 = c.driftFrequencyLow;
        f2 = c.driftFrequencyHigh;
        slow = c.driftAmplitude*(sin(2*pi*f1*t + 2*pi*rand) + ...
            0.5*sin(2*pi*f2*t + 2*pi*rand))/1.5;
        phi = c.coloredPhi;
        colored = filter(1,[1,-phi],randn(1,nSamples));
        colored = colored-mean(colored);
        colored = c.coloredNoiseSd*colored/max(std(colored),eps);
        noise(e,:) = slow + colored;
        for j = 1:numel(ids)
            cue = rows.cueIndex(ids(j));
            % A smooth trial-wise offset is calibrated from training-trial
            % baseline levels; its plateau spans the ERP interval.
            indices = cue-round(1.0*cfg.fs): ...
                cue+round(1.4*cfg.fs);
            shape = ones(1,numel(indices));
            ramp = round(0.25*cfg.fs);
            shape(1:ramp) = 0.5-0.5*cos( ...
                pi*(0:ramp-1)/(ramp-1));
            shape(end-ramp+1:end) = fliplr(shape(1:ramp));
            noise(e,indices) = noise(e,indices) + ...
                c.shiftSd*randn*shape;
            if rand < cfg.spikeProbability
                spikeIndex = cue + randi([ ...
                    -round(cfg.preSec*cfg.fs), ...
                    round(cfg.postSec*cfg.fs)]);
                noise(e,spikeIndex) = noise(e,spikeIndex) + ...
                    c.spikeAmplitude*(2*(rand > 0.5)-1);
            end
            if rand < cfg.burstProbability
                center = cue + randi([0,round(0.7*cfg.fs)]);
                halfWidth = round(c.burstDurationSec*cfg.fs/2);
                burstIndices = center-halfWidth:center+halfWidth;
                envelope = exp(-0.5*((burstIndices-center)/ ...
                    max(halfWidth/2,1)).^2);
                freq = c.burstFrequency*(0.9+0.2*rand);
                oscillation = sin(2*pi*freq* ...
                    (burstIndices-center)/cfg.fs + 2*pi*rand);
                noise(e,burstIndices) = noise(e,burstIndices) + ...
                    c.burstAmplitude*envelope.*oscillation;
            end
        end
    end
end

function selected = dominant_frequency(power,f,low,high)
    band = f >= low & f < high;
    if ~any(band)
        error('Insufficient spectral resolution for drift calibration.');
    end
    candidates = f(band);
    bandPower = power(band);
    [~,ix] = max(bandPower);
    selected = candidates(ix);
end

function curves = estimate_heldout(signal,rows,validationMask, ...
        offsets,baselineMask,cfg)
    curves = zeros(numel(offsets),3,2);
    for sideIdx = 1:2
        direction = 2*sideIdx-3;
        ids = find(validationMask & rows.cueDir == direction);
        for e = 1:3
            x = zeros(numel(ids),numel(offsets));
            for j = 1:numel(ids)
                x(j,:) = signal(e,rows.cueIndex(ids(j)) + offsets);
                x(j,:) = x(j,:) - mean(x(j,baselineMask));
            end
            curves(:,e,sideIdx) = smooth_curve( ...
                huber_location(x,cfg),cfg.smoothLambda);
        end
    end
end

function m = recovery_metrics(reference,recovered,time)
    reference = reference(:);
    recovered = recovered(:);
    difference = recovered-reference;
    window = time >= 0.25 & time < 0.5;
    baseline = time < 0;
    post = time >= 0;
    m.waveformRmse = sqrt(mean(difference.^2));
    cc = corrcoef(reference,recovered);
    if numel(cc) == 4
        m.waveformR = cc(1,2);
    else
        m.waveformR = NaN;
    end
    m.referenceMean250_500 = mean(reference(window));
    m.recoveredMean250_500 = mean(recovered(window));
    m.meanError250_500 = m.recoveredMean250_500 - ...
        m.referenceMean250_500;
    windowTimes = time(window);
    [m.referencePeak250_500,referencePeakIndex] = ...
        max(reference(window));
    [m.recoveredPeak250_500,recoveredPeakIndex] = ...
        max(recovered(window));
    m.peakError = m.recoveredPeak250_500 - ...
        m.referencePeak250_500;
    m.referenceLatencySec = windowTimes(referencePeakIndex);
    m.recoveredLatencySec = windowTimes(recoveredPeakIndex);
    m.latencyErrorSec = m.recoveredLatencySec - ...
        m.referenceLatencySec;
    m.preStimRmse = sqrt(mean(difference(baseline).^2));
    m.postStimRmse = sqrt(mean(difference(post).^2));
end

function location = huber_location(x,cfg)
    % Identical constant, scale estimate, stopping criterion and iterations
    % to the present run_q1.m estimator; no candidate estimator is involved.
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

function smooth = smooth_curve(values,lambda)
    n = numel(values);
    secondDifference = spdiags([ones(n-2,1),-2*ones(n-2,1), ...
        ones(n-2,1)],[0,1,2],n-2,n);
    smooth = (speye(n)+lambda*(secondDifference'*secondDifference)) ...
        \values(:);
end
