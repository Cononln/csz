function q1_future_event_leakage(sourceDir, outputDir, fixedEventsCsv)
%Q1_FUTURE_EVENT_LEAKAGE Bounded cue-window sensitivity to later events.
%   Replaces copied continuous EEG in target/click onset neighborhoods with
%   linear endpoint interpolation. Original MAT files and accepted flags are
%   untouched. This is a counterfactual perturbation, not artifact removal
%   or identification of a biological source.

    arguments
        sourceDir (1,1) string
        outputDir (1,1) string
        fixedEventsCsv (1,1) string
    end
    if ~isfolder(sourceDir) || ~isfile(fixedEventsCsv)
        error('Source directory or fixed event table does not exist.');
    end
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    events = readtable(fixedEventsCsv);
    required = ["record","task","trial","cueIndex","cueDir", ...
        "targetIndex","clickIndex","accepted"];
    if ~all(ismember(required,string(events.Properties.VariableNames))) || ...
            height(events) ~= 400 || sum(logical(events.accepted)) ~= 350
        error('Expected canonical 400 events with exactly 350 accepted.');
    end

    files = ["VisualCogA_Task-1.mat","VisualCogA_Task-2.mat", ...
        "VisualCogB_Task-1.mat","VisualCogB_Task-2.mat"];
    rows = table();
    for fileNo = 1:numel(files)
        file = files(fileNo);
        token = regexp(char(file),'VisualCog([AB])_Task-([12])\.mat$', ...
            'tokens','once');
        record = string(token{1});
        task = str2double(token{2});
        name = "VisualCog" + record + "_Task-" + string(task);
        ev = sortrows(events(string(events.record)==name,:), 'trial');
        if height(ev)~=100 || ~isequal(ev.trial(:),(1:100)')
            error('Canonical event order does not match %s.',file);
        end
        source = load(fullfile(sourceDir,file),'data','SampleRate');
        fs = double(source.SampleRate);
        if fs ~= 256 || size(source.data,1) < 3
            error('Unexpected signal layout in %s.',file);
        end
        raw = double(source.data(1:3,:));
        cue = double(ev.cueIndex);
        accepted = logical(ev.accepted);
        direction = double(ev.cueDir);
        if any(~ismember(direction,[-1,1])) || ...
                any(double(source.data(8,cue)) ~= direction') || ...
                any(cue(accepted)-round(0.2*fs)<1 | ...
                    cue(accepted)+round(0.8*fs)>size(raw,2))
            error('Invalid cue, direction or epoch bounds in %s.',file);
        end
        variant = ["original","both_015","both_030", ...
            "target_030","click_030"];
        width = [0,0.15,0.30,0.30,0.30];
        useTarget = [false,true,true,true,false];
        useClick = [false,true,true,false,true];
        if ~any(isfinite(ev.clickIndex))
            variant(end) = [];
            width(end) = [];
            useTarget(end) = [];
            useClick(end) = [];
        end
        [bl,al] = butter(4,30/(fs/2),'low');
        offsets = -round(0.2*fs):round(0.8*fs);
        time = offsets/fs;
        baseline = time < 0;
        mid = time >= 0.25 & time < 0.5;
        nAccepted = sum(accepted);
        smoothMatrix = local_smoothing_matrix(numel(time),12);
        channels = ["Fz","F3","F4"];
        for hp = [0.1,0.5]
            [bh,ah] = butter(2,hp/(fs/2),'high');
            estimates = nan(numel(variant),3,2);
            replaced = zeros(numel(variant),1);
            for v = 1:numel(variant)
                [copy,replaced(v)] = interpolate_events(raw,ev,fs, ...
                    width(v),useTarget(v),useClick(v));
                filtered = filtfilt(bl,al,filtfilt(bh,ah,copy'))';
                epoch = nan(nAccepted,numel(time),3);
                acceptedCue = cue(accepted);
                acceptedDir = direction(accepted);
                for k = 1:nAccepted
                    seg = filtered(:,acceptedCue(k)+offsets);
                    seg = seg-mean(seg(:,baseline),2);
                    epoch(k,:,:) = reshape(seg',1,numel(time),3);
                end
                for e = 1:3
                    left = reshape(epoch(acceptedDir==-1,:,e), ...
                        sum(acceptedDir==-1),numel(time));
                    right = reshape(epoch(acceptedDir==1,:,e), ...
                        sum(acceptedDir==1),numel(time));
                    leftHuber = local_huber(left);
                    rightHuber = local_huber(right);
                    deltaHuber = rightHuber-leftHuber;
                    leftFit = smoothMatrix\leftHuber(:);
                    rightFit = smoothMatrix\rightHuber(:);
                    estimates(v,e,1) = mean(deltaHuber(mid));
                    estimates(v,e,2) = mean(rightFit(mid)-leftFit(mid));
                end
            end
            for v = 1:numel(variant)
                for e = 1:3
                    deltaHuber = estimates(v,e,1);
                    deltaFit = estimates(v,e,2);
                    row = table(record,task,channels(e),hp,variant(v), ...
                        width(v),useTarget(v),useClick(v), ...
                        nAccepted,replaced(v),replaced(v)/size(raw,2), ...
                        deltaHuber,deltaFit, ...
                        deltaHuber-estimates(1,e,1), ...
                        deltaFit-estimates(1,e,2), ...
                        'VariableNames',{'record','task','channel', ...
                        'highpass','variant','half_width_sec', ...
                        'target_replaced','click_replaced','n_accepted', ...
                        'modified_samples','modified_fraction', ...
                        'delta_huber_250_500','delta_fit_250_500', ...
                        'change_huber','change_fit'});
                    rows = [rows;row]; %#ok<AGROW>
                end
            end
        end
    end
    outFile = fullfile(outputDir,'q1_future_event_leakage.csv');
    writetable(rows,outFile);
    fprintf('Counterfactual filter sensitivity saved to %s\n',outFile);
end

function [copy,count] = interpolate_events(raw,ev,fs,halfWidth, ...
        useTarget,useClick)
    copy = raw;
    covered = false(1,size(raw,2));
    if halfWidth == 0
        count = 0;
        return
    end
    fields = strings(0,1);
    if useTarget
        fields(end+1) = "targetIndex";
    end
    if useClick
        fields(end+1) = "clickIndex";
    end
    for fieldNo = 1:numel(fields)
        index = double(ev.(char(fields(fieldNo))));
        index = index(isfinite(index));
        for j = 1:numel(index)
            radius = round(halfWidth*fs);
            lo = max(2,index(j)-radius);
            hi = min(size(raw,2)-1,index(j)+radius);
            if any(covered(lo:hi))
                error('Counterfactual event windows overlap.');
            end
            covered(lo:hi) = true;
            alpha = (1:(hi-lo+1))/(hi-lo+2);
            copy(:,lo:hi) = copy(:,lo-1) + ...
                (copy(:,hi+1)-copy(:,lo-1))*alpha;
        end
    end
    count = sum(covered);
end

function result = local_huber(x)
    result = median(x,1);
    scale = 1.4826*median(abs(x-result),1);
    scale = max(scale,1e-9);
    for iteration = 1:15
        residual = abs(x-result);
        weight = min(1,1.345*scale./max(residual,1e-12));
        updated = sum(weight.*x,1)./sum(weight,1);
        if max(abs(updated-result),[],'all') < 1e-8
            result = updated;
            break
        end
        result = updated;
    end
end

function matrix = local_smoothing_matrix(n,lambda)
    second = spdiags([ones(n-2,1),-2*ones(n-2,1), ...
        ones(n-2,1)],[0,1,2],n-2,n);
    matrix = speye(n)+lambda*(second'*second);
end
