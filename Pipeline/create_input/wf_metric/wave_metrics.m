function meas_out = wave_metrics(mw, input, dataPath)

% get sample rate from meta file if possible, else use default created in
% main
metaFolder = fileparts(dataPath);
metaFiles = dir(fullfile(metaFolder, '*.meta'));
if isempty(metaFiles)
    % use default sample rate
    fs = input.fs;
else
    metaFilePath = fullfile(metaFolder, metaFiles(1).name);
    fs = get_sample_rate(metaFilePath, input.fs);
end

chan_pos = input.chan_pos;
[nUnit, nChan, nt] = size(mw);

pp_all = squeeze(max(mw,[],3)-min(mw,[],3));
[pp_unit, pk_chan] = max(pp_all(:,1:length(input.chan_map)),[],2);
% get background of pp -- this is always nonzero, since it includes the
% noise.
ppvalues = reshape(pp_all,[nUnit*nChan,1]);
edges = 0:2:1000;
ppDist = histcounts(ppvalues,edges);
[~,maxbin] = max(ppDist);
backVal = (edges(maxbin+1)+edges(maxbin))/2;
pp_all = pp_all - 2*backVal;
neg_val = (pp_all<0);
pp_all(neg_val) = 0;

pk_wave = zeros(nUnit,nt);
for i = 1:nUnit
    pk_wave(i,:) = squeeze(mw(i,pk_chan(i),:));
end
min_unit = min(pk_wave,[],2);

meas_out = zeros(nUnit,2);
meas_out(:,1) = pp_unit;
meas_out(:,2) = min_unit;

pk_wave_upsamp = resample(pk_wave,200,input.ts,'Dimension',2);
time_conv = (200/input.ts)*fs/1000; % to convert from points to ms
[~,nt_up] = size(pk_wave_upsamp);
timestamps = (1:nt_up)/time_conv;
window = 20;  % for slope measurement

for i = 1:nUnit
    cw = squeeze(pk_wave_upsamp(i,:));
    npts_up = numel(cw);
    [~,trough_idx] = min(cw);
    [~,peak_idx] = max(cw);
    bPos = (cw(peak_idx) > cw(trough_idx));
    back_level = mean(cw(1:5));
    
    % duration/peak to trough time
    if bPos
        % this is a postive going peak, search from peak index to the next
        % trough
        [~, loc_trough_idx] = min(cw(peak_idx:end));
        duration = timestamps(peak_idx + loc_trough_idx-1) - timestamps(peak_idx);
    else
        % search forward from trough to find peak
        [~,loc_peak_idx] = max(cw(trough_idx:end));
        duration = timestamps(trough_idx + loc_peak_idx-1) - timestamps(trough_idx);
    end
    
    % full width at half maximum (fwhm)
    fwhm = 0;
    if bPos
        threshold = cw(peak_idx) * 0.5;
        thresh_crossing_1 = min(timestamps(cw(1:peak_idx)>threshold));
        tc2_ind = peak_idx + min(find(cw(peak_idx:end)<threshold))-1;
        thresh_crossing_2 = timestamps(tc2_ind);
    else
        threshold = cw(trough_idx) * 0.5;
        thresh_crossing_1 = min(timestamps(cw(1:trough_idx)<threshold));
        tc2_ind = trough_idx + min(find(cw(trough_idx:end)>threshold))-1;
        thresh_crossing_2 = timestamps(tc2_ind);
    end
    if ~isempty(thresh_crossing_2) && ~isempty(thresh_crossing_1)
        fwhm = thresh_crossing_2 - thresh_crossing_1;
    end
    
    % PT ratio
    PT_ratio = 0;
    if cw(trough_idx) ~= 0
        PT_ratio = abs(cw(peak_idx)/cw(trough_idx));
    end
    
    % height of pre-peak, for negative going spikes
    pre_peak = 0;
    if ~bPos
        pre_peak = max(cw(1:trough_idx)) - back_level;
    end
    
    % recovery slope
    if bPos
        if trough_idx > npts_up-20
            window = npts_up - trough_idx;
        end
        
        x = timestamps(trough_idx:trough_idx+window-1);
        X = [x', ones(window,1)];
        
        
        if (peak_idx+window-1) > length(cw)
            
            X = X(1:length(cw(peak_idx:end)), :);
            lreg = regress(cw(peak_idx:length(cw))',X);
        else
            lreg = regress(cw(peak_idx:peak_idx+window-1)',X);
        end
        
        
        lreg(1) = -1*lreg(1);
    else
        % fit recover after repolarization (peak down to baseline)
        if peak_idx > npts_up-20
            window = npts_up - peak_idx;
        end
        x = timestamps(peak_idx:peak_idx+window-1);
        X = [x', ones(window,1)];
        
        
        
        if (peak_idx+window-1) > length(cw)
            X = X(1:length(peak_idx:end));
            lreg(regress(cw(peak_idx:end)',X));
        else
            lreg = regress(cw(peak_idx:peak_idx+window-1)',X);
        end
    end
    
    recovery_slope = lreg(1) * 1e-3; % convert to V/s
    
    
    
    
    % 2D waveform metrics
    pp_unit = squeeze(pp_all(i,:))';
    
    % Julien Boussard - style fit of peak-to-peak voltage vs position
    % if background sub pp_unit > 60 uV, attempt a fit of the
    % background subtracted pp_all
    if max(squeeze(pp_all(i,:))) > 60
        fitvals = fit_loc(i, pp_all, chan_pos);
        fitX = fitvals(1);
        fitZ = fitvals(2);
        fitY = fitvals(3);
    else
        fitX = chan_pos(pk_chan(i),1);
        fitZ = chan_pos(pk_chan(i),2);
        fitY = -1;      % a marker for no fit
    end
    
    % calculate spread of waveforms in z
    sp_thresh = 0.2*max(pp_unit);
    chan_above = pp_unit > sp_thresh;
    chan_above = chan_above(1:length(chan_pos));
    zmax = max(chan_pos(chan_above,2));
    zmin = min(chan_pos(chan_above,2));
    z_sp = 0;
    if ~isempty(zmax) && ~isempty(zmin)
        z_sp = zmax(1)-zmin(1);
    end
    
    
    % fill in meas_out for this unit
    meas_out(i,3) = duration;
    meas_out(i,4) = fwhm;
    meas_out(i,5) = PT_ratio;
    meas_out(i,6) = pre_peak;
    meas_out(i,7) = recovery_slope;
    meas_out(i,8) = z_sp;
    meas_out(i,9) = fitX;
    meas_out(i,10) = fitZ;
    meas_out(i,11) = fitY;
    
end

end

function fs = get_sample_rate(metaFile, default)
fs = default; % set default value and override if found
fid = fopen(metaFile, 'r');
if fid == -1
    warning(['Could not open file: ', metaFile]);
    return
end
% read file
while ~feof(fid)
    line = fgetl(fid);
    if contains(line, 'imSampRate')
        fclose(fid);
        parts = strsplit(line, '=');
        if numel(parts) == 2
            fs = str2double(parts{2});
        else
            warning('imSampRate line has more than 1 value')
        end
        return
    end
end
fclose(fid);
end
