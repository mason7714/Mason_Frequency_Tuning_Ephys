%% For frequency Tuning Analysis:
%% Make clean EPA event file for freqTuning using ComputerTimestamp_6
% concat folder needs to be placed in the data folder for this to run

clear S D

savePath = '/mnt/CL_8TB_4/Mason/Ephys/DorsalAUD_PL/Sorting/SUBJ-ID-1180/260419_concat';
infoPath = '/mnt/CL_8TB_4/Mason/Ephys/DorsalAUD_PL/Sorting/SUBJ-ID-1180/260419_concat';
DataPath = '/mnt/CL_8TB_4/Mason/Ephys/DorsalAUD_PL/Sorting/SUBJ-ID-1180/260419_concat/kilosort4';
datFile = fix_dat_path(DataPath);

S = epa.load.phy(DataPath);

% freqSessionName = string(S(4).Name);
% 
% d = dir(fullfile(infoPath, '*trialInfo.csv'));
% d = d(~contains({d.name}, 'clean'));
% d = d(~contains({d.name}, 'concat'));
% d = d(contains({d.name}, freqSessionName));

freqSessionName = string(S(4).Name);

% Search recursively
d = dir(fullfile(infoPath, '**', '*trialInfo.csv'));

% Remove generated files
d = d(~contains({d.name}, 'clean'));
d = d(~contains({d.name}, 'concat'));

% More flexible matching
d = d(contains(string({d.name}), freqSessionName) | ...
      contains(string({d.name}), 'freqTuning'));

if isempty(d)
    disp("No freqTuning trialInfo found. Available trialInfo files:")
    disp(string({dir(fullfile(infoPath, '**', '*trialInfo.csv')).name})')
    error('No trialInfo file found for S(4)');
end

trialInfoFile = fullfile(d(1).folder, d(1).name);
disp("Using trialInfo file:")
disp(trialInfoFile)

% if isempty(d)
%     error('No trialInfo file found for S(4)');
% end

T = readtable(fullfile(d(1).folder, d(1).name));

%% Use Computer Time stamps
Trial_onset  = double(T.ComputerTimestamp_6);
Trial_offset = Trial_onset + double(T.StimDuration) ./ 1000;

%% Build EPA table
T_clean = table();
T_clean.Trial_onset  = Trial_onset;
T_clean.Trial_offset = Trial_offset;
T_clean.dBSPL        = double(T.dBSPL);
T_clean.Freq         = double(T.Freq);
T_clean.TrialID      = double(T.TrialID);

%% Save
cleanInfoPath = fullfile(infoPath, 'EPA_clean_events_freqTuning');
if ~exist(cleanInfoPath, 'dir')
    mkdir(cleanInfoPath);
end

cleanFile = fullfile(cleanInfoPath, freqSessionName + "_trialInfo.csv");
writetable(T_clean, cleanFile);

disp("Saved clean file:")
disp(cleanFile)

%% Load into EPA
epa.load.events(S, cleanInfoPath);

%% Verify
sess = S(4);

E = sess.find_Event('dBSPL');
disp(string(E.Name))
disp(size(E.OnOffTimes))
disp(size(E.Values))

E2 = sess.find_Event('Freq');
disp(string(E2.Name))
disp(size(E2.OnOffTimes))
disp(size(E2.Values))

%% Launch browser if you want
D = epa.DataBrowser;

%% PSTHs (raw data no baseline subtraction)

sess = S(4);

Clusternumber = length(sess.Clusters);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

% --- PSTH parameters ---
win = [-0.1 0.3];
bin_size = 0.01;
edges = win(1):bin_size:(win(2) + bin_size);
t = edges(1:end-1) + bin_size/2;

% Columns: frequency low to high
freq_levels = sort(unique(E_freq.Values), 'ascend');

% Rows: dB SPL high to low
db_levels = sort(unique(E_db.Values), 'descend');

n_freq = length(freq_levels);
n_db   = length(db_levels);

% --- SAVE PATH ---
save_dir = fullfile(savePath, 'PSTH');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    % --- Compute PSTHs ---
    all_psth = cell(n_db, n_freq);
    global_max = 0;

    for d = 1:n_db
        
        target_db = db_levels(d);
        
        for f = 1:n_freq
            
            target_freq = freq_levels(f);
            
            trials = E_freq.Values == target_freq & E_db.Values == target_db;
            onsets = E_freq.OnOffTimes(trials, 1);
            
            all_counts = [];
            
            for k = 1:length(onsets)
                rel_spikes = spikes - onsets(k);
                rel_spikes = rel_spikes(rel_spikes >= win(1) & rel_spikes <= win(2));
                all_counts(k,:) = histcounts(rel_spikes, edges);
            end
            
            if isempty(all_counts)
                psth = zeros(1, length(t));
            else
                psth = mean(all_counts, 1) / bin_size;
            end
            
            all_psth{d,f} = psth;
            global_max = max(global_max, max(psth));
        end
    end

    if global_max == 0
        global_max = 1;
    end

    % --- Plot ---
    figure(1)
    clf
    set(gcf, 'Position', [100 100 1800 1000])

    for d = 1:n_db
        
        for f = 1:n_freq
            
            subplot(n_db, n_freq, (d-1)*n_freq + f)
            
            psth = all_psth{d,f};
            
            % BLACK PSTH
            bar(t, psth, 1, 'FaceColor', 'k', 'EdgeColor', 'none')
            
            xlim(win)
            ylim([0 global_max * 1.2])
            
            if d == 1
                title(sprintf('%d Hz', freq_levels(f)))
            end
            
            if f == 1
                ylabel(sprintf('%d dB\nFR (Hz)', db_levels(d)))
            end
            
            if d == n_db
                xlabel('Time (s)')
            end
        end
    end

    fig_title = sprintf('%s | dB x Frequency PSTH', cluster.Name);
    sgtitle(fig_title)

    % --- CLEAN FILENAME ---
    safe_name = regexprep(fig_title, '[^\w]', '_');
    save_file = fullfile(save_dir, [safe_name '.pdf']);

    % --- SAVE ---
    exportgraphics(gcf, save_file, 'Resolution', 300)

    % --- CLOSE FIGURE ---
    close(gcf)

end
%% Verify if there should be any blank PSTHs for clusters that do not have spikes within the PSTH window
sess = S(4);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

win = [-0.1 0.3];

for i = 1:length(sess.Clusters)

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    total_spikes_session = length(spikes);
    total_spikes_psth = 0;

    onsets = E_freq.OnOffTimes(:,1);

    for k = 1:length(onsets)
        total_spikes_psth = total_spikes_psth + ...
            sum(spikes >= onsets(k)+win(1) & spikes < onsets(k)+win(2));
    end

    if total_spikes_psth == 0
        fprintf('Cluster index %d | ID %d | %s | session spikes = %d | PSTH spikes = %d\n', ...
            i, cluster.ID, cluster.Name, total_spikes_session, total_spikes_psth);
    end
end
%% Heat Maps uses baseline subtraction prior to trial
sess = S(4);

Clusternumber = length(sess.Clusters);

% create color map
function cmap = redblue(m)
    if nargin < 1, m = 256; end
    
    r = [(0:m/2-1)'/(m/2); ones(m/2,1)];
    g = [(0:m/2-1)'/(m/2); (m/2-1:-1:0)'/(m/2)];
    b = [ones(m/2,1); (m/2-1:-1:0)'/(m/2)];
    
    cmap = [r g b];
end

%find global min and max across all clusters
global_min = inf;
global_max = -inf;

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for d = 1:n_db
        for f = 1:n_freq
            
            trials = E_freq.Values == freq_levels(f) & E_db.Values == db_levels(d);
            onsets = E_freq.OnOffTimes(trials,1);

            stim_FRs = [];
            base_FRs = [];

            for k = 1:length(onsets)
                stim_spikes = sum(spikes >= onsets(k) & spikes < onsets(k)+0.2);
                base_spikes = sum(spikes >= onsets(k)-0.1 & spikes < onsets(k));

                stim_FRs(end+1) = stim_spikes / 0.2;
                base_FRs(end+1) = base_spikes / 0.1;
            end

            FR_map(d,f) = mean(stim_FRs - base_FRs, 'omitnan');
        end
    end

    global_min = min(global_min, min(FR_map(:)));
    global_max = max(global_max, max(FR_map(:)));
end

max_abs = max(abs([global_min global_max]));

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

% --- Response window relative to tone onset ---
response_win = [0 0.2];      % stimulus window, 0–200 ms
baseline_win = [-0.1 0];     % optional baseline window

freq_levels = sort(unique(E_freq.Values), 'ascend');   % columns
db_levels   = sort(unique(E_db.Values), 'descend');    % rows

n_freq = length(freq_levels);
n_db   = length(db_levels);

% --- SAVE PATH ---
save_dir = fullfile(savePath, 'TuningHeatmaps');
% save_dir = 'G:\My Drive\Data\Dorsal AUD Recordings\FrequenyTuningTests\Data\260331_concat\\TuningHeatmaps\';
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for d = 1:n_db
        
        target_db = db_levels(d);
        
        for f = 1:n_freq
            
            target_freq = freq_levels(f);
            
            trials = E_freq.Values == target_freq & E_db.Values == target_db;
            onsets = E_freq.OnOffTimes(trials, 1);
            
            stim_FRs = [];
            base_FRs = [];
            
            for k = 1:length(onsets)
                
                stim_start = onsets(k) + response_win(1);
                stim_end   = onsets(k) + response_win(2);
                
                base_start = onsets(k) + baseline_win(1);
                base_end   = onsets(k) + baseline_win(2);
                
                stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
                base_spikes = sum(spikes >= base_start & spikes < base_end);
                
                stim_FRs(end+1) = stim_spikes / diff(response_win);
                base_FRs(end+1) = base_spikes / diff(baseline_win);
            end
            
            % Baseline-subtracted firing rate
            FR_map(d,f) = mean(stim_FRs - base_FRs, 'omitnan');
            
            % If you do NOT want baseline subtraction, use this instead:
            % FR_map(d,f) = mean(stim_FRs, 'omitnan');
        end
    end

    % --- Plot heatmap ---
    figure(1)
clf

imagesc(freq_levels, 1:n_db, FR_map)
clim([-max_abs max_abs])

set(gca, 'XScale', 'log')
set(gca, 'YDir', 'reverse')

yticks(1:n_db)
yticklabels(string(db_levels))

xticks(freq_levels)
xtickangle(45)

colorbar
xlabel('Frequency (Hz)')
ylabel('dB SPL')

title(sprintf('%s | Frequency tuning heatmap', cluster.Name))

colormap(redblue)
% --- Save ---
    fig_title = sprintf('%s_Frequency_Tuning_Heatmap', cluster.Name);
    safe_name = regexprep(fig_title, '[^\w]', '_');
    save_file = fullfile(save_dir, [safe_name '.pdf']);
    
    exportgraphics(gcf, save_file, 'Resolution', 300)
    close(gcf)

end

%% Heatmaps: log1p stimulus FR - log1p baseline FR

sess = S(4);

Clusternumber = length(sess.Clusters);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

% --- Response window relative to tone onset ---
response_win = [0 0.2];      % 0–200 ms after tone onset
baseline_win = [-0.1 0];     % 100 ms before tone onset

freq_levels = sort(unique(E_freq.Values), 'ascend');   % columns
db_levels   = sort(unique(E_db.Values), 'descend');    % rows

n_freq = length(freq_levels);
n_db   = length(db_levels);

% --- SAVE PATH ---
save_dir = fullfile(savePath, 'TuningHeatmaps_log1p');
% save_dir = 'G:\My Drive\Data\Dorsal AUD Recordings\FrequenyTuningTests\Data\260331_concat\TuningHeatmaps_log1p\';
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

% Find global min/max across all clusters

global_min = inf;
global_max = -inf;

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for d = 1:n_db
        for f = 1:n_freq
            
            trials = E_freq.Values == freq_levels(f) & E_db.Values == db_levels(d);
            onsets = E_freq.OnOffTimes(trials, 1);

            trial_vals = [];

            for k = 1:length(onsets)

                stim_start = onsets(k) + response_win(1);
                stim_end   = onsets(k) + response_win(2);

                base_start = onsets(k) + baseline_win(1);
                base_end   = onsets(k) + baseline_win(2);

                stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
                base_spikes = sum(spikes >= base_start & spikes < base_end);

                stim_FR = stim_spikes / diff(response_win);
                base_FR = base_spikes / diff(baseline_win);

                trial_vals(end+1) = log1p(stim_FR) - log1p(base_FR);
            end

            FR_map(d,f) = mean(trial_vals, 'omitnan');
        end
    end

    global_min = min(global_min, min(FR_map(:), [], 'omitnan'));
    global_max = max(global_max, max(FR_map(:), [], 'omitnan'));
end

max_abs = max(abs([global_min global_max]));

% Plot heatmap for each cluster

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for d = 1:n_db
        
        target_db = db_levels(d);
        
        for f = 1:n_freq
            
            target_freq = freq_levels(f);
            
            trials = E_freq.Values == target_freq & E_db.Values == target_db;
            onsets = E_freq.OnOffTimes(trials, 1);
            
            trial_vals = [];
            
            for k = 1:length(onsets)
                
                stim_start = onsets(k) + response_win(1);
                stim_end   = onsets(k) + response_win(2);
                
                base_start = onsets(k) + baseline_win(1);
                base_end   = onsets(k) + baseline_win(2);
                
                stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
                base_spikes = sum(spikes >= base_start & spikes < base_end);
                
                stim_FR = stim_spikes / diff(response_win);
                base_FR = base_spikes / diff(baseline_win);
                
                trial_vals(end+1) = log1p(stim_FR) - log1p(base_FR);
            end
            
            FR_map(d,f) = mean(trial_vals, 'omitnan');
        end
    end

    figure(1)
    clf

    imagesc(freq_levels, 1:n_db, FR_map)
    clim([-max_abs max_abs])

    set(gca, 'XScale', 'log')
    set(gca, 'YDir', 'reverse')

    yticks(1:n_db)
    yticklabels(string(db_levels))

    xticks(freq_levels)
    xtickangle(45)

    cb = colorbar;
    ylabel(cb, 'log1p(stim FR) - log1p(baseline FR)')

    xlabel('Frequency (Hz)')
    ylabel('dB SPL')

    title(sprintf('%s | log1p baseline-corrected tuning heatmap', cluster.Name))

    colormap(redblue)

    fig_title = sprintf('%s_log1p_BaselineCorrected_Frequency_Tuning_Heatmap', cluster.Name);
    safe_name = regexprep(fig_title, '[^\w]', '_');
    save_file = fullfile(save_dir, [safe_name '.pdf']);
    
    exportgraphics(gcf, save_file, 'Resolution', 300)
    close(gcf)

end


%% D' to determine responsiveness of clusters (uses baseline subtraction prior to trial)
sess = S(4);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

response_win = [0 0.2];    % sound window: 0–200 ms
baseline_win = [-0.1 0];   % baseline: 100 ms before sound

dprime_threshold = 1; % Change this threhold accordingly 

Clusternumber = length(sess.Clusters);

sound_responsive_table = table();

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;
    
    onsets = E_freq.OnOffTimes(:,1);
    
    soundFR = nan(length(onsets),1);
    baseFR  = nan(length(onsets),1);
    
    for k = 1:length(onsets)
        
        sound_start = onsets(k) + response_win(1);
        sound_end   = onsets(k) + response_win(2);
        
        base_start = onsets(k) + baseline_win(1);
        base_end   = onsets(k) + baseline_win(2);
        
        sound_spikes = sum(spikes >= sound_start & spikes < sound_end);
        base_spikes  = sum(spikes >= base_start  & spikes < base_end);
        
        soundFR(k) = sound_spikes / diff(response_win);
        baseFR(k)  = base_spikes  / diff(baseline_win);
    end
    
    mu_sound = mean(soundFR, 'omitnan');
    mu_base  = mean(baseFR, 'omitnan');
    
    var_sound = var(soundFR, 'omitnan');
    var_base  = var(baseFR, 'omitnan');
    
    pooled_sd = sqrt(0.5 * (var_sound + var_base));
    
    if pooled_sd == 0 || isnan(pooled_sd)
        dprime = NaN;
    else
        dprime = (mu_sound - mu_base) / pooled_sd;
    end
    
    sound_responsive = abs(dprime) >= dprime_threshold;
    
    if sound_responsive && dprime > 0
        direction = "enhanced";
    elseif sound_responsive && dprime < 0
        direction = "suppressed";
    else
        direction = "nonresponsive";
    end
    
    new_row = table( ...
        string(cluster.Name), ...
        i, ...
        dprime, ...
        mu_sound, ...
        mu_base, ...
        sound_responsive, ...
        direction, ...
        'VariableNames', {'ClusterName','ClusterIndex','dprime','MeanSoundFR','MeanBaselineFR','SoundResponsive','Direction'} ...
    );
    
    sound_responsive_table = [sound_responsive_table; new_row];
end

sound_responsive_table


%% D' to determine responsiveness of clusters using log1p firing rates

sess = S(4);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

response_win = [0 0.2];    % sound window: 0–200 ms
baseline_win = [-0.1 0];   % baseline: 100 ms before sound

dprime_threshold = 1; % Change this threshold accordingly 

Clusternumber = length(sess.Clusters);

sound_responsive_table = table();

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;
    
    onsets = E_freq.OnOffTimes(:,1);
    
    soundFR_log = nan(length(onsets),1);
    baseFR_log  = nan(length(onsets),1);
    
    for k = 1:length(onsets)
        
        sound_start = onsets(k) + response_win(1);
        sound_end   = onsets(k) + response_win(2);
        
        base_start = onsets(k) + baseline_win(1);
        base_end   = onsets(k) + baseline_win(2);
        
        sound_spikes = sum(spikes >= sound_start & spikes < sound_end);
        base_spikes  = sum(spikes >= base_start  & spikes < base_end);
        
        soundFR = sound_spikes / diff(response_win);
        baseFR  = base_spikes  / diff(baseline_win);
        
        % log1p transform
        soundFR_log(k) = log1p(soundFR);
        baseFR_log(k)  = log1p(baseFR);
    end
    
    mu_sound = mean(soundFR_log, 'omitnan');
    mu_base  = mean(baseFR_log, 'omitnan');
    
    var_sound = var(soundFR_log, 'omitnan');
    var_base  = var(baseFR_log, 'omitnan');
    
    pooled_sd = sqrt(0.5 * (var_sound + var_base));
    
    if pooled_sd == 0 || isnan(pooled_sd)
        dprime = NaN;
    else
        dprime = (mu_sound - mu_base) / pooled_sd;
    end
    
    sound_responsive = abs(dprime) >= dprime_threshold;
    
    if sound_responsive && dprime > 0
        direction = "enhanced";
    elseif sound_responsive && dprime < 0
        direction = "suppressed";
    else
        direction = "nonresponsive";
    end
    
    new_row = table( ...
        string(cluster.Name), ...
        i, ...
        dprime, ...
        mu_sound, ...
        mu_base, ...
        sound_responsive, ...
        direction, ...
        'VariableNames', {'ClusterName','ClusterIndex','dprime_log1p','MeanLogSoundFR','MeanLogBaselineFR','SoundResponsive','Direction'} ...
    );
    
    sound_responsive_table = [sound_responsive_table; new_row];
end

sound_responsive_table


%% Compare raw vs log1p d' responsiveness

sess = S(4);

E_freq = sess.find_Event('Freq');

response_win = [0 0.2];
baseline_win = [-0.1 0];

dprime_threshold = 1;

Clusternumber = length(sess.Clusters);

comparison_table = table();

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;
    
    onsets = E_freq.OnOffTimes(:,1);
    
    soundFR = nan(length(onsets),1);
    baseFR  = nan(length(onsets),1);
    
    for k = 1:length(onsets)
        
        stim_start = onsets(k) + response_win(1);
        stim_end   = onsets(k) + response_win(2);
        
        base_start = onsets(k) + baseline_win(1);
        base_end   = onsets(k) + baseline_win(2);
        
        stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
        base_spikes = sum(spikes >= base_start & spikes < base_end);
        
        soundFR(k) = stim_spikes / diff(response_win);
        baseFR(k)  = base_spikes / diff(baseline_win);
    end
    
    %% --- RAW d' ---
    mu_s = mean(soundFR, 'omitnan');
    mu_b = mean(baseFR, 'omitnan');
    
    var_s = var(soundFR, 'omitnan');
    var_b = var(baseFR, 'omitnan');
    
    pooled_sd = sqrt(0.5 * (var_s + var_b));
    
    if pooled_sd == 0 || isnan(pooled_sd)
        d_raw = NaN;
    else
        d_raw = (mu_s - mu_b) / pooled_sd;
    end
    
    %% --- LOG d' ---
    soundFR_log = log1p(soundFR);
    baseFR_log  = log1p(baseFR);
    
    mu_sL = mean(soundFR_log, 'omitnan');
    mu_bL = mean(baseFR_log, 'omitnan');
    
    var_sL = var(soundFR_log, 'omitnan');
    var_bL = var(baseFR_log, 'omitnan');
    
    pooled_sdL = sqrt(0.5 * (var_sL + var_bL));
    
    if pooled_sdL == 0 || isnan(pooled_sdL)
        d_log = NaN;
    else
        d_log = (mu_sL - mu_bL) / pooled_sdL;
    end
    
    %% --- Classification ---
    class_raw = classify_dprime(d_raw, dprime_threshold);
    class_log = classify_dprime(d_log, dprime_threshold);
    
    changed = class_raw ~= class_log;
    
    new_row = table( ...
        string(cluster.Name), ...
        i, ...
        d_raw, ...
        d_log, ...
        class_raw, ...
        class_log, ...
        changed, ...
        'VariableNames', {'ClusterName','ClusterIndex','d_raw','d_log','Class_raw','Class_log','Changed'} ...
    );
    
    comparison_table = [comparison_table; new_row];
end

comparison_table


% How many neurons changed
sum(comparison_table.Changed)

% Changed Neurons only
comparison_table(comparison_table.Changed == true, :)

%Visualize differences
scatter(comparison_table.d_raw, comparison_table.d_log)
xlabel('d'' raw')
ylabel('d'' log1p')
refline(1,0)



%% Find Best Frequency (BF) for all clusters (uses log1p basline corrected firing rate or just baseline subtracted)

sess = S(4);
use_log = true;   % true = log1p, false = raw baseline subtraction

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

response_win = [0 0.2];      % 0–200 ms after tone onset
baseline_win = [-0.1 0];     % 100 ms before tone onset

freq_levels = sort(unique(E_freq.Values), 'ascend');
db_levels   = sort(unique(E_db.Values), 'descend');

n_freq = length(freq_levels);
n_db   = length(db_levels);

Clusternumber = length(sess.Clusters);

BF_table = table();

for i = 1:Clusternumber
    
    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for d = 1:n_db
        
        target_db = db_levels(d);
        
        for f = 1:n_freq
            
            target_freq = freq_levels(f);

            trials = E_freq.Values == target_freq & E_db.Values == target_db;
            onsets = E_freq.OnOffTimes(trials, 1);

            trial_vals = [];

            for k = 1:length(onsets)
            
                stim_start = onsets(k) + response_win(1);
                stim_end   = onsets(k) + response_win(2);
            
                base_start = onsets(k) + baseline_win(1);
                base_end   = onsets(k) + baseline_win(2);
            
                stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
                base_spikes = sum(spikes >= base_start & spikes < base_end);
            
                stim_FR = stim_spikes / diff(response_win);
                base_FR = base_spikes / diff(baseline_win);
            
                if use_log
                    trial_vals(end+1) = log1p(stim_FR) - log1p(base_FR);
                else
                    trial_vals(end+1) = stim_FR - base_FR;
                end
            end
            
            FR_map(d,f) = mean(trial_vals, 'omitnan');
        end
    end

    % Find maximum response in full Freq x dB map
    [max_response, linear_idx] = max(FR_map(:), [], 'omitnan');
    [best_db_idx, best_freq_idx] = ind2sub(size(FR_map), linear_idx);

    BF = freq_levels(best_freq_idx);
    BF_dBSPL = db_levels(best_db_idx);

    new_row = table( ...
        string(cluster.Name), ...
        i, ...
        BF, ...
        BF_dBSPL, ...
        max_response, ...
        'VariableNames', {'ClusterName','ClusterIndex','BestFrequency_Hz','BestFrequency_dBSPL','MaxLog1pBaselineCorrectedResponse'} ...
    );

    BF_table = [BF_table; new_row];
end

BF_table
 % BestFrequency_Hz = the frequency where that unit responds most strongly, and BestFrequency_dBSPL = the sound level where that max occurred.



%% Bandwidth from frequency tuning heatmaps "bandwidth at 50% maximum"

sess = S(4);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

response_win = [0 0.2];
baseline_win = [-0.1 0];

use_log = true;   % true = log1p(stimFR)-log1p(baseFR), false = stimFR-baseFR

freq_levels = sort(unique(E_freq.Values), 'ascend');
db_levels   = sort(unique(E_db.Values), 'descend');

n_freq = length(freq_levels);
n_db   = length(db_levels);

Clusternumber = length(sess.Clusters);

bandwidth_table = table();

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for d = 1:n_db
        for f = 1:n_freq

            trials = E_freq.Values == freq_levels(f) & E_db.Values == db_levels(d);
            onsets = E_freq.OnOffTimes(trials, 1);

            trial_vals = [];

            for k = 1:length(onsets)

                stim_start = onsets(k) + response_win(1);
                stim_end   = onsets(k) + response_win(2);

                base_start = onsets(k) + baseline_win(1);
                base_end   = onsets(k) + baseline_win(2);

                stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
                base_spikes = sum(spikes >= base_start & spikes < base_end);

                stim_FR = stim_spikes / diff(response_win);
                base_FR = base_spikes / diff(baseline_win);

                if use_log
                    trial_vals(end+1) = log1p(stim_FR) - log1p(base_FR);
                else
                    trial_vals(end+1) = stim_FR - base_FR;
                end
            end

            FR_map(d,f) = mean(trial_vals, 'omitnan');
        end
    end

    % Best frequency = max response across full map
    [max_response, linear_idx] = max(FR_map(:), [], 'omitnan');
    [best_db_idx, best_freq_idx] = ind2sub(size(FR_map), linear_idx);

    BF = freq_levels(best_freq_idx);
    BF_dBSPL = db_levels(best_db_idx);

    % Use tuning curve at BF dB level
    tuning_curve = FR_map(best_db_idx, :);

    % Half-max threshold
    halfmax_thresh = 0.5 * max_response; % all frequencies where the neuron's response is at least half of its max peak:

    active_freqs = freq_levels(tuning_curve >= halfmax_thresh);

    if isempty(active_freqs) || numel(active_freqs) < 2
        bandwidth_oct = 0;
        bandwidth_low_Hz = BF;
        bandwidth_high_Hz = BF;
    else
        bandwidth_low_Hz = min(active_freqs);
        bandwidth_high_Hz = max(active_freqs);
        bandwidth_oct = log2(bandwidth_high_Hz / bandwidth_low_Hz);
    end

    new_row = table( ...
        string(cluster.Name), ...
        i, ...
        BF, ...
        BF_dBSPL, ...
        max_response, ...
        bandwidth_low_Hz, ...
        bandwidth_high_Hz, ...
        bandwidth_oct, ...
        'VariableNames', {'ClusterName','ClusterIndex','BF_Hz','BF_dBSPL','MaxResponse','BandwidthLow_Hz','BandwidthHigh_Hz','Bandwidth_octaves'} ...
    );

    bandwidth_table = [bandwidth_table; new_row];

end

bandwidth_table


% 0 octaves       = only one tested frequency above half-max
% ~1 octave       = response spans a 2-fold frequency range
% ~2 octaves      = response spans a 4-fold frequency range
% larger values   = broader tuning

%% Q20 bandwidth using d' >= 1 as response criterion

sess = S(4);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

response_win = [0 0.2];
baseline_win = [-0.1 0];

use_log = true;          % true = log1p FR d', false = raw FR d'
dprime_criterion = 1;    % response criterion

freq_levels = sort(unique(E_freq.Values), 'ascend');
db_levels   = sort(unique(E_db.Values), 'descend');

n_freq = length(freq_levels);
n_db   = length(db_levels);

Clusternumber = length(sess.Clusters);

Q20_table = table();

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);
    dprime_map = nan(n_db, n_freq);

    for d = 1:n_db
        for f = 1:n_freq

            trials = E_freq.Values == freq_levels(f) & E_db.Values == db_levels(d);
            onsets = E_freq.OnOffTimes(trials, 1);

            stim_vals = [];
            base_vals = [];

            for k = 1:length(onsets)

                stim_start = onsets(k) + response_win(1);
                stim_end   = onsets(k) + response_win(2);

                base_start = onsets(k) + baseline_win(1);
                base_end   = onsets(k) + baseline_win(2);

                stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
                base_spikes = sum(spikes >= base_start & spikes < base_end);

                stim_FR = stim_spikes / diff(response_win);
                base_FR = base_spikes / diff(baseline_win);

                if use_log
                    stim_vals(end+1) = log1p(stim_FR);
                    base_vals(end+1) = log1p(base_FR);
                else
                    stim_vals(end+1) = stim_FR;
                    base_vals(end+1) = base_FR;
                end
            end

            % Response magnitude map
            FR_map(d,f) = mean(stim_vals - base_vals, 'omitnan');

            % d' map
            mu_stim = mean(stim_vals, 'omitnan');
            mu_base = mean(base_vals, 'omitnan');

            var_stim = var(stim_vals, 'omitnan');
            var_base = var(base_vals, 'omitnan');

            pooled_sd = sqrt(0.5 * (var_stim + var_base));

            if pooled_sd == 0 || isnan(pooled_sd)
                dprime_map(d,f) = NaN;
            else
                dprime_map(d,f) = (mu_stim - mu_base) / pooled_sd;
            end
        end
    end

    %% BF = frequency with maximum response magnitude
    [max_response, linear_idx] = max(FR_map(:), [], 'omitnan');
    [best_db_idx, best_freq_idx] = ind2sub(size(FR_map), linear_idx);

    BF = freq_levels(best_freq_idx);
    BF_dBSPL = db_levels(best_db_idx);

    %% Threshold = lowest dB at BF where d' >= 1
    db_ascending = sort(db_levels, 'ascend');

    threshold_dB = NaN;

    for dd = 1:length(db_ascending)

        this_db = db_ascending(dd);
        row_idx = find(db_levels == this_db);

        if dprime_map(row_idx, best_freq_idx) >= dprime_criterion
            threshold_dB = this_db;
            break
        end
    end

    %% Q20 = threshold + 20 dB
    if isnan(threshold_dB)

        Q20_dB = NaN;
        Q20_Low_Hz = NaN;
        Q20_High_Hz = NaN;
        Q20_Bandwidth_octaves = NaN;

    else

        Q20_target = threshold_dB + 20;

        available_above = db_levels(db_levels >= Q20_target);

        if isempty(available_above)

            Q20_dB = NaN;
            Q20_Low_Hz = NaN;
            Q20_High_Hz = NaN;
            Q20_Bandwidth_octaves = NaN;

        else

            Q20_dB = min(available_above);
            q20_row_idx = find(db_levels == Q20_dB);

            responsive_freqs = freq_levels(dprime_map(q20_row_idx,:) >= dprime_criterion);

            if numel(responsive_freqs) < 2
                Q20_Low_Hz = BF;
                Q20_High_Hz = BF;
                Q20_Bandwidth_octaves = 0;
            else
                Q20_Low_Hz = min(responsive_freqs);
                Q20_High_Hz = max(responsive_freqs);
                Q20_Bandwidth_octaves = log2(Q20_High_Hz / Q20_Low_Hz);
            end
        end
    end

    new_row = table( ...
        string(cluster.Name), ...
        i, ...
        BF, ...
        BF_dBSPL, ...
        max_response, ...
        threshold_dB, ...
        Q20_dB, ...
        Q20_Low_Hz, ...
        Q20_High_Hz, ...
        Q20_Bandwidth_octaves, ...
        'VariableNames', {'ClusterName','ClusterIndex','BF_Hz','BF_dBSPL','MaxResponse','Threshold_dBSPL','Q20_dBSPL','Q20_Low_Hz','Q20_High_Hz','Q20_Bandwidth_octaves'} ...
    );

    Q20_table = [Q20_table; new_row];

end

Q20_table


%% Plot heatmap for one cluster with Q20 bandwidth overlay

cluster_to_plot = 15;   % <-- change this to the ClusterIndex you want

sess = S(4);

E_freq = sess.find_Event('Freq');
E_db   = sess.find_Event('dBSPL');

response_win = [0 0.2];
baseline_win = [-0.1 0];
use_log = true;

freq_levels = sort(unique(E_freq.Values), 'ascend');
db_levels   = sort(unique(E_db.Values), 'descend');

n_freq = length(freq_levels);
n_db   = length(db_levels);

% Get row from Q20 table
row = Q20_table(Q20_table.ClusterIndex == cluster_to_plot, :);

cluster = sess.Clusters(cluster_to_plot);
spikes = cluster.SpikeTimes;

FR_map = nan(n_db, n_freq);

for d = 1:n_db
    
    target_db = db_levels(d);
    
    for f = 1:n_freq
        
        target_freq = freq_levels(f);
        
        trials = E_freq.Values == target_freq & E_db.Values == target_db;
        onsets = E_freq.OnOffTimes(trials, 1);
        
        trial_vals = [];
        
        for k = 1:length(onsets)
            
            stim_start = onsets(k) + response_win(1);
            stim_end   = onsets(k) + response_win(2);
            
            base_start = onsets(k) + baseline_win(1);
            base_end   = onsets(k) + baseline_win(2);
            
            stim_spikes = sum(spikes >= stim_start & spikes < stim_end);
            base_spikes = sum(spikes >= base_start & spikes < base_end);
            
            stim_FR = stim_spikes / diff(response_win);
            base_FR = base_spikes / diff(baseline_win);
            
            if use_log
                trial_vals(end+1) = log1p(stim_FR) - log1p(base_FR);
            else
                trial_vals(end+1) = stim_FR - base_FR;
            end
        end
        
        FR_map(d,f) = mean(trial_vals, 'omitnan');
    end
end

% Plot heatmap

figure
imagesc(freq_levels, 1:n_db, FR_map)

set(gca, 'XScale', 'log')
set(gca, 'YDir', 'reverse')

yticks(1:n_db)
yticklabels(string(db_levels))

xticks(freq_levels)
xtickangle(45)

xlabel('Frequency (Hz)')
ylabel('dB SPL')

cb = colorbar;
ylabel(cb, 'log1p(stim FR) - log1p(baseline FR)')

colormap(redblue)

title(sprintf('%s | Q20 bandwidth = %.2f octaves', ...
    cluster.Name, row.Q20_Bandwidth_octaves))

hold on

% Overlay Q20 bandwidth

q20_db = row.Q20_dBSPL;
low_hz = row.Q20_Low_Hz;
high_hz = row.Q20_High_Hz;
bf_hz = row.BF_Hz;

q20_row = find(db_levels == q20_db);

% horizontal Q20 bandwidth line
plot([low_hz high_hz], [q20_row q20_row], 'c-', 'LineWidth', 4)

% mark low/high edges
plot(low_hz, q20_row, 'co', 'MarkerFaceColor', 'c')
plot(high_hz, q20_row, 'co', 'MarkerFaceColor', 'c')

% mark BF
plot(bf_hz, q20_row, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 8)

legend({'Q20 bandwidth','Q20 edges','BF'}, 'Location', 'bestoutside')

% Restrict to Q20 frequency range

low_hz  = row.Q20_Low_Hz;
high_hz = row.Q20_High_Hz;

freq_mask = freq_levels >= low_hz & freq_levels <= high_hz;

freq_subset = freq_levels(freq_mask);
FR_subset   = FR_map(:, freq_mask);

% Plot restricted heatmap

figure
imagesc(freq_subset, 1:n_db, FR_subset)

set(gca, 'XScale', 'log')
set(gca, 'YDir', 'reverse')

yticks(1:n_db)
yticklabels(string(db_levels))

xticks(freq_subset)
xtickangle(45)

xlabel('Frequency (Hz)')
ylabel('dB SPL')

cb = colorbar;
ylabel(cb, 'log1p(stim FR) - log1p(baseline FR)')

colormap(redblue)

title(sprintf('%s | Q20 BW = %.2f oct (%.0f–%.0f Hz)', ...
    cluster.Name, ...
    row.Q20_Bandwidth_octaves, ...
    low_hz, high_hz))

hold on

% Overlay Q20 bandwidth

q20_db = row.Q20_dBSPL;
low_hz = row.Q20_Low_Hz;
high_hz = row.Q20_High_Hz;
bf_hz = row.BF_Hz;

q20_row = find(db_levels == q20_db);

% horizontal Q20 bandwidth line
plot([low_hz high_hz], [q20_row q20_row], 'c-', 'LineWidth', 4)

% mark low/high edges
plot(low_hz, q20_row, 'co', 'MarkerFaceColor', 'c')
plot(high_hz, q20_row, 'co', 'MarkerFaceColor', 'c')

% mark BF
plot(bf_hz, q20_row, 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 8)

legend({'Q20 bandwidth','Q20 edges','BF'}, 'Location', 'bestoutside')










