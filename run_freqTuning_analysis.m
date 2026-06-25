function run_freqTuning_analysis(subjectPath)
% ============================================================
% Run frequency tuning analysis across all *_concat folders
%
% Example:
%
% run_freqTuning_analysis( ...
% '/mnt/CL_8TB_4/Mason/Ephys/DorsalAUD_PL/Sorting/SUBJ-ID-1180')
%
% ============================================================

tic

clearvars -except subjectPath

concatDirs = dir(fullfile(subjectPath, '*_concat'));
concatDirs = concatDirs([concatDirs.isdir]);

for c = 1:length(concatDirs)

    clear S D

    concatName = concatDirs(c).name;
    savePath = fullfile(concatDirs(c).folder, concatName);
    infoPath = savePath;
    DataPath = fullfile(savePath, 'kilosort4');

    fprintf('\n\n==============================\n')
    fprintf('Running concat folder: %s\n', savePath)
    fprintf('==============================\n')

    excelFile = fullfile(savePath, [concatName '_freqTuning_tables.xlsx']);

    if ~exist(DataPath, 'dir')
        warning('No kilosort4 folder found for %s. Skipping.', concatName)
        continue
    end

    try
        datFile = fix_dat_path(DataPath);
        S = epa.load.phy(DataPath);
    catch ME
        warning('Could not load phy data for %s. Skipping.\n%s', concatName, ME.message)
        continue
    end

    % if length(S) < 4
    %     warning('%s has fewer than 4 sessions. Skipping.', concatName)
    %     continue
    % end
    % 
    % freqSessionName = string(S(4).Name);
%% Find frequency tuning session dynamically

sessionNames = string({S.Name});

freqSessionIdx = find(contains(sessionNames, 'freqTuning', 'IgnoreCase', true), 1);

if isempty(freqSessionIdx)
    warning('No freqTuning session found for %s. Skipping.', concatName)
    disp("Available sessions:")
    disp(sessionNames')
    continue
end

freqSessionName = string(S(freqSessionIdx).Name);

fprintf('Using freqTuning session index %d: %s\n', ...
    freqSessionIdx, freqSessionName)


    %% ============================================================
    %% Find ORIGINAL trialInfo file
    %% ============================================================

    d = dir(fullfile(infoPath, '**', '*trialInfo.csv'));

    fullnames = string(fullfile({d.folder}, {d.name}));

    d = d(~contains(fullnames, 'EPA_clean_events_freqTuning'));
    d = d(~contains(string({d.name}), 'clean'));
    d = d(~contains(string({d.name}), 'concat'));

    d = d(contains(string({d.name}), freqSessionName) | ...
          contains(string({d.name}), 'freqTuning'));

    if isempty(d)
        disp("No ORIGINAL freqTuning trialInfo found.")
        warning('No original trialInfo file found for %s. Skipping.', concatName)
        continue
    end

    trialInfoFile = fullfile(d(1).folder, d(1).name);

    disp("Using trialInfo file:")
    disp(trialInfoFile)

    T = readtable(trialInfoFile);

    if ~ismember('ComputerTimestamp_6', T.Properties.VariableNames)
        warning('Selected trialInfo file does not contain ComputerTimestamp_6.')
        continue
    end

    %% ============================================================
    %% Build EPA clean events table
    %% ============================================================

    Trial_onset  = double(T.ComputerTimestamp_6);
    Trial_offset = Trial_onset + double(T.StimDuration) ./ 1000;

    T_clean = table();

    T_clean.Trial_onset  = Trial_onset;
    T_clean.Trial_offset = Trial_offset;
    T_clean.dBSPL        = double(T.dBSPL);
    T_clean.Freq         = double(T.Freq);
    T_clean.TrialID      = double(T.TrialID);

    cleanInfoPath = fullfile(infoPath, 'EPA_clean_events_freqTuning');

    if ~exist(cleanInfoPath, 'dir')
        mkdir(cleanInfoPath);
    end

    cleanFile = fullfile( ...
        cleanInfoPath, ...
        freqSessionName + "_trialInfo.csv");

    writetable(T_clean, cleanFile);

    writetable(T_clean, excelFile, 'Sheet', 'T_clean');

    disp("Saved clean file:")
    disp(cleanFile)

    %% ============================================================
    %% Load EPA events
    %% ============================================================

    epa.load.events(S, cleanInfoPath);

    % sess = S(4);
    sess = S(freqSessionIdx);

    E_freq = sess.find_Event('Freq');
    E_db   = sess.find_Event('dBSPL');

    Clusternumber = length(sess.Clusters);

    freq_levels = sort(unique(E_freq.Values), 'ascend');
    db_levels   = sort(unique(E_db.Values), 'descend');

    n_freq = length(freq_levels);
    n_db   = length(db_levels);

    response_win = [0 0.2];
    baseline_win = [-0.1 0];

    %% ============================================================
    %% PSTHs
    %% ============================================================

    make_PSTHs( ...
        sess, ...
        E_freq, ...
        E_db, ...
        freq_levels, ...
        db_levels, ...
        savePath, ...
        excelFile)

    %% ============================================================
    %% Raw baseline-subtracted heatmaps
    %% ============================================================

    make_raw_heatmaps( ...
        sess, ...
        E_freq, ...
        E_db, ...
        freq_levels, ...
        db_levels, ...
        response_win, ...
        baseline_win, ...
        savePath)

    %% ============================================================
    %% log1p heatmaps
    %% ============================================================

    make_log1p_heatmaps( ...
        sess, ...
        E_freq, ...
        E_db, ...
        freq_levels, ...
        db_levels, ...
        response_win, ...
        baseline_win, ...
        savePath)

    %% ============================================================
    %% Dprime raw
    %% ============================================================

    sound_responsive_table = compute_dprime_raw( ...
        sess, ...
        E_freq, ...
        response_win, ...
        baseline_win);

    writetable(sound_responsive_table, ...
        excelFile, ...
        'Sheet', ...
        'Dprime_raw')

    %% ============================================================
    %% Dprime log1p
    %% ============================================================

    sound_responsive_table_log1p = compute_dprime_log1p( ...
        sess, ...
        E_freq, ...
        response_win, ...
        baseline_win);

    writetable(sound_responsive_table_log1p, ...
        excelFile, ...
        'Sheet', ...
        'Dprime_log1p')

    %% ============================================================
    %% Dprime comparison
    %% ============================================================

    comparison_table = compare_dprime_raw_vs_log( ...
        sess, ...
        E_freq, ...
        response_win, ...
        baseline_win, ...
        savePath);

    writetable(comparison_table, ...
        excelFile, ...
        'Sheet', ...
        'Dprime_comparison')

    %% ============================================================
    %% Best Frequency
    %% ============================================================

    BF_table = compute_best_frequency( ...
        sess, ...
        E_freq, ...
        E_db, ...
        freq_levels, ...
        db_levels, ...
        response_win, ...
        baseline_win);

    writetable(BF_table, ...
        excelFile, ...
        'Sheet', ...
        'BestFrequency')

    %% ============================================================
    %% Bandwidth
    %% ============================================================

    bandwidth_table = compute_bandwidth( ...
        sess, ...
        E_freq, ...
        E_db, ...
        freq_levels, ...
        db_levels, ...
        response_win, ...
        baseline_win);

    writetable(bandwidth_table, ...
        excelFile, ...
        'Sheet', ...
        'Bandwidth')

    %% ============================================================
    %% Q20
    %% ============================================================

    Q20_table = compute_Q20( ...
        sess, ...
        E_freq, ...
        E_db, ...
        freq_levels, ...
        db_levels, ...
        response_win, ...
        baseline_win);

    writetable(Q20_table, ...
        excelFile, ...
        'Sheet', ...
        'Q20')

    fprintf('\nFinished concat folder: %s\n', concatName)
    fprintf('Saved Excel workbook:\n%s\n', excelFile)

end

toc

end

%% ============================================================
%% Helper functions
%% ============================================================

function cmap = redblue(m)

if nargin < 1
    m = 256;
end

if mod(m,2) ~= 0
    m = m + 1;
end

r = [(0:m/2-1)'/(m/2); ones(m/2,1)];
g = [(0:m/2-1)'/(m/2); (m/2-1:-1:0)'/(m/2)];
b = [ones(m/2,1); (m/2-1:-1:0)'/(m/2)];

cmap = [r g b];

end


function class_label = classify_dprime_local(dprime, threshold)

if isnan(dprime)
    class_label = "NaN";
elseif abs(dprime) < threshold
    class_label = "nonresponsive";
elseif dprime >= threshold
    class_label = "enhanced";
elseif dprime <= -threshold
    class_label = "suppressed";
else
    class_label = "nonresponsive";
end

end

% makePSTHs
function make_PSTHs(sess, E_freq, E_db, freq_levels, db_levels, savePath, excelFile)

Clusternumber = length(sess.Clusters);

n_freq = length(freq_levels);
n_db   = length(db_levels);

win = [-0.1 0.3];
bin_size = 0.01;
edges = win(1):bin_size:(win(2) + bin_size);
t = edges(1:end-1) + bin_size/2;

save_dir = fullfile(savePath, 'PSTH');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

blank_psth_table = table();

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    all_psth = cell(n_db, n_freq);
    global_max = 0;

    for dB_i = 1:n_db

        target_db = db_levels(dB_i);

        for f_i = 1:n_freq

            target_freq = freq_levels(f_i);

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

            all_psth{dB_i,f_i} = psth;
            global_max = max(global_max, max(psth));
        end
    end

    if global_max == 0
        global_max = 1;
    end

    figure(1)
    clf
    set(gcf, 'Position', [100 100 1800 1000])

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            subplot(n_db, n_freq, (dB_i-1)*n_freq + f_i)

            psth = all_psth{dB_i,f_i};

            bar(t, psth, 1, 'FaceColor', 'k', 'EdgeColor', 'none')

            xlim(win)
            ylim([0 global_max * 1.2])

            if dB_i == 1
                title(sprintf('%d Hz', freq_levels(f_i)))
            end

            if f_i == 1
                ylabel(sprintf('%d dB\nFR (Hz)', db_levels(dB_i)))
            end

            if dB_i == n_db
                xlabel('Time (s)')
            end
        end
    end

    fig_title = sprintf('%s | dB x Frequency PSTH', cluster.Name);
    sgtitle(fig_title)

    safe_name = regexprep(fig_title, '[^\w]', '_');
    save_file = fullfile(save_dir, [safe_name '.pdf']);

    exportgraphics(gcf, save_file, 'Resolution', 300)
    close(gcf)

    total_spikes_session = length(spikes);
    total_spikes_psth = 0;
    onsets_all = E_freq.OnOffTimes(:,1);

    for k = 1:length(onsets_all)
        total_spikes_psth = total_spikes_psth + ...
            sum(spikes >= onsets_all(k)+win(1) & spikes < onsets_all(k)+win(2));
    end

    if total_spikes_psth == 0

        new_row = table( ...
            string(cluster.Name), ...
            i, ...
            cluster.ID, ...
            total_spikes_session, ...
            total_spikes_psth, ...
            'VariableNames', {'ClusterName','ClusterIndex','ClusterID','SessionSpikes','PSTHSpikes'} ...
        );

        blank_psth_table = [blank_psth_table; new_row];

        fprintf('Blank PSTH: Cluster index %d | ID %d | %s | session spikes = %d | PSTH spikes = %d\n', ...
            i, cluster.ID, cluster.Name, total_spikes_session, total_spikes_psth);
    end
end

writetable(blank_psth_table, excelFile, 'Sheet', 'Blank_PSTHs');

end


function make_raw_heatmaps(sess, E_freq, E_db, freq_levels, db_levels, response_win, baseline_win, savePath)

Clusternumber = length(sess.Clusters);
n_freq = length(freq_levels);
n_db   = length(db_levels);

save_dir = fullfile(savePath, 'TuningHeatmaps');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

global_min = inf;
global_max = -inf;

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            trials = E_freq.Values == freq_levels(f_i) & E_db.Values == db_levels(dB_i);
            onsets = E_freq.OnOffTimes(trials, 1);

            stim_FRs = [];
            base_FRs = [];

            for k = 1:length(onsets)

                stim_spikes = sum(spikes >= onsets(k)+response_win(1) & spikes < onsets(k)+response_win(2));
                base_spikes = sum(spikes >= onsets(k)+baseline_win(1) & spikes < onsets(k)+baseline_win(2));

                stim_FRs(end+1) = stim_spikes / diff(response_win);
                base_FRs(end+1) = base_spikes / diff(baseline_win);
            end

            FR_map(dB_i,f_i) = mean(stim_FRs - base_FRs, 'omitnan');
        end
    end

    global_min = min(global_min, min(FR_map(:), [], 'omitnan'));
    global_max = max(global_max, max(FR_map(:), [], 'omitnan'));
end

max_abs = max(abs([global_min global_max]));

if isempty(max_abs) || isnan(max_abs) || max_abs == 0
    max_abs = 1;
end

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            trials = E_freq.Values == freq_levels(f_i) & E_db.Values == db_levels(dB_i);
            onsets = E_freq.OnOffTimes(trials, 1);

            stim_FRs = [];
            base_FRs = [];

            for k = 1:length(onsets)

                stim_spikes = sum(spikes >= onsets(k)+response_win(1) & spikes < onsets(k)+response_win(2));
                base_spikes = sum(spikes >= onsets(k)+baseline_win(1) & spikes < onsets(k)+baseline_win(2));

                stim_FRs(end+1) = stim_spikes / diff(response_win);
                base_FRs(end+1) = base_spikes / diff(baseline_win);
            end

            FR_map(dB_i,f_i) = mean(stim_FRs - base_FRs, 'omitnan');
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

    colorbar
    xlabel('Frequency (Hz)')
    ylabel('dB SPL')

    title(sprintf('%s | Frequency tuning heatmap', cluster.Name))

    colormap(redblue)

    fig_title = sprintf('%s_Frequency_Tuning_Heatmap', cluster.Name);
    safe_name = regexprep(fig_title, '[^\w]', '_');
    save_file = fullfile(save_dir, [safe_name '.pdf']);

    exportgraphics(gcf, save_file, 'Resolution', 300)
    close(gcf)
end

end

function make_log1p_heatmaps(sess, E_freq, E_db, freq_levels, db_levels, response_win, baseline_win, savePath)

Clusternumber = length(sess.Clusters);
n_freq = length(freq_levels);
n_db   = length(db_levels);

save_dir = fullfile(savePath, 'TuningHeatmaps_log1p');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

global_min = inf;
global_max = -inf;

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            trials = E_freq.Values == freq_levels(f_i) & E_db.Values == db_levels(dB_i);
            onsets = E_freq.OnOffTimes(trials, 1);

            trial_vals = [];

            for k = 1:length(onsets)

                stim_spikes = sum(spikes >= onsets(k)+response_win(1) & spikes < onsets(k)+response_win(2));
                base_spikes = sum(spikes >= onsets(k)+baseline_win(1) & spikes < onsets(k)+baseline_win(2));

                stim_FR = stim_spikes / diff(response_win);
                base_FR = base_spikes / diff(baseline_win);

                trial_vals(end+1) = log1p(stim_FR) - log1p(base_FR);
            end

            FR_map(dB_i,f_i) = mean(trial_vals, 'omitnan');
        end
    end

    global_min = min(global_min, min(FR_map(:), [], 'omitnan'));
    global_max = max(global_max, max(FR_map(:), [], 'omitnan'));
end

max_abs = max(abs([global_min global_max]));

if isempty(max_abs) || isnan(max_abs) || max_abs == 0
    max_abs = 1;
end

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            trials = E_freq.Values == freq_levels(f_i) & E_db.Values == db_levels(dB_i);
            onsets = E_freq.OnOffTimes(trials, 1);

            trial_vals = [];

            for k = 1:length(onsets)

                stim_spikes = sum(spikes >= onsets(k)+response_win(1) & spikes < onsets(k)+response_win(2));
                base_spikes = sum(spikes >= onsets(k)+baseline_win(1) & spikes < onsets(k)+baseline_win(2));

                stim_FR = stim_spikes / diff(response_win);
                base_FR = base_spikes / diff(baseline_win);

                trial_vals(end+1) = log1p(stim_FR) - log1p(base_FR);
            end

            FR_map(dB_i,f_i) = mean(trial_vals, 'omitnan');
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

end


function sound_responsive_table = compute_dprime_raw(sess, E_freq, response_win, baseline_win)

dprime_threshold = 1;
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

end


function sound_responsive_table_log1p = compute_dprime_log1p(sess, E_freq, response_win, baseline_win)

dprime_threshold = 1;
Clusternumber = length(sess.Clusters);

sound_responsive_table_log1p = table();

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

    sound_responsive_table_log1p = [sound_responsive_table_log1p; new_row];
end

end


function comparison_table = compare_dprime_raw_vs_log(sess, E_freq, response_win, baseline_win, savePath)

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

    %% Raw d'
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

    %% Log1p d'
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

    %% Classification
    class_raw = classify_dprime_local(d_raw, dprime_threshold);
    class_log = classify_dprime_local(d_log, dprime_threshold);

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

%% Save comparison plot
save_dir = fullfile(savePath, 'DprimeComparison');

if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

figure(1)
clf

scatter(comparison_table.d_raw, comparison_table.d_log)

xlabel('d'' raw')
ylabel('d'' log1p')

refline(1,0)

exportgraphics(gcf, fullfile(save_dir, 'raw_vs_log1p_dprime.pdf'), 'Resolution', 300)

close(gcf)

end


function BF_table = compute_best_frequency(sess, E_freq, E_db, freq_levels, db_levels, response_win, baseline_win)

use_log = true;

Clusternumber = length(sess.Clusters);

n_freq = length(freq_levels);
n_db   = length(db_levels);

BF_table = table();

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            trials = E_freq.Values == freq_levels(f_i) & ...
                     E_db.Values   == db_levels(dB_i);

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

            FR_map(dB_i,f_i) = mean(trial_vals, 'omitnan');
        end
    end

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

end


function bandwidth_table = compute_bandwidth(sess, E_freq, E_db, freq_levels, db_levels, response_win, baseline_win)

use_log = true;

Clusternumber = length(sess.Clusters);
n_freq = length(freq_levels);
n_db   = length(db_levels);

bandwidth_table = table();

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            trials = E_freq.Values == freq_levels(f_i) & ...
                     E_db.Values   == db_levels(dB_i);

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

            FR_map(dB_i,f_i) = mean(trial_vals, 'omitnan');
        end
    end

    [max_response, linear_idx] = max(FR_map(:), [], 'omitnan');
    [best_db_idx, best_freq_idx] = ind2sub(size(FR_map), linear_idx);

    BF = freq_levels(best_freq_idx);
    BF_dBSPL = db_levels(best_db_idx);

    tuning_curve = FR_map(best_db_idx, :);

    halfmax_thresh = 0.5 * max_response;
    active_freqs = freq_levels(tuning_curve >= halfmax_thresh);

    if isempty(active_freqs) || numel(active_freqs) < 2
        bandwidth_low_Hz = BF;
        bandwidth_high_Hz = BF;
        bandwidth_oct = 0;
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

end


function Q20_table = compute_Q20(sess, E_freq, E_db, freq_levels, db_levels, response_win, baseline_win)

use_log = true;
dprime_criterion = 1;

Clusternumber = length(sess.Clusters);
n_freq = length(freq_levels);
n_db   = length(db_levels);

Q20_table = table();

for i = 1:Clusternumber

    cluster = sess.Clusters(i);
    spikes = cluster.SpikeTimes;

    FR_map = nan(n_db, n_freq);
    dprime_map = nan(n_db, n_freq);

    for dB_i = 1:n_db
        for f_i = 1:n_freq

            trials = E_freq.Values == freq_levels(f_i) & ...
                     E_db.Values   == db_levels(dB_i);

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

            FR_map(dB_i,f_i) = mean(stim_vals - base_vals, 'omitnan');

            mu_stim = mean(stim_vals, 'omitnan');
            mu_base = mean(base_vals, 'omitnan');

            var_stim = var(stim_vals, 'omitnan');
            var_base = var(base_vals, 'omitnan');

            pooled_sd = sqrt(0.5 * (var_stim + var_base));

            if pooled_sd == 0 || isnan(pooled_sd)
                dprime_map(dB_i,f_i) = NaN;
            else
                dprime_map(dB_i,f_i) = (mu_stim - mu_base) / pooled_sd;
            end
        end
    end

    [max_response, linear_idx] = max(FR_map(:), [], 'omitnan');
    [best_db_idx, best_freq_idx] = ind2sub(size(FR_map), linear_idx);

    BF = freq_levels(best_freq_idx);
    BF_dBSPL = db_levels(best_db_idx);

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

end


