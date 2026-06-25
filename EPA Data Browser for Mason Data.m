%% Clean trialInfo using session-matching filename (Works for active)

clear S D

infoPath = 'G:\My Drive\Data\Dorsal AUD Recordings\FrequenyTuningTests\Data\260331_concat\Info files';
DataPath = 'G:\My Drive\Data\Dorsal AUD Recordings\FrequenyTuningTests\Data\260331_concat\kilosort4';
datFile = fix_dat_path(DataPath);

d = dir(fullfile(infoPath, '*trialInfo.csv'));
d = d(~contains({d.name}, 'clean'));
d = d(~contains({d.name}, 'concat'));

trialInfoFile = fullfile(d(1).folder, d(1).name);
T = readtable(trialInfoFile);

T_clean = T(:, {
    'Trial_onset', ...
    'Trial_offset', ...
    'AMdepth', ...
    'dBSPL', ...
    'TrialID', ...
    'TrialType', ...
    'Hit', ...
    'Miss', ...
    'CR', ...
    'FA'
});

T_clean.Trial_onset  = double(T_clean.Trial_onset);
T_clean.Trial_offset = double(T_clean.Trial_offset);

cleanInfoPath = fullfile(infoPath, 'EPA_clean_events');
if ~exist(cleanInfoPath, 'dir')
    mkdir(cleanInfoPath);
end

% IMPORTANT: filename must contain session name
sessionName = string(T.Session_id{1});
cleanFile = fullfile(cleanInfoPath, sessionName + "_trialInfo.csv");

writetable(T_clean, cleanFile);

disp("Saved clean event file:")
disp(cleanFile)

%% Load

% datFile = fix_dat_path(DataPath);

S = epa.load.phy(DataPath); 
epa.load.events(S, cleanInfoPath);

%% Check correct session

for i = 1:length(S)
    fprintf('%d: %s\n', i, S(i).Name)
end
%% Verify event timing

sess = S(2);  % change to selected session if needed

E = sess.find_Event('dBSPL');
disp(string(E.Name))
disp(size(E.OnOffTimes))
disp(size(E.Values))


E2 = sess.find_Event('Freq');
disp(string(E2.Name))
disp(size(E2.OnOffTimes))
disp(size(E2.Values))
%% Open DataBrowser

D = epa.DataBrowser;