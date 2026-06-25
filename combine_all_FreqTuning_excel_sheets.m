function combine_all_FreqTuning_excel_sheets(subjectPath)
sheetsToCombine = {
    'Dprime_raw'
    'Dprime_log1p'
    'Dprime_comparison'
    'BestFrequency'
    'Bandwidth'
    'Q20'
};

concatDirs = dir(fullfile(subjectPath, '*_concat'));
concatDirs = concatDirs([concatDirs.isdir]);

% combinedExcelFile = fullfile(subjectPath, 'SUBJ-ID-1180_combined_freqTuning_tables.xlsx');
[~, subjectID] = fileparts(subjectPath);

combinedExcelFile = fullfile( ...
    subjectPath, ...
    [subjectID '_combined_freqTuning_tables.xlsx']);

for s = 1:length(sheetsToCombine)

    sheetName = sheetsToCombine{s};
    combinedTable = table();

    for c = 1:length(concatDirs)

        concatName = concatDirs(c).name;
        concatPath = fullfile(concatDirs(c).folder, concatName);

        excelFile = fullfile(concatPath, [concatName '_freqTuning_tables.xlsx']);

        if ~exist(excelFile, 'file')
            warning('Missing Excel file: %s', excelFile)
            continue
        end

        try
            T = readtable(excelFile, 'Sheet', sheetName);
        catch
            warning('Missing sheet %s in %s', sheetName, excelFile)
            continue
        end

        % Add concat number as first column
        ConcatNumber = repmat(string(concatName), height(T), 1);
        T = addvars(T, ConcatNumber, 'Before', 1);

        combinedTable = [combinedTable; T];
    end

    if ~isempty(combinedTable)
        writetable(combinedTable, combinedExcelFile, 'Sheet', sheetName);
    end
end

disp("Saved combined Excel workbook:")
disp(combinedExcelFile)
end