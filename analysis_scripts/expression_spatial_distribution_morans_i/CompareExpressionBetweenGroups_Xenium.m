%% Compare gene expression between WT and 5xFAD mice by hippocampal region
% Each CSV contains transcript-level rows with at least:
%   gene
%   region_name
%
% Optional:
%   cell_type
%
% The script:
%   1. Reads CSV files for WT and 5xFAD mice
%   2. Standardizes gene names before counting
%   3. Merges known gene aliases before counting
%   4. Removes negative-control and excluded genes before counting
%   5. Standardizes hippocampal region names before counting
%   6. Excludes unwanted regions before counting
%   7. Merges multiple CSVs that belong to the same mouse
%   8. Counts expression per gene and hippocampal region
%   9. Fills missing gene-region-mouse combinations with zero
%   10. Applies normalization for non-rnaseqde tests: TMM, quantile, or none
%   11. Optionally filters genes before testing by:
%        - expression quantile
%        - fold change
%        - minimum mean expression
%        - sample support
%        - variance across samples
%        - effect above approximate Poisson noise
%   12. Compares 5xFAD vs WT within each hippocampal region
%   13. Calculates P values using rnaseqde, Welch t-test, or permutation testing
%   14. Applies FDR correction

clear; clc; close all;

%% =============================
% USER SETTINGS
% =============================

scriptDir = fileparts(mfilename('fullpath'));
inputFolder  = fullfile(scriptDir, 'data', 'xenium');   % edit if needed
outputFolder = fullfile(scriptDir, 'results', 'xenium_expression');

% false = compare per gene and hippocampal region
% true  = compare per gene, hippocampal region, and cell type
useCellTypes = false;

%% =============================
% Normalization settings for filtering / t-test / permutation
% =============================

% Options:
% "TMM"      = TMM normalization across all mice
% "quantile" = quantile normalization across all mice
% "none"     = use raw counts
%
% Note:
% rnaseqde always uses raw integer counts internally, regardless of this flag.
normalizationMethod = "TMM";

% TMM parameters, used only if normalizationMethod = "TMM"
Mtrim = 0.30;
Atrim = 0.05;

%% =============================
% Statistical test settings
% =============================

% Options:
% "rnaseqde"                  = negative-binomial RNA-seq DE test on raw integer counts
% "welch_ttest_log2_plus1"    = Welch t-test on log2(normalized expression + 1)
% "welch_ttest_raw"           = Welch t-test on normalized expression values
% "permutation"               = permutation test on normalized expression values
statisticalTest = "rnaseqde";

% rnaseqde settings, used only if statisticalTest = "rnaseqde"
rnaseqdeVarianceLink = "local"; % "constant"
rnaseqdeFDRMethod = "bh";

% Permutation settings, used only if statisticalTest = "permutation"
number_of_realizations = 1000;
useExhaustivePermutations = true;

% P-value method for permutation only:
% "normal_fit" = fit null permutation distribution using mean and SD
% "empirical"  = direct empirical two-sided permutation P value
pValueMode = "normal_fit";

%% =============================
% Reporting thresholds
% =============================

cutoff_p = 0.05;
cutoff_FDR = 0.15;

% If true, significant table requires only P < cutoff_p.
% If false, significant table requires only FDR < cutoff_FDR.
% FDR is calculated globally by bhFDR.
applyPCutoffForReporting = false;

%% =============================
% Pre-testing filtering settings
% =============================

useExpressionQuantileFiltering = false;
expression_quantile_cutoff = 0.5;

useFoldChangeFiltering = true;
foldChangeCutoff = 1.1;

useMinimumMeanExpressionFiltering = true;
minimumMeanExpression = 100;

useSampleSupportFiltering = false;
minimumNumberOfExpressingSamples = 3;
sampleSupportExpressionThreshold = 100;

useVarianceFiltering = false;
minimumVarianceAcrossSamples = 1;

% Approximate Poisson-noise-aware effect-size filter:
%   expected_noise = sqrt(mean_WT + mean_5xFAD)
%   keep if abs(mean_5xFAD - mean_WT) >= poissonNoiseZ * expected_noise
usePoissonNoiseFiltering = false;
poissonNoiseZ = 1.5;

% If true, terminal display of significant results includes raw and
% normalized counts for each mouse. If false, it only shows group means.
printIndividualSampleCountsInTerminal = false;

% Fold change is calculated using:
%   (mean_5xFAD_analysis + foldChangePseudoCount) /
%   (mean_WT_analysis    + foldChangePseudoCount)
foldChangePseudoCount = 1;


if ~exist(outputFolder, "dir")
    mkdir(outputFolder);
end

rng(1);

%% =============================
% EXAMPLE sample list
% Replace filenames with your own CSV basenames in inputFolder.
% Required columns: gene, region_name (and cell_type if useCellTypes=true).
% Groups must be "WT" or "FAD". Keep equal counts in each group (here 2+2).
% Multiple files listed for one sample are merged before counting.
% =============================

samples = struct();

samples(1).name = "WT1";
samples(1).group = "WT";
samples(1).files = {
    "WT_animal1_sectionA.csv"
    "WT_animal1_sectionB.csv"
    };

samples(2).name = "WT2";
samples(2).group = "WT";
samples(2).files = {
    "WT_animal2.csv"
    };

samples(3).name = "FAD1";
samples(3).group = "FAD";
samples(3).files = {
    "FAD_animal1_sectionA.csv"
    "FAD_animal1_sectionB.csv"
    };

samples(4).name = "FAD2";
samples(4).group = "FAD";
samples(4).files = {
    "FAD_animal2.csv"
    };

sampleNames = string({samples.name});
sampleGroups = string({samples.group});

wtSampleNames = sampleNames(sampleGroups == "WT");
fadSampleNames = sampleNames(sampleGroups == "FAD");

fprintf("Number of WT samples: %d\n", numel(wtSampleNames));
fprintf("Number of 5xFAD samples: %d\n", numel(fadSampleNames));

if numel(wtSampleNames) ~= numel(fadSampleNames)
    error("The number of WT samples and 5xFAD samples must be equal.");
end

if numel(wtSampleNames) < 2 || numel(fadSampleNames) < 2
    error("At least two samples per group are required for statistical testing.");
end

%% =============================
% Diagnostic: gene-name variants before counting
% =============================

geneNameDiagnostics = generateGeneNameDiagnostics(samples, inputFolder, outputFolder);

fprintf("\n=============================\n");
fprintf("Gene-name variant diagnostics, raw variants > 1\n");
fprintf("=============================\n\n");

disp(geneNameDiagnostics(geneNameDiagnostics.number_of_raw_variants > 1, :));

%% =============================
% Read and count expression per mouse
% =============================

sampleTables = cell(numel(samples), 1);

for i = 1:numel(samples)

    fprintf("Reading sample %s...\n", samples(i).name);

    sampleTables{i} = readAndCountMouseCSV( ...
        samples(i).files, ...
        inputFolder, ...
        samples(i).name, ...
        useCellTypes);
end

%% =============================
% Build combined raw count table
% =============================

if useCellTypes
    keyVars = {'gene', 'region_name', 'cell_type'};
else
    keyVars = {'gene', 'region_name'};
end

rawExpressionTable = buildCohortExpressionTable(sampleTables, sampleNames, keyVars);

% Fill missing sample values with zero
for i = 1:numel(sampleNames)
    currentSample = char(sampleNames(i));
    rawExpressionTable.(currentSample)(ismissing(rawExpressionTable.(currentSample))) = 0;
end

rawExpressionTable = sortrows(rawExpressionTable, keyVars);

writetable(rawExpressionTable, ...
    fullfile(outputFolder, "raw_counts_gene_region_per_mouse.xlsx"));

%% =============================
% Diagnostic: number of unique standardized genes per region
% =============================

geneCountByRegion = summarizeUniqueGenesByRegion(rawExpressionTable, useCellTypes);


fprintf("\n=============================\n");
fprintf("Number of unique standardized genes per region\n");
fprintf("=============================\n\n");

disp(geneCountByRegion);

writetable(geneCountByRegion, ...
    fullfile(outputFolder, "diagnostic_unique_genes_per_region.xlsx"));

if any(geneCountByRegion.number_unique_genes > 100)
    warning("Some regions contain more than 100 standardized analyzed genes. Check diagnostics.");
end

%% =============================
% Diagnostic: raw region names after standardization/exclusion
% =============================

fprintf("\n=============================\n");
fprintf("Region names present in rawExpressionTable after standardization/exclusion\n");
fprintf("=============================\n\n");

regionNamesRaw = unique(rawExpressionTable.region_name);
disp(regionNamesRaw);

rawRegionTable = table(regionNamesRaw, 'VariableNames', {'region_name'});

writetable(rawRegionTable, ...
    fullfile(outputFolder, "diagnostic_region_names_raw.xlsx"));

rawRegionCellTable = table();

if useCellTypes
    rawRegionCellTable = unique(rawExpressionTable(:, {'region_name', 'cell_type'}));
    writetable(rawRegionCellTable, ...
        fullfile(outputFolder, "diagnostic_region_celltype_names_raw.xlsx"));
end

%% =============================
% Normalization for analysisExpressionTable
% =============================

rawMatrix = rawExpressionTable{:, cellstr(sampleNames)};

switch string(normalizationMethod)

    case "TMM"

        [analysisMatrix, normalizationTable] = performTMMNormalization( ...
            rawMatrix, sampleNames, Mtrim, Atrim);

        fprintf("\nTMM normalization was applied across all mice.\n");

    case "quantile"

        analysisMatrix = performQuantileNormalization(rawMatrix);

        normalizationTable = table( ...
            sampleNames(:), ...
            sum(rawMatrix, 1)', ...
            repmat("quantile", numel(sampleNames), 1), ...
            'VariableNames', {'sample', ...
                              'library_size_raw', ...
                              'normalization_method'});

        fprintf("\nQuantile normalization was applied across all mice.\n");

    case "none"

        analysisMatrix = rawMatrix;

        normalizationTable = table( ...
            sampleNames(:), ...
            sum(rawMatrix, 1)', ...
            repmat("none", numel(sampleNames), 1), ...
            'VariableNames', {'sample', ...
                              'library_size_raw', ...
                              'normalization_method'});

        fprintf("\nNo normalization was applied. Raw counts are used for analysisExpressionTable.\n");

    otherwise

        error("Unknown normalizationMethod: %s. Use 'TMM', 'quantile', or 'none'.", normalizationMethod);
end

analysisExpressionTable = rawExpressionTable(:, keyVars);

for i = 1:numel(sampleNames)
    analysisExpressionTable.(char(sampleNames(i))) = analysisMatrix(:, i);
end

writetable(analysisExpressionTable, ...
    fullfile(outputFolder, "analysis_expression_gene_region_per_mouse.xlsx"));

writetable(normalizationTable, ...
    fullfile(outputFolder, "normalization_factors_or_summary.xlsx"));

%% =============================
% Compare WT vs 5xFAD within each region
% =============================

[allResultsTable, filterDiagnosticsTable] = compareGroupsByRegion( ...
    analysisExpressionTable, ...
    rawExpressionTable, ...
    sampleNames, ...
    wtSampleNames, ...
    fadSampleNames, ...
    useCellTypes, ...
    useExpressionQuantileFiltering, ...
    expression_quantile_cutoff, ...
    useFoldChangeFiltering, ...
    foldChangeCutoff, ...
    foldChangePseudoCount, ...
    useMinimumMeanExpressionFiltering, ...
    minimumMeanExpression, ...
    useSampleSupportFiltering, ...
    minimumNumberOfExpressingSamples, ...
    sampleSupportExpressionThreshold, ...
    useVarianceFiltering, ...
    minimumVarianceAcrossSamples, ...
    usePoissonNoiseFiltering, ...
    poissonNoiseZ, ...
    statisticalTest, ...
    rnaseqdeVarianceLink, ...
    rnaseqdeFDRMethod, ...
    useExhaustivePermutations, ...
    number_of_realizations, ...
    pValueMode);

if height(allResultsTable) == 0
    error("No gene-region pairs passed the filters. Relax filtering thresholds and rerun.");
end

%% =============================
% Diagnostic: filtering summary
% =============================

fprintf("\n=============================\n");
fprintf("Filter diagnostics by region\n");
fprintf("=============================\n\n");

disp(filterDiagnosticsTable);

writetable(filterDiagnosticsTable, ...
    fullfile(outputFolder, "diagnostic_filter_counts_by_region.xlsx"));

%% =============================
% FDR
% =============================

if statisticalTest == "rnaseqde"

    % rnaseqde is run separately per region, so rnaseqde_AdjustedPValue
    % is region-level FDR. For global reporting, recalculate FDR across
    % all tested gene-region pairs using the pooled raw P values.

    allResultsTable.rnaseqde_region_level_AdjustedPValue = ...
        allResultsTable.rnaseqde_AdjustedPValue;

    if rnaseqdeFDRMethod == "storey"
        allResultsTable.FDR = mafdr(allResultsTable.pvalue);
    elseif rnaseqdeFDRMethod == "bh"
        allResultsTable.FDR = bhFDR(allResultsTable.pvalue);
    else
        error("Unknown rnaseqdeFDRMethod: %s", rnaseqdeFDRMethod);
    end

else

    allResultsTable.FDR = bhFDR(allResultsTable.pvalue);

end

if useCellTypes
    allResultsTable = sortrows(allResultsTable, ...
        {'region_name', 'cell_type', 'FDR'}, ...
        {'ascend', 'ascend', 'ascend'});
else
    allResultsTable = sortrows(allResultsTable, ...
        {'region_name', 'FDR'}, ...
        {'ascend', 'ascend'});
end

%% =============================
% Significant results
% =============================

if applyPCutoffForReporting

    significantLocations = ...
        isfinite(allResultsTable.pvalue) & ...
        allResultsTable.pvalue < cutoff_p;

    reportingMode = sprintf("P-value only: P < %.4g", cutoff_p);

else

    significantLocations = ...
        isfinite(allResultsTable.FDR) & ...
        allResultsTable.FDR < cutoff_FDR;

    reportingMode = sprintf("FDR only: FDR < %.4g", cutoff_FDR);

end

significantResultsTable = allResultsTable(significantLocations, :);

if applyPCutoffForReporting
    significantResultsTable = sortrows(significantResultsTable, 'pvalue', 'ascend');
else
    significantResultsTable = sortrows(significantResultsTable, 'FDR', 'ascend');
end

fprintf("\n=============================\n");
fprintf("Reporting mode\n");
fprintf("=============================\n");
fprintf("%s\n", reportingMode);

fprintf("\nNormalization method for filtering/t-test/permutation: %s\n", normalizationMethod);
fprintf("Statistical test: %s\n", statisticalTest);

if statisticalTest == "rnaseqde"
    fprintf("rnaseqde VarianceLink: %s\n", rnaseqdeVarianceLink);
    fprintf("rnaseqde FDRMethod: %s\n", rnaseqdeFDRMethod);
    fprintf("rnaseqde used raw integer counts for testing.\n");
end

fprintf("\nExpression quantile filtering enabled: %d\n", useExpressionQuantileFiltering);
fprintf("Expression quantile cutoff: %.3f\n", expression_quantile_cutoff);

fprintf("\nFold-change filtering enabled: %d\n", useFoldChangeFiltering);
fprintf("Fold-change cutoff: %.3f\n", foldChangeCutoff);
fprintf("Fold-change pseudo-count: %.3f\n", foldChangePseudoCount);

fprintf("\nMinimum mean expression filtering enabled: %d\n", useMinimumMeanExpressionFiltering);
fprintf("Minimum mean expression: %.3f\n", minimumMeanExpression);

fprintf("\nSample support filtering enabled: %d\n", useSampleSupportFiltering);
fprintf("Minimum number of expressing samples: %d\n", minimumNumberOfExpressingSamples);
fprintf("Sample support expression threshold: %.3f\n", sampleSupportExpressionThreshold);

fprintf("\nVariance filtering enabled: %d\n", useVarianceFiltering);
fprintf("Minimum variance across samples: %.3f\n", minimumVarianceAcrossSamples);

fprintf("\nPoisson-noise filtering enabled: %d\n", usePoissonNoiseFiltering);
fprintf("Poisson noise Z cutoff: %.3f\n", poissonNoiseZ);

fprintf("\nNumber of tested gene-region pairs: %d\n", height(allResultsTable));
fprintf("Number of finite P values: %d\n", sum(isfinite(allResultsTable.pvalue)));
fprintf("Number of P < cutoff_p: %d\n", sum(isfinite(allResultsTable.pvalue) & allResultsTable.pvalue < cutoff_p));
fprintf("Number of FDR < cutoff_FDR: %d\n", sum(isfinite(allResultsTable.FDR) & allResultsTable.FDR < cutoff_FDR));
fprintf("Number of reported significant results: %d\n", height(significantResultsTable));

fprintf("\n=============================\n");
fprintf("Significant results, compact view\n");
fprintf("=============================\n\n");

significantResultsDisplayTable = makeCompactSignificantResultsTable( ...
    significantResultsTable, ...
    sampleNames, ...
    useCellTypes, ...
    printIndividualSampleCountsInTerminal);

disp(significantResultsDisplayTable);

%% =============================
% Statistics on significant genes per region
% =============================

significantGenesPerRegion = summarizeSignificantGenesPerRegion( ...
    significantResultsTable, ...
    useCellTypes);

fprintf("\n=============================\n");
fprintf("Significant genes per region\n");
fprintf("=============================\n\n");

disp(significantGenesPerRegion);

writetable(significantGenesPerRegion, ...
    fullfile(outputFolder, "significant_genes_per_region_summary.xlsx"));

writetable(significantResultsDisplayTable, ...
    fullfile(outputFolder, "significant_WT_vs_5xFAD_region_results_compact.xlsx"));

writetable(allResultsTable, ...
    fullfile(outputFolder, "all_WT_vs_5xFAD_region_results.xlsx"));

writetable(significantResultsTable, ...
    fullfile(outputFolder, "significant_WT_vs_5xFAD_region_results.xlsx"));

%% =============================
% Count tests per region
% =============================

if useCellTypes
    testCounts = groupsummary(allResultsTable, ...
        ["region_name", "cell_type"], ...
        @(x) sum(isfinite(x)), ...
        "pvalue");
else
    testCounts = groupsummary(allResultsTable, ...
        "region_name", ...
        @(x) sum(isfinite(x)), ...
        "pvalue");
end

testCounts.Properties.VariableNames(end) = {'number_of_tests'};
testCounts = removevars(testCounts, "GroupCount");

%% =============================
% Save workspace
% =============================

save(fullfile(outputFolder, "WT_vs_5xFAD_region_expression_workspace.mat"), ...
    "rawExpressionTable", ...
    "analysisExpressionTable", ...
    "normalizationTable", ...
    "allResultsTable", ...
    "significantResultsTable", ...
    "filterDiagnosticsTable", ...
    "testCounts", ...
    "rawRegionTable", ...
    "rawRegionCellTable", ...
    "geneNameDiagnostics", ...
    "geneCountByRegion", ...
    "samples", ...
    "sampleNames", ...
    "sampleGroups", ...
    "wtSampleNames", ...
    "fadSampleNames", ...
    "useCellTypes", ...
    "normalizationMethod", ...
    "Mtrim", ...
    "Atrim", ...
    "statisticalTest", ...
    "rnaseqdeVarianceLink", ...
    "rnaseqdeFDRMethod", ...
    "useExpressionQuantileFiltering", ...
    "expression_quantile_cutoff", ...
    "useFoldChangeFiltering", ...
    "foldChangeCutoff", ...
    "foldChangePseudoCount", ...
    "useMinimumMeanExpressionFiltering", ...
    "minimumMeanExpression", ...
    "useSampleSupportFiltering", ...
    "minimumNumberOfExpressingSamples", ...
    "sampleSupportExpressionThreshold", ...
    "useVarianceFiltering", ...
    "minimumVarianceAcrossSamples", ...
    "usePoissonNoiseFiltering", ...
    "poissonNoiseZ", ...
    "useExhaustivePermutations", ...
    "number_of_realizations", ...
    "pValueMode", ...
    "cutoff_p", ...
    "cutoff_FDR", ...
    "applyPCutoffForReporting", ...
    "significantResultsDisplayTable", ...
    "significantGenesPerRegion", ...
    "printIndividualSampleCountsInTerminal", ...
    "reportingMode");

fprintf("\nAnalysis complete. Results saved to folder:\n%s\n", outputFolder);

%% ========================================================================
% Local functions
% ========================================================================

function sampleTable = readAndCountMouseCSV(fileList, inputFolder, sampleName, useCellTypes)

    if useCellTypes
        requiredCols = {'gene', 'region_name', 'cell_type'};
    else
        requiredCols = {'gene', 'region_name'};
    end

    allRows = table();

    for f = 1:numel(fileList)

        currentFile = string(fileList{f});
        currentPath = fullfile(inputFolder, currentFile);

        if ~isfile(currentPath)
            error("File not found: %s", currentPath);
        end

        T = readtable(currentPath, "TextType", "string");

        missingCols = setdiff(string(requiredCols), string(T.Properties.VariableNames));

        if ~isempty(missingCols)
            error("File %s is missing required column(s): %s", ...
                currentFile, strjoin(missingCols, ", "));
        end

        T = T(:, requiredCols);
        T = cleanInputTable(T, useCellTypes);

        allRows = [allRows; T];
    end

    countTable = makeCountTable(allRows, useCellTypes);

    countTable = renameTableVariableRobust(countTable, ...
        "copy_count", ...
        char(sampleName));

    sampleTable = countTable;
end


function T = cleanInputTable(T, useCellTypes)

    %% Preserve original gene names for diagnostics
    T.gene_original = string(T.gene);

    %% Standardize gene names before counting
    T.gene = standardizeGeneSymbol(T.gene);

    %% Merge known gene aliases
    T.gene(T.gene == "IRG1")  = "ACOD1";
    T.gene(T.gene == "INOS")  = "NOS1";
    T.gene(T.gene == "CCRLG") = "CCRL2";
    T.gene(T.gene == "IDB3")  = "ID3";

    % Handle HIF-1 variants.
    % Target name is HIF-1Α with Greek uppercase alpha.
    T.gene(T.gene == "HIF-1?") = "HIF-1Α";
    T.gene(T.gene == "HIF-1A") = "HIF-1Α";
    T.gene(T.gene == "HIF1A")  = "HIF-1Α";

    %% Exclude genes that should not be analyzed
    excludeGenes = ["NEG1", "NEG2", "NEG3", "FCGR3"];

    excludeGeneRows = ismember(T.gene, excludeGenes);

    if any(excludeGeneRows)
        fprintf("Excluding %d rows from genes: NEG1 / NEG2 / NEG3 / FCGR3\n", ...
            sum(excludeGeneRows));
    end

    T = T(~excludeGeneRows, :);

    %% Clean region names
    T.region_name = strtrim(string(T.region_name));

    %% Remove rows with missing gene or missing region
    keepRows = ~(ismissing(T.gene) | strlength(T.gene) == 0 | ...
                 ismissing(T.region_name) | strlength(T.region_name) == 0);

    T = T(keepRows, :);

    %% Standardize hippocampal region names
    % Rules:
    %   SM        -> DG-CA1
    %   inner_DG  -> Hilus
    %   upper_CA1 -> SO
    %   unassigned, under_DG, and CA2 are excluded from the analysis

    regionLower = lower(strtrim(T.region_name));

    T.region_name(regionLower == "sm") = "DG-CA1";
    T.region_name(regionLower == "inner_dg") = "Hilus";
    T.region_name(regionLower == "upper_ca1") = "SO";

    %% Exclude regions that should not be tested
    regionLower = lower(strtrim(T.region_name));

    excludeRegionRows = ...
        regionLower == "unassigned" | ...
        regionLower == "under_dg" | ...
        regionLower == "ca2";

    if any(excludeRegionRows)
        fprintf("Excluding %d rows from regions: unassigned / under_DG / CA2\n", ...
            sum(excludeRegionRows));
    end

    T = T(~excludeRegionRows, :);

    %% Clean cell type if needed
    if useCellTypes
        T.cell_type = string(T.cell_type);
        T.cell_type(ismissing(T.cell_type) | strlength(strtrim(T.cell_type)) == 0) = "unknown";
        T.cell_type = strtrim(T.cell_type);
    end
end


function geneClean = standardizeGeneSymbol(geneRaw)

    geneClean = string(geneRaw);

    % Replace non-breaking spaces with regular spaces
    geneClean = replace(geneClean, char(160), " ");

    % Remove control characters
    geneClean = regexprep(geneClean, '[\x00-\x1F\x7F]', '');

    % Remove leading/trailing whitespace
    geneClean = strtrim(geneClean);

    % Remove surrounding single or double quotes
    geneClean = regexprep(geneClean, '^["'']+|["'']+$', '');

    % Remove internal whitespace
    geneClean = regexprep(geneClean, '\s+', '');

    % Standardize common alpha representations
    geneClean = replace(geneClean, "α", "Α");  % Greek lowercase alpha -> Greek uppercase alpha

    % Make grouping case-insensitive
    geneClean = upper(geneClean);
end


function summaryTable = makeCountTable(T, useCellTypes)

    if height(T) == 0

        if useCellTypes
            summaryTable = table( ...
                strings(0,1), strings(0,1), strings(0,1), zeros(0,1), ...
                'VariableNames', {'gene', 'region_name', 'cell_type', 'copy_count'});
        else
            summaryTable = table( ...
                strings(0,1), strings(0,1), zeros(0,1), ...
                'VariableNames', {'gene', 'region_name', 'copy_count'});
        end

        return;
    end

    if useCellTypes

        [G, geneNames, regionNames, cellTypes] = findgroups( ...
            T.gene, ...
            T.region_name, ...
            T.cell_type);

        copyCounts = splitapply(@numel, T.gene, G);

        summaryTable = table( ...
            geneNames, ...
            regionNames, ...
            cellTypes, ...
            copyCounts, ...
            'VariableNames', {'gene', 'region_name', 'cell_type', 'copy_count'});

    else

        [G, geneNames, regionNames] = findgroups( ...
            T.gene, ...
            T.region_name);

        copyCounts = splitapply(@numel, T.gene, G);

        summaryTable = table( ...
            geneNames, ...
            regionNames, ...
            copyCounts, ...
            'VariableNames', {'gene', 'region_name', 'copy_count'});
    end
end


function cohortTable = buildCohortExpressionTable(sampleTables, sampleNames, keyVars)

    cohortTable = sampleTables{1};

    for i = 2:numel(sampleTables)

        cohortTable = outerjoin( ...
            cohortTable, ...
            sampleTables{i}, ...
            "Keys", keyVars, ...
            "MergeKeys", true);
    end

    for i = 1:numel(sampleNames)

        sampleName = char(sampleNames(i));

        if ~ismember(sampleName, cohortTable.Properties.VariableNames)
            cohortTable.(sampleName) = zeros(height(cohortTable), 1);
        end
    end
end


function T = renameTableVariableRobust(T, oldName, newName)

    oldName = char(oldName);
    newName = char(newName);

    idx = strcmp(T.Properties.VariableNames, oldName);

    if ~any(idx)
        error("Variable '%s' was not found in the table.", oldName);
    end

    T.Properties.VariableNames(idx) = {newName};
end


function geneNameDiagnostics = generateGeneNameDiagnostics(samples, inputFolder, outputFolder)

    %% Compact diagnostic, does not write all transcript-level rows to Excel

    rawVariantCountsByFile = table();

    for i = 1:numel(samples)

        sampleName = string(samples(i).name);

        for f = 1:numel(samples(i).files)

            currentFile = string(samples(i).files{f});
            currentPath = fullfile(inputFolder, currentFile);

            if ~isfile(currentPath)
                error("File not found: %s", currentPath);
            end

            T = readtable(currentPath, "TextType", "string");

            if ~ismember("gene", string(T.Properties.VariableNames))
                error("File %s is missing required column: gene", currentFile);
            end

            rawGene = string(T.gene);
            standardizedGene = standardizeGeneSymbol(rawGene);

            keepRows = ~(ismissing(standardizedGene) | strlength(standardizedGene) == 0);
            rawGene = rawGene(keepRows);
            standardizedGene = standardizedGene(keepRows);

            [G, rawGeneUnique, standardizedGeneUnique] = findgroups(rawGene, standardizedGene);
            variantCounts = splitapply(@numel, rawGene, G);

            tempTable = table( ...
                repmat(sampleName, numel(rawGeneUnique), 1), ...
                repmat(currentFile, numel(rawGeneUnique), 1), ...
                rawGeneUnique, ...
                standardizedGeneUnique, ...
                variantCounts, ...
                'VariableNames', {'sample', ...
                                  'file', ...
                                  'raw_gene', ...
                                  'standardized_gene', ...
                                  'number_of_occurrences'});

            rawVariantCountsByFile = [rawVariantCountsByFile; tempTable];
        end
    end

    [Gmap, rawGeneMap, standardizedGeneMap] = findgroups( ...
        rawVariantCountsByFile.raw_gene, ...
        rawVariantCountsByFile.standardized_gene);

    totalOccurrences = splitapply( ...
        @sum, ...
        rawVariantCountsByFile.number_of_occurrences, ...
        Gmap);

    rawToStandardizedMap = table( ...
        rawGeneMap, ...
        standardizedGeneMap, ...
        totalOccurrences, ...
        'VariableNames', {'raw_gene', ...
                          'standardized_gene', ...
                          'total_occurrences'});

    rawToStandardizedMap = sortrows(rawToStandardizedMap, ...
        {'standardized_gene', 'raw_gene'});

    standardizedGenes = unique(rawToStandardizedMap.standardized_gene);

    geneNameDiagnostics = table();

    for i = 1:numel(standardizedGenes)

        currentGene = standardizedGenes(i);

        rows = rawToStandardizedMap.standardized_gene == currentGene;

        rawVariants = rawToStandardizedMap.raw_gene(rows);
        rawVariantsJoined = strjoin(rawVariants, " | ");

        totalCountsForGene = sum(rawToStandardizedMap.total_occurrences(rows));

        newRow = table( ...
            currentGene, ...
            numel(rawVariants), ...
            totalCountsForGene, ...
            rawVariantsJoined, ...
            'VariableNames', {'standardized_gene', ...
                              'number_of_raw_variants', ...
                              'total_occurrences', ...
                              'raw_gene_variants'});

        geneNameDiagnostics = [geneNameDiagnostics; newRow];
    end

    geneNameDiagnostics = sortrows(geneNameDiagnostics, ...
        {'number_of_raw_variants', 'total_occurrences'}, ...
        {'descend', 'descend'});

    writetable(geneNameDiagnostics, ...
        fullfile(outputFolder, "diagnostic_gene_name_variants.xlsx"));

    writetable(rawToStandardizedMap, ...
        fullfile(outputFolder, "diagnostic_raw_to_standardized_gene_map.xlsx"));

    writetable(rawVariantCountsByFile, ...
        fullfile(outputFolder, "diagnostic_raw_gene_variant_counts_by_file.xlsx"));

    save(fullfile(outputFolder, "diagnostic_gene_name_variants.mat"), ...
        "geneNameDiagnostics", ...
        "rawToStandardizedMap", ...
        "rawVariantCountsByFile");
end


function geneCountByRegion = summarizeUniqueGenesByRegion(rawExpressionTable, useCellTypes)

    geneCountByRegion = table();

    if useCellTypes
        [G, regionList, cellTypeList] = findgroups( ...
            rawExpressionTable.region_name, ...
            rawExpressionTable.cell_type);
    else
        [G, regionList] = findgroups(rawExpressionTable.region_name);
    end

    for i = 1:max(G)

        rows = G == i;

        numberUniqueGenes = numel(unique(rawExpressionTable.gene(rows)));

        if useCellTypes
            newRow = table( ...
                regionList(i), ...
                cellTypeList(i), ...
                numberUniqueGenes, ...
                'VariableNames', {'region_name', ...
                                  'cell_type', ...
                                  'number_unique_genes'});
        else
            newRow = table( ...
                regionList(i), ...
                numberUniqueGenes, ...
                'VariableNames', {'region_name', ...
                                  'number_unique_genes'});
        end

        geneCountByRegion = [geneCountByRegion; newRow];
    end
end


function [normalizedMatrix, normalizationTable] = performTMMNormalization(rawMatrix, sampleNames, Mtrim, Atrim)

    rawMatrix = double(rawMatrix);

    nSamples = size(rawMatrix, 2);

    librarySizes = sum(rawMatrix, 1);

    if any(librarySizes == 0)
        error("At least one sample has total library size zero. TMM normalization cannot be performed.");
    end

    % Choose reference sample as the one with library size closest to median.
    medianLibrarySize = median(librarySizes);
    [~, refIdx] = min(abs(librarySizes - medianLibrarySize));

    rawFactors = ones(1, nSamples);

    for s = 1:nSamples

        if s == refIdx
            rawFactors(s) = 1;
        else
            rawFactors(s) = calculateTMMFactor( ...
                rawMatrix(:, refIdx), ...
                rawMatrix(:, s), ...
                librarySizes(refIdx), ...
                librarySizes(s), ...
                Mtrim, ...
                Atrim);
        end
    end

    % Center factors so geometric mean is 1.
    positiveFactors = rawFactors(rawFactors > 0 & isfinite(rawFactors));

    if isempty(positiveFactors)
        warning("No valid TMM factors. Using factors of 1.");
        tmmFactors = ones(size(rawFactors));
    else
        geometricMeanFactor = exp(mean(log(positiveFactors)));
        tmmFactors = rawFactors ./ geometricMeanFactor;
    end

    tmmFactors(~isfinite(tmmFactors) | tmmFactors <= 0) = 1;

    effectiveLibrarySizes = librarySizes .* tmmFactors;

    meanEffectiveLibrarySize = mean(effectiveLibrarySizes);

    normalizedMatrix = rawMatrix ./ effectiveLibrarySizes .* meanEffectiveLibrarySize;

    normalizationTable = table( ...
        sampleNames(:), ...
        librarySizes(:), ...
        rawFactors(:), ...
        tmmFactors(:), ...
        effectiveLibrarySizes(:), ...
        repmat(sampleNames(refIdx), nSamples, 1), ...
        repmat("TMM", nSamples, 1), ...
        'VariableNames', {'sample', ...
                          'library_size_raw', ...
                          'TMM_factor_raw_vs_reference', ...
                          'TMM_factor_centered', ...
                          'effective_library_size', ...
                          'TMM_reference_sample', ...
                          'normalization_method'});
end


function tmmFactor = calculateTMMFactor(xRef, yTest, libRef, libTest, Mtrim, Atrim)

    xRef = double(xRef);
    yTest = double(yTest);

    valid = xRef > 0 & yTest > 0;

    if sum(valid) < 5
        tmmFactor = 1;
        return;
    end

    x = xRef(valid);
    y = yTest(valid);

    M = log2((y ./ libTest) ./ (x ./ libRef));
    A = 0.5 .* log2((y ./ libTest) .* (x ./ libRef));

    V = ((libTest - y) ./ (libTest .* y)) + ...
        ((libRef  - x) ./ (libRef  .* x));

    weights = 1 ./ V;

    validMA = isfinite(M) & isfinite(A) & isfinite(weights) & weights > 0;

    M = M(validMA);
    A = A(validMA);
    weights = weights(validMA);

    if numel(M) < 5
        tmmFactor = 1;
        return;
    end

    lowerM = quantile(M, Mtrim / 2);
    upperM = quantile(M, 1 - Mtrim / 2);

    keepM = M >= lowerM & M <= upperM;

    lowerA = quantile(A, Atrim / 2);
    upperA = quantile(A, 1 - Atrim / 2);

    keepA = A >= lowerA & A <= upperA;

    keep = keepM & keepA;

    if sum(keep) < 3
        tmmFactor = 1;
        return;
    end

    Mtrimmed = M(keep);
    weightsTrimmed = weights(keep);

    weightedMeanM = sum(weightsTrimmed .* Mtrimmed) ./ sum(weightsTrimmed);

    tmmFactor = 2 .^ weightedMeanM;

    if ~isfinite(tmmFactor) || tmmFactor <= 0
        tmmFactor = 1;
    end
end


function normalizedMatrix = performQuantileNormalization(rawMatrix)

    rawMatrix = double(rawMatrix);

    if exist("quantilenorm", "file") == 2
        normalizedMatrix = quantilenorm(rawMatrix);
    else
        warning("quantilenorm was not found. Using local quantile normalization function.");
        normalizedMatrix = localQuantileNormalize(rawMatrix);
    end
end


function normalizedMatrix = localQuantileNormalize(rawMatrix)

    rawMatrix = double(rawMatrix);

    [sortedValues, sortIdx] = sort(rawMatrix, 1, "ascend");

    meanSortedValues = mean(sortedValues, 2, "omitnan");

    normalizedMatrix = zeros(size(rawMatrix));

    for col = 1:size(rawMatrix, 2)
        normalizedMatrix(sortIdx(:, col), col) = meanSortedValues;
    end
end


function [allResultsTable, filterDiagnosticsTable] = compareGroupsByRegion( ...
    analysisExpressionTable, ...
    rawExpressionTable, ...
    sampleNames, ...
    wtSampleNames, ...
    fadSampleNames, ...
    useCellTypes, ...
    useExpressionQuantileFiltering, ...
    expression_quantile_cutoff, ...
    useFoldChangeFiltering, ...
    foldChangeCutoff, ...
    foldChangePseudoCount, ...
    useMinimumMeanExpressionFiltering, ...
    minimumMeanExpression, ...
    useSampleSupportFiltering, ...
    minimumNumberOfExpressingSamples, ...
    sampleSupportExpressionThreshold, ...
    useVarianceFiltering, ...
    minimumVarianceAcrossSamples, ...
    usePoissonNoiseFiltering, ...
    poissonNoiseZ, ...
    statisticalTest, ...
    rnaseqdeVarianceLink, ...
    rnaseqdeFDRMethod, ...
    useExhaustivePermutations, ...
    number_of_realizations, ...
    pValueMode)

    allResultsTable = table();
    filterDiagnosticsTable = table();

    if useCellTypes
        [G, regionList, cellTypeList] = findgroups( ...
            analysisExpressionTable.region_name, ...
            analysisExpressionTable.cell_type);
    else
        [G, regionList] = findgroups(analysisExpressionTable.region_name);
    end

    numberOfGroups = max(G);

    for groupIndex = 1:numberOfGroups

        groupRows = G == groupIndex;

        regionName = regionList(groupIndex);

        if useCellTypes
            cellTypeName = cellTypeList(groupIndex);
        end

        groupAnalysisTable = analysisExpressionTable(groupRows, :);
        groupRawTable = rawExpressionTable(groupRows, :);

        analysisMatrix = groupAnalysisTable{:, cellstr(sampleNames)};
        rawMatrix = groupRawTable{:, cellstr(sampleNames)};

        wtAnalysis = groupAnalysisTable{:, cellstr(wtSampleNames)};
        fadAnalysis = groupAnalysisTable{:, cellstr(fadSampleNames)};

        meanWT = mean(wtAnalysis, 2);
        meanFAD = mean(fadAnalysis, 2);

        diffFADminusWT = meanFAD - meanWT;

        foldChangeRatio_plusPseudo = ...
            (meanFAD + foldChangePseudoCount) ./ ...
            (meanWT  + foldChangePseudoCount);

        absFoldChange_plusPseudo = max( ...
            foldChangeRatio_plusPseudo, ...
            1 ./ foldChangeRatio_plusPseudo);

        log2FC_plusPseudo = log2(foldChangeRatio_plusPseudo);

        %% Additional expression statistics before filtering

        meanAllSamples = mean(analysisMatrix, 2);
        varianceAcrossSamples = var(analysisMatrix, 0, 2);

        numberExpressingSamples = ...
            sum(analysisMatrix >= sampleSupportExpressionThreshold, 2);

        absoluteDifference = abs(meanFAD - meanWT);

        poissonNoiseEstimate = sqrt(meanWT + meanFAD);
        poissonNoiseZObserved = absoluteDifference ./ poissonNoiseEstimate;
        poissonNoiseZObserved(~isfinite(poissonNoiseZObserved)) = NaN;

        %% Expression quantile filter

        wtThreshold = quantile(meanWT, expression_quantile_cutoff);
        fadThreshold = quantile(meanFAD, expression_quantile_cutoff);

        expressionFilterCandidate = ...
            (meanWT >= wtThreshold) | ...
            (meanFAD >= fadThreshold);

        if useExpressionQuantileFiltering
            expressionFilterPass = expressionFilterCandidate;
        else
            expressionFilterPass = true(height(groupAnalysisTable), 1);
        end

        %% Fold-change filter

        foldChangeFilterPass = true(height(groupAnalysisTable), 1);

        if useFoldChangeFiltering
            foldChangeFilterPass = absFoldChange_plusPseudo >= foldChangeCutoff;
        end

        %% Minimum mean expression filter

        minimumMeanExpressionFilterPass = true(height(groupAnalysisTable), 1);

        if useMinimumMeanExpressionFiltering
            minimumMeanExpressionFilterPass = meanAllSamples >= minimumMeanExpression;
        end

        %% Sample support filter

        sampleSupportFilterPass = true(height(groupAnalysisTable), 1);

        if useSampleSupportFiltering
            sampleSupportFilterPass = ...
                numberExpressingSamples >= minimumNumberOfExpressingSamples;
        end

        %% Variance filter

        varianceFilterPass = true(height(groupAnalysisTable), 1);

        if useVarianceFiltering
            varianceFilterPass = ...
                varianceAcrossSamples >= minimumVarianceAcrossSamples;
        end

        %% Poisson-noise-aware effect-size filter

        poissonNoiseFilterPass = true(height(groupAnalysisTable), 1);

        if usePoissonNoiseFiltering
            poissonNoiseFilterPass = ...
                poissonNoiseEstimate > 0 & ...
                absoluteDifference >= poissonNoiseZ .* poissonNoiseEstimate;
        end

        %% Combined filter

        keepGenes = ...
            expressionFilterPass & ...
            foldChangeFilterPass & ...
            minimumMeanExpressionFilterPass & ...
            sampleSupportFilterPass & ...
            varianceFilterPass & ...
            poissonNoiseFilterPass;

        %% Diagnostic row before subsetting

        totalBeforeFilter = height(groupAnalysisTable);

        numberWouldPassExpressionQuantileFilter = sum(expressionFilterCandidate);
        numberPassExpressionQuantile = sum(expressionFilterPass);
        numberPassFoldChange = sum(foldChangeFilterPass);
        numberPassMinimumMeanExpression = sum(minimumMeanExpressionFilterPass);
        numberPassSampleSupport = sum(sampleSupportFilterPass);
        numberPassVariance = sum(varianceFilterPass);
        numberPassPoissonNoise = sum(poissonNoiseFilterPass);
        numberPassAllFilters = sum(keepGenes);

        if useCellTypes
            diagnosticRow = table( ...
                regionName, ...
                cellTypeName, ...
                totalBeforeFilter, ...
                numberWouldPassExpressionQuantileFilter, ...
                numberPassExpressionQuantile, ...
                numberPassFoldChange, ...
                numberPassMinimumMeanExpression, ...
                numberPassSampleSupport, ...
                numberPassVariance, ...
                numberPassPoissonNoise, ...
                numberPassAllFilters, ...
                wtThreshold, ...
                fadThreshold, ...
                'VariableNames', {'region_name', ...
                                  'cell_type', ...
                                  'number_before_filter', ...
                                  'number_would_pass_expression_quantile_filter', ...
                                  'number_pass_expression_quantile_filter', ...
                                  'number_pass_fold_change_filter', ...
                                  'number_pass_minimum_mean_expression_filter', ...
                                  'number_pass_sample_support_filter', ...
                                  'number_pass_variance_filter', ...
                                  'number_pass_poisson_noise_filter', ...
                                  'number_pass_all_filters', ...
                                  'WT_expression_quantile_threshold', ...
                                  'FAD_expression_quantile_threshold'});
        else
            diagnosticRow = table( ...
                regionName, ...
                totalBeforeFilter, ...
                numberWouldPassExpressionQuantileFilter, ...
                numberPassExpressionQuantile, ...
                numberPassFoldChange, ...
                numberPassMinimumMeanExpression, ...
                numberPassSampleSupport, ...
                numberPassVariance, ...
                numberPassPoissonNoise, ...
                numberPassAllFilters, ...
                wtThreshold, ...
                fadThreshold, ...
                'VariableNames', {'region_name', ...
                                  'number_before_filter', ...
                                  'number_would_pass_expression_quantile_filter', ...
                                  'number_pass_expression_quantile_filter', ...
                                  'number_pass_fold_change_filter', ...
                                  'number_pass_minimum_mean_expression_filter', ...
                                  'number_pass_sample_support_filter', ...
                                  'number_pass_variance_filter', ...
                                  'number_pass_poisson_noise_filter', ...
                                  'number_pass_all_filters', ...
                                  'WT_expression_quantile_threshold', ...
                                  'FAD_expression_quantile_threshold'});
        end

        filterDiagnosticsTable = [filterDiagnosticsTable; diagnosticRow];

        if ~any(keepGenes)
            continue;
        end

        %% Subset to genes that pass filters

        groupAnalysisTable = groupAnalysisTable(keepGenes, :);
        groupRawTable = groupRawTable(keepGenes, :);

        analysisMatrix = analysisMatrix(keepGenes, :);
        rawMatrix = rawMatrix(keepGenes, :);

        meanWT = meanWT(keepGenes);
        meanFAD = meanFAD(keepGenes);
        diffFADminusWT = diffFADminusWT(keepGenes);

        foldChangeRatio_plusPseudo = foldChangeRatio_plusPseudo(keepGenes);
        absFoldChange_plusPseudo = absFoldChange_plusPseudo(keepGenes);
        log2FC_plusPseudo = log2FC_plusPseudo(keepGenes);

        meanAllSamples = meanAllSamples(keepGenes);
        varianceAcrossSamples = varianceAcrossSamples(keepGenes);
        numberExpressingSamples = numberExpressingSamples(keepGenes);
        absoluteDifference = absoluteDifference(keepGenes);
        poissonNoiseEstimate = poissonNoiseEstimate(keepGenes);
        poissonNoiseZObserved = poissonNoiseZObserved(keepGenes);

        expressionFilterCandidate = expressionFilterCandidate(keepGenes);
        expressionFilterPass = expressionFilterPass(keepGenes);
        foldChangeFilterPass = foldChangeFilterPass(keepGenes);
        minimumMeanExpressionFilterPass = minimumMeanExpressionFilterPass(keepGenes);
        sampleSupportFilterPass = sampleSupportFilterPass(keepGenes);
        varianceFilterPass = varianceFilterPass(keepGenes);
        poissonNoiseFilterPass = poissonNoiseFilterPass(keepGenes);

        wtRaw = groupRawTable{:, cellstr(wtSampleNames)};
        fadRaw = groupRawTable{:, cellstr(fadSampleNames)};

        meanRawWT = mean(wtRaw, 2);
        meanRawFAD = mean(fadRaw, 2);

        wtAnalysisFiltered = groupAnalysisTable{:, cellstr(wtSampleNames)};
        fadAnalysisFiltered = groupAnalysisTable{:, cellstr(fadSampleNames)};

        %% Statistical testing

        switch string(statisticalTest)

            case "rnaseqde"

                rawCountTableForRnaseqde = buildRnaseqdeCountTable( ...
                    groupRawTable, ...
                    sampleNames);

                rnaseqdeOutput = rnaseqde( ...
                    rawCountTableForRnaseqde, ...
                    wtSampleNames, ...
                    fadSampleNames, ...
                    VarianceLink = rnaseqdeVarianceLink, ...
                    FDRMethod = rnaseqdeFDRMethod, ...
                    IDColumns = "gene");

                pValues = rnaseqdeOutput.PValue;
                tStats = nan(size(pValues));
                degreesOfFreedom = nan(size(pValues));
                numberOfPermutationsUsed = NaN;

                rnaseqdeMean1 = rnaseqdeOutput.Mean1;
                rnaseqdeMean2 = rnaseqdeOutput.Mean2;
                rnaseqdeLog2FoldChange = rnaseqdeOutput.Log2FoldChange;
                rnaseqdeAdjustedPValue = rnaseqdeOutput.AdjustedPValue;

            case "welch_ttest_log2_plus1"

                wtForTest = log2(wtAnalysisFiltered + 1);
                fadForTest = log2(fadAnalysisFiltered + 1);

                [pValues, tStats, degreesOfFreedom] = ...
                    welchTTestRows(wtForTest, fadForTest);

                numberOfPermutationsUsed = NaN;

                rnaseqdeMean1 = nan(size(pValues));
                rnaseqdeMean2 = nan(size(pValues));
                rnaseqdeLog2FoldChange = nan(size(pValues));
                rnaseqdeAdjustedPValue = nan(size(pValues));

            case "welch_ttest_raw"

                wtForTest = wtAnalysisFiltered;
                fadForTest = fadAnalysisFiltered;

                [pValues, tStats, degreesOfFreedom] = ...
                    welchTTestRows(wtForTest, fadForTest);

                numberOfPermutationsUsed = NaN;

                rnaseqdeMean1 = nan(size(pValues));
                rnaseqdeMean2 = nan(size(pValues));
                rnaseqdeLog2FoldChange = nan(size(pValues));
                rnaseqdeAdjustedPValue = nan(size(pValues));

            case "permutation"

                permDiffs = calculatePermutationDiffs( ...
                    analysisMatrix, ...
                    numel(wtSampleNames), ...
                    useExhaustivePermutations, ...
                    number_of_realizations);

                pValues = calculatePermutationPValues( ...
                    diffFADminusWT, ...
                    permDiffs, ...
                    pValueMode);

                tStats = nan(size(pValues));
                degreesOfFreedom = nan(size(pValues));

                if useExhaustivePermutations
                    numberOfPermutationsUsed = size(permDiffs, 2);
                else
                    numberOfPermutationsUsed = number_of_realizations;
                end

                rnaseqdeMean1 = nan(size(pValues));
                rnaseqdeMean2 = nan(size(pValues));
                rnaseqdeLog2FoldChange = nan(size(pValues));
                rnaseqdeAdjustedPValue = nan(size(pValues));

            otherwise

                error("Unknown statisticalTest: %s", statisticalTest);
        end

        %% Build result table

        resultTable = groupAnalysisTable(:, {'gene', 'region_name'});

        if useCellTypes
            resultTable.cell_type = groupAnalysisTable.cell_type;
        end

        for s = 1:numel(sampleNames)
            sampleName = char(sampleNames(s));
            resultTable.([sampleName '_raw']) = rawMatrix(:, s);
        end

        for s = 1:numel(sampleNames)
            sampleName = char(sampleNames(s));
            resultTable.([sampleName '_analysis']) = analysisMatrix(:, s);
        end

        resultTable.mean_WT_raw = meanRawWT;
        resultTable.mean_5xFAD_raw = meanRawFAD;

        resultTable.mean_WT_analysis = meanWT;
        resultTable.mean_5xFAD_analysis = meanFAD;

        resultTable.diff_5xFAD_minus_WT = diffFADminusWT;

        resultTable.foldChange_5xFAD_over_WT_plusPseudo = foldChangeRatio_plusPseudo;
        resultTable.absFoldChange_plusPseudo = absFoldChange_plusPseudo;
        resultTable.log2FC_5xFAD_over_WT_plusPseudo = log2FC_plusPseudo;
        resultTable.foldChangePseudoCount = repmat(foldChangePseudoCount, height(resultTable), 1);

        resultTable.mean_expression_all_samples = meanAllSamples;
        resultTable.variance_across_samples = varianceAcrossSamples;
        resultTable.number_expressing_samples = numberExpressingSamples;

        resultTable.absolute_difference_5xFAD_vs_WT = absoluteDifference;
        resultTable.poisson_noise_estimate = poissonNoiseEstimate;
        resultTable.poisson_noise_z_observed = poissonNoiseZObserved;

        resultTable.expression_filter_candidate = expressionFilterCandidate;
        resultTable.expression_filter_pass = expressionFilterPass;
        resultTable.fold_change_filter_pass = foldChangeFilterPass;
        resultTable.minimum_mean_expression_filter_pass = minimumMeanExpressionFilterPass;
        resultTable.sample_support_filter_pass = sampleSupportFilterPass;
        resultTable.variance_filter_pass = varianceFilterPass;
        resultTable.poisson_noise_filter_pass = poissonNoiseFilterPass;

        resultTable.expression_quantile_filter_enabled = ...
            repmat(useExpressionQuantileFiltering, height(resultTable), 1);
        resultTable.expression_quantile_cutoff = ...
            repmat(expression_quantile_cutoff, height(resultTable), 1);
        resultTable.WT_expression_quantile_threshold = ...
            repmat(wtThreshold, height(resultTable), 1);
        resultTable.FAD_expression_quantile_threshold = ...
            repmat(fadThreshold, height(resultTable), 1);

        resultTable.fold_change_filter_enabled = ...
            repmat(useFoldChangeFiltering, height(resultTable), 1);
        resultTable.foldChangeCutoff = ...
            repmat(foldChangeCutoff, height(resultTable), 1);

        resultTable.minimum_mean_expression_filter_enabled = ...
            repmat(useMinimumMeanExpressionFiltering, height(resultTable), 1);
        resultTable.minimumMeanExpression = ...
            repmat(minimumMeanExpression, height(resultTable), 1);

        resultTable.sample_support_filter_enabled = ...
            repmat(useSampleSupportFiltering, height(resultTable), 1);
        resultTable.minimumNumberOfExpressingSamples = ...
            repmat(minimumNumberOfExpressingSamples, height(resultTable), 1);
        resultTable.sampleSupportExpressionThreshold = ...
            repmat(sampleSupportExpressionThreshold, height(resultTable), 1);

        resultTable.variance_filter_enabled = ...
            repmat(useVarianceFiltering, height(resultTable), 1);
        resultTable.minimumVarianceAcrossSamples = ...
            repmat(minimumVarianceAcrossSamples, height(resultTable), 1);

        resultTable.poisson_noise_filter_enabled = ...
            repmat(usePoissonNoiseFiltering, height(resultTable), 1);
        resultTable.poissonNoiseZCutoff = ...
            repmat(poissonNoiseZ, height(resultTable), 1);

        resultTable.pvalue = pValues;
        resultTable.t_statistic = tStats;
        resultTable.degrees_of_freedom = degreesOfFreedom;

        resultTable.rnaseqde_Mean1_WT = rnaseqdeMean1;
        resultTable.rnaseqde_Mean2_5xFAD = rnaseqdeMean2;
        resultTable.rnaseqde_Log2FoldChange = rnaseqdeLog2FoldChange;
        resultTable.rnaseqde_AdjustedPValue = rnaseqdeAdjustedPValue;

        resultTable.number_of_WT_mice = ...
            repmat(numel(wtSampleNames), height(resultTable), 1);
        resultTable.number_of_5xFAD_mice = ...
            repmat(numel(fadSampleNames), height(resultTable), 1);

        resultTable.statistical_test = ...
            repmat(statisticalTest, height(resultTable), 1);
        resultTable.rnaseqdeVarianceLink = ...
            repmat(rnaseqdeVarianceLink, height(resultTable), 1);
        resultTable.rnaseqdeFDRMethod = ...
            repmat(rnaseqdeFDRMethod, height(resultTable), 1);
        resultTable.pvalue_mode = ...
            repmat(pValueMode, height(resultTable), 1);
        resultTable.number_of_permutations = ...
            repmat(numberOfPermutationsUsed, height(resultTable), 1);

        allResultsTable = [allResultsTable; resultTable];
    end
end


function rawCountTableForRnaseqde = buildRnaseqdeCountTable(groupRawTable, sampleNames)

    % rnaseqde requires integer count data in a table.
    % This function uses raw counts and forces sample count columns to
    % nonnegative integer doubles.

    rawCountTableForRnaseqde = table();
    rawCountTableForRnaseqde.gene = groupRawTable.gene;

    for s = 1:numel(sampleNames)

        sampleName = char(sampleNames(s));
        counts = groupRawTable.(sampleName);

        counts = double(counts);
        counts(~isfinite(counts)) = 0;
        counts(counts < 0) = 0;

        % Should already be integers because these are transcript counts.
        counts = round(counts);

        rawCountTableForRnaseqde.(sampleName) = counts;
    end
end


function [pValues, tStats, degreesOfFreedom] = welchTTestRows(group1, group2)

    group1 = double(group1);
    group2 = double(group2);

    n1 = size(group1, 2);
    n2 = size(group2, 2);

    mean1 = mean(group1, 2, "omitnan");
    mean2 = mean(group2, 2, "omitnan");

    var1 = var(group1, 0, 2, "omitnan");
    var2 = var(group2, 0, 2, "omitnan");

    se = sqrt(var1 ./ n1 + var2 ./ n2);

    tStats = nan(size(mean1));
    degreesOfFreedom = nan(size(mean1));
    pValues = nan(size(mean1));

    valid = se > 0 & isfinite(se);

    tStats(valid) = (mean2(valid) - mean1(valid)) ./ se(valid);

    numerator = (var1 ./ n1 + var2 ./ n2) .^ 2;

    denominator = ...
        ((var1 ./ n1) .^ 2) ./ (n1 - 1) + ...
        ((var2 ./ n2) .^ 2) ./ (n2 - 1);

    degreesOfFreedom(valid) = numerator(valid) ./ denominator(valid);

    pValues(valid) = twoSidedTPValue(tStats(valid), degreesOfFreedom(valid));

    zeroVariance = se == 0 & isfinite(mean1) & isfinite(mean2);

    sameMean = zeroVariance & mean1 == mean2;
    differentMean = zeroVariance & mean1 ~= mean2;

    pValues(sameMean) = 1;
    tStats(sameMean) = 0;
    degreesOfFreedom(sameMean) = n1 + n2 - 2;

    pValues(differentMean) = eps;
    tStats(differentMean) = sign(mean2(differentMean) - mean1(differentMean)) .* Inf;
    degreesOfFreedom(differentMean) = n1 + n2 - 2;
end


function p = twoSidedTPValue(tStats, degreesOfFreedom)

    tStats = double(tStats);
    degreesOfFreedom = double(degreesOfFreedom);

    p = nan(size(tStats));

    valid = isfinite(tStats) & isfinite(degreesOfFreedom) & degreesOfFreedom > 0;

    if ~any(valid)
        return;
    end

    tAbs = abs(tStats(valid));
    df = degreesOfFreedom(valid);

    if exist("tcdf", "file") == 2
        p(valid) = 2 .* (1 - tcdf(tAbs, df));
    else
        x = df ./ (df + tAbs.^2);
        p(valid) = betainc(x, df ./ 2, 0.5);
    end

    p(p < 0) = 0;
    p(p > 1) = 1;
end


function permDiffs = calculatePermutationDiffs( ...
    expressionMatrix, ...
    nWT, ...
    useExhaustivePermutations, ...
    number_of_realizations)

    nSamples = size(expressionMatrix, 2);

    if useExhaustivePermutations

        wtCombinations = nchoosek(1:nSamples, nWT);
        nPermutations = size(wtCombinations, 1);

        permDiffs = zeros(size(expressionMatrix, 1), nPermutations);

        for p = 1:nPermutations

            wtIdx = wtCombinations(p, :);
            fadIdx = setdiff(1:nSamples, wtIdx);

            permWT = expressionMatrix(:, wtIdx);
            permFAD = expressionMatrix(:, fadIdx);

            permDiffs(:, p) = mean(permFAD, 2) - mean(permWT, 2);
        end

    else

        permDiffs = zeros(size(expressionMatrix, 1), number_of_realizations);

        for p = 1:number_of_realizations

            shuffledIdx = randperm(nSamples);

            wtIdx = shuffledIdx(1:nWT);
            fadIdx = shuffledIdx(nWT + 1:end);

            permWT = expressionMatrix(:, wtIdx);
            permFAD = expressionMatrix(:, fadIdx);

            permDiffs(:, p) = mean(permFAD, 2) - mean(permWT, 2);
        end
    end
end


function pValues = calculatePermutationPValues(observedDiff, permDiffs, pValueMode)

    observedDiff = double(observedDiff);

    pValues = nan(size(observedDiff));

    switch string(pValueMode)

        case "normal_fit"

            nullMean = mean(permDiffs, 2);
            nullStd = std(permDiffs, 0, 2);

            valid = nullStd > 0 & isfinite(nullStd);

            zScores = nan(size(observedDiff));
            zScores(valid) = (observedDiff(valid) - nullMean(valid)) ./ nullStd(valid);

            pValues(valid) = erfc(abs(zScores(valid)) ./ sqrt(2));

        case "empirical"

            nPermutations = size(permDiffs, 2);

            for i = 1:numel(observedDiff)

                if ~isfinite(observedDiff(i))
                    pValues(i) = NaN;
                    continue;
                end

                pValues(i) = ...
                    (sum(abs(permDiffs(i, :)) >= abs(observedDiff(i))) + 1) ./ ...
                    (nPermutations + 1);
            end

        otherwise

            error("Unknown pValueMode: %s. Use 'normal_fit' or 'empirical'.", pValueMode);
    end
end


function FDR = bhFDR(pValues)

    pValues = double(pValues);
    FDR = nan(size(pValues));

    valid = isfinite(pValues);

    if ~any(valid)
        return;
    end

    validIdx = find(valid);
    p = pValues(valid);

    [pSorted, sortOrder] = sort(p, "ascend");

    m = numel(pSorted);

    qSorted = pSorted .* m ./ transpose(1:m);

    qSorted = flipud(cummin(flipud(qSorted)));

    qSorted(qSorted > 1) = 1;

    q = nan(size(p));
    q(sortOrder) = qSorted;

    FDR(validIdx) = q;
end

function compactTable = makeCompactSignificantResultsTable( ...
    significantResultsTable, ...
    sampleNames, ...
    useCellTypes, ...
    printIndividualSampleCountsInTerminal)

if height(significantResultsTable) == 0
    compactTable = significantResultsTable;
    return;
end

%% Core identifier columns
if useCellTypes
    baseCols = {'gene', 'region_name', 'cell_type'};
else
    baseCols = {'gene', 'region_name'};
end

%% Main summary columns
% Removed:
%   diff_5xFAD_minus_WT
%   rnaseqde_Mean1_WT and all rnaseqde_* columns
summaryCols = { ...
    'mean_WT_raw', ...
    'mean_5xFAD_raw', ...
    'mean_WT_analysis', ...
    'mean_5xFAD_analysis', ...
    'foldChange_5xFAD_over_WT_plusPseudo', ...
    'log2FC_5xFAD_over_WT_plusPseudo', ...
    'pvalue', ...
    'FDR'};

%% Optional individual mouse counts
individualCountCols = {};

if printIndividualSampleCountsInTerminal

    for s = 1:numel(sampleNames)

        sampleName = char(sampleNames(s));

        rawCol = [sampleName '_raw'];
        analysisCol = [sampleName '_analysis'];

        if ismember(rawCol, significantResultsTable.Properties.VariableNames)
            individualCountCols = [individualCountCols, {rawCol}];
        end

        if ismember(analysisCol, significantResultsTable.Properties.VariableNames)
            individualCountCols = [individualCountCols, {analysisCol}];
        end
    end
end

requestedCols = [baseCols, individualCountCols, summaryCols];

existingCols = requestedCols( ...
    ismember(requestedCols, significantResultsTable.Properties.VariableNames));

compactTable = significantResultsTable(:, existingCols);

if ismember('FDR', compactTable.Properties.VariableNames)
    compactTable = sortrows(compactTable, 'FDR', 'ascend');
elseif ismember('pvalue', compactTable.Properties.VariableNames)
    compactTable = sortrows(compactTable, 'pvalue', 'ascend');
end
end

function significantGenesPerRegion = summarizeSignificantGenesPerRegion( ...
    significantResultsTable, ...
    useCellTypes)

    if height(significantResultsTable) == 0

        if useCellTypes
            significantGenesPerRegion = table( ...
                strings(0,1), ...
                strings(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                nan(0,1), ...
                nan(0,1), ...
                'VariableNames', {'region_name', ...
                                  'cell_type', ...
                                  'number_significant_rows', ...
                                  'number_unique_significant_genes', ...
                                  'number_higher_in_5xFAD', ...
                                  'number_higher_in_WT', ...
                                  'minimum_pvalue', ...
                                  'minimum_FDR'});
        else
            significantGenesPerRegion = table( ...
                strings(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                zeros(0,1), ...
                nan(0,1), ...
                nan(0,1), ...
                'VariableNames', {'region_name', ...
                                  'number_significant_rows', ...
                                  'number_unique_significant_genes', ...
                                  'number_higher_in_5xFAD', ...
                                  'number_higher_in_WT', ...
                                  'minimum_pvalue', ...
                                  'minimum_FDR'});
        end

        return;
    end

    significantGenesPerRegion = table();

    if useCellTypes
        [G, regionList, cellTypeList] = findgroups( ...
            significantResultsTable.region_name, ...
            significantResultsTable.cell_type);
    else
        [G, regionList] = findgroups(significantResultsTable.region_name);
    end

    for i = 1:max(G)

        rows = G == i;

        regionTable = significantResultsTable(rows, :);

        numberSignificantRows = height(regionTable);
        numberUniqueSignificantGenes = numel(unique(regionTable.gene));

        if ismember("log2FC_5xFAD_over_WT_plusPseudo", ...
                string(regionTable.Properties.VariableNames))

            log2fc = regionTable.log2FC_5xFAD_over_WT_plusPseudo;

            numberHigherIn5xFAD = sum(log2fc > 0);
            numberHigherInWT = sum(log2fc < 0);

        else

            numberHigherIn5xFAD = NaN;
            numberHigherInWT = NaN;

        end

        minimumPValue = min(regionTable.pvalue, [], "omitnan");
        minimumFDR = min(regionTable.FDR, [], "omitnan");

        if useCellTypes

            newRow = table( ...
                regionList(i), ...
                cellTypeList(i), ...
                numberSignificantRows, ...
                numberUniqueSignificantGenes, ...
                numberHigherIn5xFAD, ...
                numberHigherInWT, ...
                minimumPValue, ...
                minimumFDR, ...
                'VariableNames', {'region_name', ...
                                  'cell_type', ...
                                  'number_significant_rows', ...
                                  'number_unique_significant_genes', ...
                                  'number_higher_in_5xFAD', ...
                                  'number_higher_in_WT', ...
                                  'minimum_pvalue', ...
                                  'minimum_FDR'});

        else

            newRow = table( ...
                regionList(i), ...
                numberSignificantRows, ...
                numberUniqueSignificantGenes, ...
                numberHigherIn5xFAD, ...
                numberHigherInWT, ...
                minimumPValue, ...
                minimumFDR, ...
                'VariableNames', {'region_name', ...
                                  'number_significant_rows', ...
                                  'number_unique_significant_genes', ...
                                  'number_higher_in_5xFAD', ...
                                  'number_higher_in_WT', ...
                                  'minimum_pvalue', ...
                                  'minimum_FDR'});
        end

        significantGenesPerRegion = [significantGenesPerRegion; newRow];
    end

    significantGenesPerRegion = sortrows(significantGenesPerRegion, ...
        "number_unique_significant_genes", ...
        "descend");
end
