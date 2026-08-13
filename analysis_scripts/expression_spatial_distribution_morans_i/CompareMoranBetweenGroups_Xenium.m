%% Compare animal-averaged Morans I between WT and 5xFAD mice by hippocampal region
% Each CSV contains transcript-level rows with at least:
%   gene
%   region_name
%   global_x_um
%   global_y_um
%   Z_um
%
% The script:
%   1. Reads the user-defined sample list (see EXAMPLE sample list below)
%   2. Uses shared gene-name aliases and exclusions
%   3. Uses hippocampal region remapping and exclusions
%   4. Calculates Morans I per gene, region, and physical section
%   5. Averages consecutive sections from the same animal before inference
%   6. Aggregates RNA molecules into 10 um cubic voxels
%   7. Uses six-neighbor voxel connectivity: +/-x, +/-y, +/-z
%   8. Optionally quantile-normalizes Morans I values across animals
%   9. Compares WT vs 5xFAD per gene-region pair
%   10. Supports:
%        - empirical-Bayes moderated t-test, limma-like
%        - Welch t-test
%        - permutation test
%   11. Applies global BH-FDR correction across all tested gene-region pairs
%
% Important moderated t-test detail:
%   The empirical-Bayes prior variance is estimated from all valid
%   gene-region pairs, before applying the effect-size filters.
%   The statistical test itself is then performed only on the filtered
%   gene-region pairs.
%
% Direction convention:
%   diff_5xFAD_minus_WT = mean(5xFAD) - mean(WT)

clear; clc; close all;

%% =============================
% USER SETTINGS
% =============================

scriptDir = fileparts(mfilename('fullpath'));
inputFolder  = fullfile(scriptDir, '..', '..', '..', 'data', 'xenium');   % repo data/xenium (from Zenodo)
outputFolder = fullfile(scriptDir, 'results', 'xenium_moran');

if ~exist(outputFolder, "dir")
    mkdir(outputFolder);
end

%% =============================
% Morans I settings
% =============================

voxelSizeUm = 10;

% Minimum number of RNA molecules for a gene in a region in a given animal
% required to calculate Morans I. Below this, Morans I is set to NaN.
minimumMoleculesPerGeneRegionAnimal = 3;

% Safety limit for the number of voxels in one region of one animal.
% If this is exceeded, increase voxelSizeUm or check coordinate scale.
maximumNumberOfVoxelsPerRegion = 2000000;

%% =============================
% Normalization settings
% =============================

% Quantile normalization of Morans I values across animals
useQuantileNormalization = false;

%% =============================
% Statistical test settings
% =============================

% Options:
% "moderated_ttest" = empirical-Bayes moderated t-test, limma-like
% "welch_ttest"     = Welch t-test
% "permutation"     = permutation test with normal-fit p value
statisticalTest = "moderated_ttest";

% Moderated t-test settings
moderatedTtestPriorDF = 10;

% Options:
% "median" = prior variance is median of observed gene-region variances
% "mean"   = prior variance is mean of observed gene-region variances
moderatedTtestPriorVarianceMode = "median";

% Optional transform applied before statistical testing.
% For Morans I, "none" is usually safest.
% Options:
% "none"
% "atanh"  = Fisher-like transform after clipping values to (-1,1)
moranTestTransform = "none";

% Permutation settings, used only if statisticalTest = "permutation"
number_of_realizations = 1000;

% With 2 WT and 2 5xFAD biological animals after section averaging,
% exhaustive permutation gives nchoosek(4,2)=6.
% true  = use all possible label permutations
% false = use random permutations
useExhaustivePermutations = true;

%% =============================
% Filtering before statistical testing
% =============================

% Require at least this many valid Morans I values per group
minimumValidSamplesPerGroup = 2;

% Filter by absolute mean difference between 5xFAD and WT.
% The threshold is calculated from all valid gene-region pairs.
useMeanDifferenceQuantileFilter = false;
meanDifferenceQuantile = 0.80;

% Absolute Morans I difference filter:
% Only gene-region pairs with abs(mean_5xFAD - mean_WT) >= this cutoff
% after Morans I quantile normalization will be tested.
useAbsDiffNormalizedFilter = true;
absDiffNormalizedCutoff = 0.02;

% Mean Morans I magnitude filter:
% Only gene-region pairs where either group has mean normalized Morans I
% above this threshold will be tested.
%
% Rule:
%   mean_5xFAD_Moran_normalized >= 0.1 OR mean_WT_Moran_normalized >= 0.1
useMeanMoranNormalizedFilter = false;
meanMoranNormalizedCutoff = 0.15;

%% =============================
% Reporting thresholds
% =============================

cutoff_p = 0.05;
cutoff_FDR = 0.15;

% If true, significant table requires P < cutoff_p.
% If false, significant table requires FDR < cutoff_FDR.
applyPCutoffForReporting = false;

%% =============================
% Sample list (published Xenium names)
% Required columns: gene, region_name, global_x_um, global_y_um, Z_um
% (aliases such as x/y/z or region are also accepted).
%
% Morans I is calculated per physical section listed below, then sections
% from the same animal are averaged in averageSectionsToAnimals before
% group comparison. Groups must be "WT" or "FAD". After averaging, keep
% equal animal counts in each group (here 2+2).
% Place renamed Archive CSV files in data/xenium/.
% =============================

samples = struct();

% WT_6: two sections averaged to animal "WT_6"
samples(1).name = "WT_6_section1";
samples(1).group = "WT";
samples(1).files = {
    "WT_6_section1.csv"
    };

samples(2).name = "WT_6_section2";
samples(2).group = "WT";
samples(2).files = {
    "WT_6_section2.csv"
    };

% WT_5: single section
samples(3).name = "WT_5";
samples(3).group = "WT";
samples(3).files = {
    "WT_5.csv"
    };

% 5xFAD_5: two sections averaged to animal "5xFAD_5"
samples(4).name = "5xFAD_5_section1";
samples(4).group = "FAD";
samples(4).files = {
    "5xFAD_5_section1.csv"
    };

samples(5).name = "5xFAD_5_section2";
samples(5).group = "FAD";
samples(5).files = {
    "5xFAD_5_section2.csv"
    };

% 5xFAD_6: single section
samples(6).name = "5xFAD_6";
samples(6).group = "FAD";
samples(6).files = {
    "5xFAD_6.csv"
    };

sectionSampleNames = string({samples.name});

% Animal-level names after section averaging (must match averageSectionsToAnimals)
sampleNames = ["WT_6", "WT_5", "5xFAD_5", "5xFAD_6"];
sampleGroups = ["WT", "WT", "FAD", "FAD"];

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
% Calculate Morans I per animal
% =============================

allMoranLong = table();

for i = 1:numel(samples)

    fprintf("\n=============================\n");
    fprintf("Processing sample %s\n", samples(i).name);
    fprintf("=============================\n");

    moleculeTable = readAndCleanMoranInputFiles( ...
        samples(i).files, ...
        inputFolder);

    sampleMoranTable = calculateMoranForOneAnimal( ...
        moleculeTable, ...
        samples(i).name, ...
        samples(i).group, ...
        voxelSizeUm, ...
        minimumMoleculesPerGeneRegionAnimal, ...
        maximumNumberOfVoxelsPerRegion);

    allMoranLong = [allMoranLong; sampleMoranTable];
end

writetable(allMoranLong, ...
    fullfile(outputFolder, "Moran_I_per_gene_region_animal_long.xlsx"));

%% =============================
% Build wide Moran table
% =============================

rawMoranSectionWideTable = buildWideMoranTable(allMoranLong, sectionSampleNames);

writetable(rawMoranSectionWideTable, ...
    fullfile(outputFolder, "Moran_I_raw_wide_by_section.xlsx"));

rawMoranWideTable = averageSectionsToAnimals(rawMoranSectionWideTable);

writetable(rawMoranWideTable, ...
    fullfile(outputFolder, "Moran_I_raw_wide_by_animal.xlsx"));

fprintf("\n=============================\n");
fprintf("Raw Morans I table\n");
fprintf("=============================\n\n");

fprintf("Number of gene-region pairs: %d\n", height(rawMoranWideTable));

%% =============================
% Quantile-normalize Morans I across animals
% =============================

rawMoranMatrix = rawMoranWideTable{:, cellstr(sampleNames)};

if useQuantileNormalization

    normalizedMoranMatrix = quantileNormalizeMatrixWithNaNs(rawMoranMatrix);

    fprintf("\nQuantile normalization was applied to Morans I values across animals.\n");

else

    normalizedMoranMatrix = rawMoranMatrix;

    fprintf("\nQuantile normalization was not applied.\n");
end

normalizedMoranWideTable = rawMoranWideTable(:, {'gene', 'region_name'});

for s = 1:numel(sampleNames)
    normalizedMoranWideTable.(char(sampleNames(s))) = normalizedMoranMatrix(:, s);
end

writetable(normalizedMoranWideTable, ...
    fullfile(outputFolder, "Moran_I_normalized_wide.xlsx"));

%% =============================
% Compare WT vs 5xFAD
% =============================

allResultsTable = compareMoranBetweenGroups( ...
    rawMoranWideTable, ...
    normalizedMoranWideTable, ...
    sampleNames, ...
    wtSampleNames, ...
    fadSampleNames, ...
    statisticalTest, ...
    moderatedTtestPriorDF, ...
    moderatedTtestPriorVarianceMode, ...
    moranTestTransform, ...
    useExhaustivePermutations, ...
    number_of_realizations, ...
    minimumValidSamplesPerGroup, ...
    useMeanDifferenceQuantileFilter, ...
    meanDifferenceQuantile, ...
    useAbsDiffNormalizedFilter, ...
    absDiffNormalizedCutoff, ...
    useMeanMoranNormalizedFilter, ...
    meanMoranNormalizedCutoff);

if height(allResultsTable) == 0
    error("No gene-region pairs passed the filters.");
end

allResultsTable.FDR = bhFDR(allResultsTable.pvalue);

if applyPCutoffForReporting
    significantRows = isfinite(allResultsTable.pvalue) & ...
        allResultsTable.pvalue < cutoff_p;
    reportingMode = sprintf("P-value only: P < %.4g", cutoff_p);
else
    significantRows = isfinite(allResultsTable.FDR) & ...
        allResultsTable.FDR < cutoff_FDR;
    reportingMode = sprintf("FDR only: FDR < %.4g", cutoff_FDR);
end

significantResultsTable = allResultsTable(significantRows, :);

if applyPCutoffForReporting
    significantResultsTable = sortrows(significantResultsTable, "pvalue", "ascend");
else
    significantResultsTable = sortrows(significantResultsTable, "FDR", "ascend");
end

writetable(allResultsTable, ...
    fullfile(outputFolder, "all_WT_vs_5xFAD_Moran_results.xlsx"));

writetable(significantResultsTable, ...
    fullfile(outputFolder, "significant_WT_vs_5xFAD_Moran_results.xlsx"));

%% =============================
% Compact terminal display
% =============================

compactSignificantTable = makeCompactSignificantMoranTable(significantResultsTable);

fprintf("\n=============================\n");
fprintf("Reporting mode\n");
fprintf("=============================\n");
fprintf("%s\n", reportingMode);

fprintf("\nStatistical test: %s\n", statisticalTest);

if statisticalTest == "moderated_ttest"
    fprintf("Moderated t-test prior DF: %.3f\n", moderatedTtestPriorDF);
    fprintf("Prior variance mode: %s\n", moderatedTtestPriorVarianceMode);
end

fprintf("Moran test transform: %s\n", moranTestTransform);
fprintf("Minimum valid samples per group: %d\n", minimumValidSamplesPerGroup);

fprintf("\nNumber of tested gene-region pairs: %d\n", height(allResultsTable));
fprintf("Number of finite P values: %d\n", sum(isfinite(allResultsTable.pvalue)));
fprintf("Number of P < %.4g: %d\n", cutoff_p, sum(isfinite(allResultsTable.pvalue) & allResultsTable.pvalue < cutoff_p));
fprintf("Number of FDR < %.4g: %d\n", cutoff_FDR, sum(isfinite(allResultsTable.FDR) & allResultsTable.FDR < cutoff_FDR));
fprintf("Number of reported significant results: %d\n", height(significantResultsTable));

fprintf("\n=============================\n");
fprintf("Significant Morans I results, compact view\n");
fprintf("=============================\n\n");

disp(compactSignificantTable);

writetable(compactSignificantTable, ...
    fullfile(outputFolder, "significant_WT_vs_5xFAD_Moran_results_compact.xlsx"));

%% =============================
% Significant genes per region
% =============================

significantGenesPerRegion = summarizeSignificantMoranGenesPerRegion(significantResultsTable);

fprintf("\n=============================\n");
fprintf("Significant Morans I genes per region\n");
fprintf("=============================\n\n");

disp(significantGenesPerRegion);

writetable(significantGenesPerRegion, ...
    fullfile(outputFolder, "significant_Moran_genes_per_region_summary.xlsx"));

%% =============================
% Save workspace
% =============================

save(fullfile(outputFolder, "WT_vs_5xFAD_Moran_workspace.mat"), ...
    "allMoranLong", ...
    "rawMoranWideTable", ...
    "normalizedMoranWideTable", ...
    "allResultsTable", ...
    "significantResultsTable", ...
    "compactSignificantTable", ...
    "significantGenesPerRegion", ...
    "samples", ...
    "sampleNames", ...
    "sampleGroups", ...
    "wtSampleNames", ...
    "fadSampleNames", ...
    "voxelSizeUm", ...
    "minimumMoleculesPerGeneRegionAnimal", ...
    "useQuantileNormalization", ...
    "statisticalTest", ...
    "moderatedTtestPriorDF", ...
    "moderatedTtestPriorVarianceMode", ...
    "moranTestTransform", ...
    "useExhaustivePermutations", ...
    "number_of_realizations", ...
    "minimumValidSamplesPerGroup", ...
    "useMeanDifferenceQuantileFilter", ...
    "meanDifferenceQuantile", ...
    "useAbsDiffNormalizedFilter", ...
    "absDiffNormalizedCutoff", ...
    "useMeanMoranNormalizedFilter", ...
    "meanMoranNormalizedCutoff", ...
    "cutoff_p", ...
    "cutoff_FDR", ...
    "applyPCutoffForReporting", ...
    "reportingMode");

fprintf("\nAnalysis complete. Results saved to folder:\n%s\n", outputFolder);

%% ========================================================================
% Local functions
% ========================================================================

function T = readAndCleanMoranInputFiles(fileList, inputFolder)

    allRows = table();

    for f = 1:numel(fileList)

        currentFile = string(fileList{f});
        currentPath = fullfile(inputFolder, currentFile);

        if ~isfile(currentPath)
            error("File not found: %s", currentPath);
        end

        fprintf("Reading %s\n", currentFile);

        Traw = readtable(currentPath, ...
            "TextType", "string", ...
            "VariableNamingRule", "preserve");

        varNames = string(Traw.Properties.VariableNames);

        geneCol = findFirstColumn(varNames, ["gene", "genes", "gene_name", "Gene"]);
        regionCol = findFirstColumn(varNames, ["region_name", "region", "hippocampal_region"]);
        xCol = findFirstColumn(varNames, ["global_x_um", "global_x", "x_um", "x"]);
        yCol = findFirstColumn(varNames, ["global_y_um", "global_y", "y_um", "y"]);
        zCol = findFirstColumn(varNames, ["Z_um", "Z", "global_z_um", "z_um", "z"]);

        if strlength(geneCol) == 0
            error("Could not find gene column in %s", currentFile);
        end

        if strlength(regionCol) == 0
            error("Could not find region column in %s", currentFile);
        end

        if strlength(xCol) == 0 || strlength(yCol) == 0 || strlength(zCol) == 0
            error("Could not find global_x_um, global_y_um, and Z_um columns in %s", currentFile);
        end

        genes = standardizeGeneSymbol(Traw.(geneCol));
        [regions, excludeRegionRows] = standardizeRegionName(Traw.(regionCol));

        x = toDoubleVector(Traw.(xCol));
        y = toDoubleVector(Traw.(yCol));
        z = toDoubleVector(Traw.(zCol));

        excludeGeneRows = isExcludedGene(genes);

        keepRows = ...
            ~ismissing(genes) & strlength(genes) > 0 & ...
            ~ismissing(regions) & strlength(regions) > 0 & ...
            ~excludeGeneRows & ...
            ~excludeRegionRows & ...
            isfinite(x) & isfinite(y) & isfinite(z);

        Tclean = table( ...
            genes(keepRows), ...
            regions(keepRows), ...
            x(keepRows), ...
            y(keepRows), ...
            z(keepRows), ...
            'VariableNames', {'gene', ...
                              'region_name', ...
                              'global_x_um', ...
                              'global_y_um', ...
                              'Z_um'});

        allRows = [allRows; Tclean];
    end

    T = allRows;
end


function sampleMoranTable = calculateMoranForOneAnimal( ...
    moleculeTable, ...
    sampleName, ...
    sampleGroup, ...
    voxelSizeUm, ...
    minimumMoleculesPerGeneRegionAnimal, ...
    maximumNumberOfVoxelsPerRegion)

    sampleMoranTable = table();

    if height(moleculeTable) == 0
        warning("No valid molecules for sample %s", sampleName);
        return;
    end

    regionList = unique(moleculeTable.region_name);

    for r = 1:numel(regionList)

        regionName = regionList(r);

        regionRows = moleculeTable.region_name == regionName;
        regionTable = moleculeTable(regionRows, :);

        fprintf("  Region %s: %d molecules\n", regionName, height(regionTable));

        regionMoranTable = calculateMoranForOneRegion( ...
            regionTable, ...
            sampleName, ...
            sampleGroup, ...
            regionName, ...
            voxelSizeUm, ...
            minimumMoleculesPerGeneRegionAnimal, ...
            maximumNumberOfVoxelsPerRegion);

        sampleMoranTable = [sampleMoranTable; regionMoranTable];
    end
end


function regionMoranTable = calculateMoranForOneRegion( ...
    regionTable, ...
    sampleName, ...
    sampleGroup, ...
    regionName, ...
    voxelSizeUm, ...
    minimumMoleculesPerGeneRegionAnimal, ...
    maximumNumberOfVoxelsPerRegion)

    coords = [regionTable.global_x_um, regionTable.global_y_um, regionTable.Z_um];

    minCoords = min(coords, [], 1);
    voxelCoords = floor((coords - minCoords) ./ voxelSizeUm) + 1;

    voxelCoords = max(voxelCoords, 1);

    gridDims = max(voxelCoords, [], 1);
    gridDims = double(gridDims);

    numberOfVoxels = prod(gridDims);

    if numberOfVoxels > maximumNumberOfVoxelsPerRegion
        warning("Skipping region %s in sample %s because it requires %d voxels, above maximumNumberOfVoxelsPerRegion = %d.", ...
            regionName, sampleName, numberOfVoxels, maximumNumberOfVoxelsPerRegion);
        regionMoranTable = table();
        return;
    end

    linearVoxelIndex = sub2ind(gridDims, ...
        voxelCoords(:,1), ...
        voxelCoords(:,2), ...
        voxelCoords(:,3));

    geneList = unique(regionTable.gene);

    regionMoranTable = table();

    for g = 1:numel(geneList)

        geneName = geneList(g);

        geneRows = regionTable.gene == geneName;

        numberOfMolecules = sum(geneRows);

        if numberOfMolecules < minimumMoleculesPerGeneRegionAnimal

            moranI = NaN;
            numberNonzeroVoxels = NaN;

        else

            geneVoxelIndex = linearVoxelIndex(geneRows);

            counts = accumarray( ...
                geneVoxelIndex, ...
                1, ...
                [numberOfVoxels, 1], ...
                @sum, ...
                0);

            numberNonzeroVoxels = sum(counts > 0);

            moranI = computeMoranIFromVoxelCounts(counts, gridDims);
        end

        newRow = table( ...
            string(sampleName), ...
            string(sampleGroup), ...
            geneName, ...
            regionName, ...
            moranI, ...
            numberOfMolecules, ...
            numberOfVoxels, ...
            numberNonzeroVoxels, ...
            voxelSizeUm, ...
            'VariableNames', {'sample', ...
                              'group', ...
                              'gene', ...
                              'region_name', ...
                              'Moran_I', ...
                              'number_molecules', ...
                              'number_voxels', ...
                              'number_nonzero_voxels', ...
                              'voxel_size_um'});

        regionMoranTable = [regionMoranTable; newRow];
    end
end


function moranI = computeMoranIFromVoxelCounts(counts, gridDims)

    counts = double(counts);

    N = numel(counts);

    if N < 2
        moranI = NaN;
        return;
    end

    xbar = mean(counts);
    dev = counts - xbar;

    denominator = sum(dev .^ 2);

    if denominator <= 0 || ~isfinite(denominator)
        moranI = NaN;
        return;
    end

    D = reshape(dev, gridDims);

    numerator = 0;
    numberDirectedWeights = 0;

    if gridDims(1) > 1
        pairProducts = D(1:end-1,:,:) .* D(2:end,:,:);
        numerator = numerator + 2 .* sum(pairProducts(:));
        numberDirectedWeights = numberDirectedWeights + 2 .* numel(pairProducts);
    end

    if gridDims(2) > 1
        pairProducts = D(:,1:end-1,:) .* D(:,2:end,:);
        numerator = numerator + 2 .* sum(pairProducts(:));
        numberDirectedWeights = numberDirectedWeights + 2 .* numel(pairProducts);
    end

    if gridDims(3) > 1
        pairProducts = D(:,:,1:end-1) .* D(:,:,2:end);
        numerator = numerator + 2 .* sum(pairProducts(:));
        numberDirectedWeights = numberDirectedWeights + 2 .* numel(pairProducts);
    end

    if numberDirectedWeights == 0
        moranI = NaN;
        return;
    end

    moranI = (N ./ numberDirectedWeights) .* (numerator ./ denominator);
end


function rawMoranWideTable = buildWideMoranTable(allMoranLong, sampleNames)

    keyVars = {'gene', 'region_name'};

    rawMoranWideTable = table();

    for s = 1:numel(sampleNames)

        sampleName = sampleNames(s);

        sampleRows = allMoranLong.sample == sampleName;

        sampleTable = allMoranLong(sampleRows, {'gene', 'region_name', 'Moran_I'});

        sampleTable = renamevars(sampleTable, "Moran_I", sampleName);

        if s == 1
            rawMoranWideTable = sampleTable;
        else
            rawMoranWideTable = outerjoin( ...
                rawMoranWideTable, ...
                sampleTable, ...
                "Keys", keyVars, ...
                "MergeKeys", true);
        end
    end

    rawMoranWideTable = sortrows(rawMoranWideTable, keyVars);
end


function animalWideTable = averageSectionsToAnimals(sectionWideTable)
% Average multi-section animals to one column each.
% Edit this helper if the sample list uses different section names.

    animalWideTable = sectionWideTable(:, {'gene', 'region_name'});

    animalWideTable.("WT_6") = averageAvailableColumns( ...
        sectionWideTable, ["WT_6_section1", "WT_6_section2"]);
    animalWideTable.("WT_5") = sectionWideTable.("WT_5");
    animalWideTable.("5xFAD_5") = averageAvailableColumns( ...
        sectionWideTable, ["5xFAD_5_section1", "5xFAD_5_section2"]);
    animalWideTable.("5xFAD_6") = sectionWideTable.("5xFAD_6");
end


function values = averageAvailableColumns(T, columnNames)

    matrix = T{:, cellstr(columnNames)};
    values = mean(matrix, 2, "omitnan");
    values(all(~isfinite(matrix), 2)) = NaN;
end


function normalizedMatrix = quantileNormalizeMatrixWithNaNs(rawMatrix)

    rawMatrix = double(rawMatrix);

    nanMask = ~isfinite(rawMatrix);

    filledMatrix = rawMatrix;

    for c = 1:size(filledMatrix, 2)

        col = filledMatrix(:, c);
        finiteCol = col(isfinite(col));

        if isempty(finiteCol)
            fillValue = 0;
        else
            fillValue = median(finiteCol);
        end

        col(~isfinite(col)) = fillValue;
        filledMatrix(:, c) = col;
    end

    [sortedValues, sortIdx] = sort(filledMatrix, 1, "ascend");

    meanSortedValues = mean(sortedValues, 2, "omitnan");

    normalizedMatrix = zeros(size(filledMatrix));

    for c = 1:size(filledMatrix, 2)
        normalizedMatrix(sortIdx(:, c), c) = meanSortedValues;
    end

    normalizedMatrix(nanMask) = NaN;
end


function allResultsTable = compareMoranBetweenGroups( ...
    rawMoranWideTable, ...
    normalizedMoranWideTable, ...
    sampleNames, ...
    wtSampleNames, ...
    fadSampleNames, ...
    statisticalTest, ...
    moderatedTtestPriorDF, ...
    moderatedTtestPriorVarianceMode, ...
    moranTestTransform, ...
    useExhaustivePermutations, ...
    number_of_realizations, ...
    minimumValidSamplesPerGroup, ...
    useMeanDifferenceQuantileFilter, ...
    meanDifferenceQuantile, ...
    useAbsDiffNormalizedFilter, ...
    absDiffNormalizedCutoff, ...
    useMeanMoranNormalizedFilter, ...
    meanMoranNormalizedCutoff)

    sampleCols = cellstr(sampleNames);

    rawMatrix = rawMoranWideTable{:, sampleCols};
    normalizedMatrix = normalizedMoranWideTable{:, sampleCols};

    wtIdx = find(ismember(sampleNames, wtSampleNames));
    fadIdx = find(ismember(sampleNames, fadSampleNames));

    wtRaw = rawMatrix(:, wtIdx);
    fadRaw = rawMatrix(:, fadIdx);

    wtNorm = normalizedMatrix(:, wtIdx);
    fadNorm = normalizedMatrix(:, fadIdx);

    nWTvalid = sum(isfinite(wtNorm), 2);
    nFADvalid = sum(isfinite(fadNorm), 2);

    meanWTRaw = mean(wtRaw, 2, "omitnan");
    meanFADRaw = mean(fadRaw, 2, "omitnan");

    meanWTNorm = mean(wtNorm, 2, "omitnan");
    meanFADNorm = mean(fadNorm, 2, "omitnan");

    diffNorm = meanFADNorm - meanWTNorm;
    absDiffNorm = abs(diffNorm);

    validForTesting = ...
        nWTvalid >= minimumValidSamplesPerGroup & ...
        nFADvalid >= minimumValidSamplesPerGroup & ...
        isfinite(diffNorm);

    %% Mean-difference quantile filter

    if useMeanDifferenceQuantileFilter

        diffThreshold = quantile(absDiffNorm(validForTesting), meanDifferenceQuantile);

        passMeanDifferenceQuantile = ...
            validForTesting & ...
            absDiffNorm >= diffThreshold;

    else

        diffThreshold = NaN;
        passMeanDifferenceQuantile = true(size(validForTesting));
    end

    %% Absolute normalized Morans I difference filter

    numberTestsBeforeAbsDiffFilter = sum(validForTesting);

    if useAbsDiffNormalizedFilter

        passAbsDiffNormalizedFilter = ...
            validForTesting & ...
            absDiffNorm >= absDiffNormalizedCutoff;

    else

        passAbsDiffNormalizedFilter = true(size(validForTesting));
    end

    numberTestsAfterAbsDiffFilter = sum(validForTesting & passAbsDiffNormalizedFilter);

    %% Mean normalized Morans I magnitude filter

    if useMeanMoranNormalizedFilter

        passMeanMoranNormalizedFilter = ...
            validForTesting & ...
            (meanFADNorm >= meanMoranNormalizedCutoff | ...
             meanWTNorm  >= meanMoranNormalizedCutoff);

    else

        passMeanMoranNormalizedFilter = true(size(validForTesting));
    end

    %% Combined filter

    keepRows = ...
        validForTesting & ...
        passAbsDiffNormalizedFilter & ...
        passMeanDifferenceQuantile & ...
        passMeanMoranNormalizedFilter;

    %% Print filtering diagnostics

    fprintf("\n=============================\n");
    fprintf("Morans I filtering diagnostics\n");
    fprintf("=============================\n");

    fprintf("Valid gene-region pairs before effect-size filters: %d\n", ...
        sum(validForTesting));

    if useAbsDiffNormalizedFilter
        fprintf("Passing abs(diff_5xFAD_minus_WT_normalized) >= %.4f: %d\n", ...
            absDiffNormalizedCutoff, ...
            numberTestsAfterAbsDiffFilter);

        fprintf("Removed by abs-diff filter alone: %d\n", ...
            numberTestsBeforeAbsDiffFilter - numberTestsAfterAbsDiffFilter);
    else
        fprintf("Absolute diff filter disabled.\n");
    end

    if useMeanDifferenceQuantileFilter
        fprintf("Passing mean-difference quantile filter, quantile %.2f, threshold %.6f: %d\n", ...
            meanDifferenceQuantile, ...
            diffThreshold, ...
            sum(validForTesting & passMeanDifferenceQuantile));
    else
        fprintf("Mean-difference quantile filter disabled.\n");
    end

    if useMeanMoranNormalizedFilter
        fprintf("Passing mean normalized Morans I >= %.4f in WT or 5xFAD: %d\n", ...
            meanMoranNormalizedCutoff, ...
            sum(validForTesting & passMeanMoranNormalizedFilter));
    else
        fprintf("Mean normalized Morans I filter disabled.\n");
    end

    fprintf("Passing all active pre-testing filters: %d\n", ...
        sum(keepRows));

    if ~any(keepRows)
        allResultsTable = table();
        return;
    end

    %% Transform Morans I values for testing

    testMatrix = normalizedMatrix;

    testMatrix = transformMoranMatrixForTesting(testMatrix, moranTestTransform);

    % Prior is estimated from all valid gene-region pairs.
    testMatrixForPrior = testMatrix(validForTesting, :);

    % Statistical testing is performed only on filtered gene-region pairs.
    testMatrixFiltered = testMatrix(keepRows, :);

    switch string(statisticalTest)

        case "moderated_ttest"

            priorVarianceGlobal = estimateModeratedPriorVariance( ...
                testMatrixForPrior, ...
                wtIdx, ...
                fadIdx, ...
                moderatedTtestPriorVarianceMode);

            [pValues, tStats, degreesOfFreedom, ordinaryVariance, moderatedVariance, priorVariance] = ...
                moderatedTTestRowsWithFixedPrior( ...
                    testMatrixFiltered, ...
                    wtIdx, ...
                    fadIdx, ...
                    moderatedTtestPriorDF, ...
                    priorVarianceGlobal);

            numberOfPermutationsUsed = NaN(size(pValues));

        case "welch_ttest"

            [pValues, tStats, degreesOfFreedom] = ...
                welchTTestRows(testMatrixFiltered(:, wtIdx), testMatrixFiltered(:, fadIdx));

            ordinaryVariance = nan(size(pValues));
            moderatedVariance = nan(size(pValues));
            priorVariance = nan(size(pValues));
            numberOfPermutationsUsed = NaN(size(pValues));

        case "permutation"

            [pValues, numberOfPermutationsUsed] = permutationPValuesRows( ...
                testMatrixFiltered, ...
                wtIdx, ...
                fadIdx, ...
                useExhaustivePermutations, ...
                number_of_realizations);

            tStats = nan(size(pValues));
            degreesOfFreedom = nan(size(pValues));
            ordinaryVariance = nan(size(pValues));
            moderatedVariance = nan(size(pValues));
            priorVariance = nan(size(pValues));

        otherwise

            error("Unknown statisticalTest: %s", statisticalTest);
    end

    keptRawMatrix = rawMatrix(keepRows, :);
    keptNormMatrix = normalizedMatrix(keepRows, :);

    allResultsTable = rawMoranWideTable(keepRows, {'gene', 'region_name'});

    for s = 1:numel(sampleNames)
        sampleName = char(sampleNames(s));
        allResultsTable.([sampleName '_Moran_raw']) = keptRawMatrix(:, s);
    end

    for s = 1:numel(sampleNames)
        sampleName = char(sampleNames(s));
        allResultsTable.([sampleName '_Moran_normalized']) = keptNormMatrix(:, s);
    end

    allResultsTable.mean_WT_Moran_raw = meanWTRaw(keepRows);
    allResultsTable.mean_5xFAD_Moran_raw = meanFADRaw(keepRows);

    allResultsTable.mean_WT_Moran_normalized = meanWTNorm(keepRows);
    allResultsTable.mean_5xFAD_Moran_normalized = meanFADNorm(keepRows);

    allResultsTable.diff_5xFAD_minus_WT_normalized = diffNorm(keepRows);
    allResultsTable.abs_diff_5xFAD_minus_WT_normalized = absDiffNorm(keepRows);

    allResultsTable.n_WT_valid = nWTvalid(keepRows);
    allResultsTable.n_5xFAD_valid = nFADvalid(keepRows);

    allResultsTable.mean_difference_quantile_filter_enabled = ...
        repmat(useMeanDifferenceQuantileFilter, height(allResultsTable), 1);

    allResultsTable.mean_difference_quantile = ...
        repmat(meanDifferenceQuantile, height(allResultsTable), 1);

    allResultsTable.mean_difference_quantile_threshold = ...
        repmat(diffThreshold, height(allResultsTable), 1);

    allResultsTable.abs_diff_normalized_filter_enabled = ...
        repmat(useAbsDiffNormalizedFilter, height(allResultsTable), 1);

    allResultsTable.abs_diff_normalized_cutoff = ...
        repmat(absDiffNormalizedCutoff, height(allResultsTable), 1);

    allResultsTable.mean_moran_normalized_filter_enabled = ...
        repmat(useMeanMoranNormalizedFilter, height(allResultsTable), 1);

    allResultsTable.mean_moran_normalized_cutoff = ...
        repmat(meanMoranNormalizedCutoff, height(allResultsTable), 1);

    allResultsTable.pvalue = pValues;
    allResultsTable.t_statistic = tStats;
    allResultsTable.degrees_of_freedom = degreesOfFreedom;

    allResultsTable.ordinary_variance = ordinaryVariance;
    allResultsTable.moderated_variance = moderatedVariance;
    allResultsTable.prior_variance = priorVariance;
    allResultsTable.moderated_prior_DF = repmat(moderatedTtestPriorDF, height(allResultsTable), 1);

    allResultsTable.statistical_test = repmat(string(statisticalTest), height(allResultsTable), 1);
    allResultsTable.moran_test_transform = repmat(string(moranTestTransform), height(allResultsTable), 1);
    allResultsTable.number_of_permutations = numberOfPermutationsUsed;

    allResultsTable = sortrows(allResultsTable, "pvalue", "ascend");
end


function testMatrix = transformMoranMatrixForTesting(testMatrix, moranTestTransform)

    switch string(moranTestTransform)

        case "none"

            return;

        case "atanh"

            epsilon = 1e-6;

            testMatrix(testMatrix >= 1) = 1 - epsilon;
            testMatrix(testMatrix <= -1) = -1 + epsilon;

            testMatrix = atanh(testMatrix);

        otherwise

            error("Unknown moranTestTransform: %s", moranTestTransform);
    end
end


function priorVariance = estimateModeratedPriorVariance( ...
    X, ...
    wtIdx, ...
    fadIdx, ...
    priorVarianceMode)

    Xwt = X(:, wtIdx);
    Xfad = X(:, fadIdx);

    [~, varWT, nWT] = rowMeanVarNaN(Xwt);
    [~, varFAD, nFAD] = rowMeanVarNaN(Xfad);

    residualDF = nWT + nFAD - 2;

    ordinaryVariance = ...
        ((nWT - 1) .* varWT + (nFAD - 1) .* varFAD) ./ residualDF;

    validVariance = ...
        isfinite(ordinaryVariance) & ...
        ordinaryVariance > 0 & ...
        residualDF > 0;

    switch string(priorVarianceMode)

        case "median"

            priorVariance = median(ordinaryVariance(validVariance), "omitnan");

        case "mean"

            priorVariance = mean(ordinaryVariance(validVariance), "omitnan");

        otherwise

            error("Unknown priorVarianceMode: %s", priorVarianceMode);
    end

    if ~isfinite(priorVariance) || priorVariance <= 0
        priorVariance = eps;
    end
end


function [pValues, tStats, degreesOfFreedom, ordinaryVariance, moderatedVariance, priorVarianceVector] = ...
    moderatedTTestRowsWithFixedPrior(X, wtIdx, fadIdx, priorDF, priorVariance)

    Xwt = X(:, wtIdx);
    Xfad = X(:, fadIdx);

    [meanWT, varWT, nWT] = rowMeanVarNaN(Xwt);
    [meanFAD, varFAD, nFAD] = rowMeanVarNaN(Xfad);

    residualDF = nWT + nFAD - 2;

    ordinaryVariance = ...
        ((nWT - 1) .* varWT + (nFAD - 1) .* varFAD) ./ residualDF;

    moderatedVariance = ...
        (priorDF .* priorVariance + residualDF .* ordinaryVariance) ./ ...
        (priorDF + residualDF);

    standardError = sqrt(moderatedVariance .* (1 ./ nWT + 1 ./ nFAD));

    tStats = (meanFAD - meanWT) ./ standardError;

    degreesOfFreedom = residualDF + priorDF;

    pValues = twoSidedTPValue(tStats, degreesOfFreedom);

    invalid = ...
        ~isfinite(tStats) | ...
        ~isfinite(degreesOfFreedom) | ...
        degreesOfFreedom <= 0;

    pValues(invalid) = NaN;

    priorVarianceVector = repmat(priorVariance, size(pValues));
end


function [pValues, tStats, degreesOfFreedom] = welchTTestRows(group1, group2)

    [mean1, var1, n1] = rowMeanVarNaN(group1);
    [mean2, var2, n2] = rowMeanVarNaN(group2);

    standardError = sqrt(var1 ./ n1 + var2 ./ n2);

    tStats = (mean2 - mean1) ./ standardError;

    numerator = (var1 ./ n1 + var2 ./ n2) .^ 2;

    denominator = ...
        ((var1 ./ n1) .^ 2) ./ (n1 - 1) + ...
        ((var2 ./ n2) .^ 2) ./ (n2 - 1);

    degreesOfFreedom = numerator ./ denominator;

    pValues = twoSidedTPValue(tStats, degreesOfFreedom);

    invalid = ...
        ~isfinite(tStats) | ...
        ~isfinite(degreesOfFreedom) | ...
        degreesOfFreedom <= 0;

    pValues(invalid) = NaN;
end


function [pValues, numberOfPermutationsUsed] = permutationPValuesRows( ...
    X, ...
    wtIdx, ...
    fadIdx, ...
    useExhaustivePermutations, ...
    number_of_realizations)

    nSamples = size(X, 2);
    nWT = numel(wtIdx);

    observedDiff = mean(X(:, fadIdx), 2, "omitnan") - ...
                   mean(X(:, wtIdx), 2, "omitnan");

    if useExhaustivePermutations

        wtCombinations = nchoosek(1:nSamples, nWT);
        nPermutations = size(wtCombinations, 1);

    else

        nPermutations = number_of_realizations;
        wtCombinations = [];
    end

    permDiffs = zeros(size(X, 1), nPermutations);

    for p = 1:nPermutations

        if useExhaustivePermutations
            permWTidx = wtCombinations(p, :);
        else
            shuffledIdx = randperm(nSamples);
            permWTidx = shuffledIdx(1:nWT);
        end

        permFADidx = setdiff(1:nSamples, permWTidx);

        permDiffs(:, p) = ...
            mean(X(:, permFADidx), 2, "omitnan") - ...
            mean(X(:, permWTidx), 2, "omitnan");
    end

    nullMean = mean(permDiffs, 2, "omitnan");
    nullStd = std(permDiffs, 0, 2, "omitnan");

    valid = nullStd > 0 & isfinite(nullStd);

    zScores = nan(size(observedDiff));
    zScores(valid) = (observedDiff(valid) - nullMean(valid)) ./ nullStd(valid);

    pValues = nan(size(observedDiff));
    pValues(valid) = erfc(abs(zScores(valid)) ./ sqrt(2));

    numberOfPermutationsUsed = repmat(nPermutations, size(pValues));
end


function [rowMean, rowVar, rowN] = rowMeanVarNaN(X)

    X = double(X);

    finiteMask = isfinite(X);

    rowN = sum(finiteMask, 2);

    Xzero = X;
    Xzero(~finiteMask) = 0;

    rowMean = sum(Xzero, 2) ./ rowN;

    rowMean(rowN == 0) = NaN;

    centered = X - rowMean;
    centered(~finiteMask) = 0;

    rowVar = sum(centered .^ 2, 2) ./ (rowN - 1);

    rowVar(rowN < 2) = NaN;
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
        x = df ./ (df + tAbs .^ 2);
        p(valid) = betainc(x, df ./ 2, 0.5);
    end

    p(p < 0) = 0;
    p(p > 1) = 1;
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


function compactTable = makeCompactSignificantMoranTable(significantResultsTable)

    if height(significantResultsTable) == 0
        compactTable = significantResultsTable;
        return;
    end

    keepCols = { ...
        'gene', ...
        'region_name', ...
        'mean_WT_Moran_raw', ...
        'mean_5xFAD_Moran_raw', ...
        'mean_WT_Moran_normalized', ...
        'mean_5xFAD_Moran_normalized', ...
        'diff_5xFAD_minus_WT_normalized', ...
        'n_WT_valid', ...
        'n_5xFAD_valid', ...
        'pvalue', ...
        'FDR'};

    keepCols = keepCols(ismember(keepCols, significantResultsTable.Properties.VariableNames));

    compactTable = significantResultsTable(:, keepCols);

    compactTable = sortrows(compactTable, "FDR", "ascend");
end


function significantGenesPerRegion = summarizeSignificantMoranGenesPerRegion(significantResultsTable)

    if height(significantResultsTable) == 0

        significantGenesPerRegion = table( ...
            strings(0,1), ...
            zeros(0,1), ...
            zeros(0,1), ...
            zeros(0,1), ...
            nan(0,1), ...
            nan(0,1), ...
            'VariableNames', {'region_name', ...
                              'number_significant_genes', ...
                              'number_higher_in_5xFAD', ...
                              'number_higher_in_WT', ...
                              'minimum_pvalue', ...
                              'minimum_FDR'});
        return;
    end

    [G, regionList] = findgroups(significantResultsTable.region_name);

    significantGenesPerRegion = table();

    for i = 1:max(G)

        rows = G == i;

        regionTable = significantResultsTable(rows, :);

        numberSignificantGenes = numel(unique(regionTable.gene));

        numberHigherIn5xFAD = sum(regionTable.diff_5xFAD_minus_WT_normalized > 0);
        numberHigherInWT = sum(regionTable.diff_5xFAD_minus_WT_normalized < 0);

        minimumP = min(regionTable.pvalue, [], "omitnan");
        minimumFDR = min(regionTable.FDR, [], "omitnan");

        newRow = table( ...
            regionList(i), ...
            numberSignificantGenes, ...
            numberHigherIn5xFAD, ...
            numberHigherInWT, ...
            minimumP, ...
            minimumFDR, ...
            'VariableNames', {'region_name', ...
                              'number_significant_genes', ...
                              'number_higher_in_5xFAD', ...
                              'number_higher_in_WT', ...
                              'minimum_pvalue', ...
                              'minimum_FDR'});

        significantGenesPerRegion = [significantGenesPerRegion; newRow];
    end

    significantGenesPerRegion = sortrows(significantGenesPerRegion, ...
        "number_significant_genes", ...
        "descend");
end


function geneClean = standardizeGeneSymbol(geneRaw)

    geneClean = string(geneRaw);

    geneClean = replace(geneClean, char(160), " ");
    geneClean = regexprep(geneClean, '[\x00-\x1F\x7F]', '');
    geneClean = strtrim(geneClean);
    geneClean = regexprep(geneClean, '^["'']+|["'']+$', '');
    geneClean = regexprep(geneClean, '\s+', '');

    % Greek lowercase alpha -> Greek uppercase alpha
    geneClean = replace(geneClean, string(char(945)), string(char(913)));

    geneClean = upper(geneClean);

    %% Known aliases
    geneClean(geneClean == "IRG1") = "ACOD1";
    geneClean(geneClean == "INOS") = "NOS1";
    geneClean(geneClean == "CCRLG") = "CCRL2";
    geneClean(geneClean == "IDB3") = "ID3";

    %% HIF-1 alpha variants
    hifTarget = "HIF-1" + string(char(913));

    hifVariants = [
        "HIF-1?"
        "HIF-1A"
        "HIF1A"
        "HIF-1" + string(char(913))
        "HIF1" + string(char(913))
    ];

    geneClean(ismember(geneClean, hifVariants)) = hifTarget;
end


function excluded = isExcludedGene(gene)

    gene = string(gene);

    excludedGenes = ["NEG1", "NEG2", "NEG3", "FCGR3"];

    excluded = ismember(gene, excludedGenes);
end


function [regionClean, excludeRegion] = standardizeRegionName(regionRaw)

    regionClean = string(regionRaw);

    regionClean = replace(regionClean, char(160), " ");
    regionClean = regexprep(regionClean, '[\x00-\x1F\x7F]', '');
    regionClean = strtrim(regionClean);

    regionKey = lower(regionClean);
    regionKey = regexprep(regionKey, '\s+', '_');
    regionKey = replace(regionKey, "-", "_");

    excludeRegion = ...
        regionKey == "unassigned" | ...
        regionKey == "under_dg" | ...
        regionKey == "ca2";

    regionClean(regionKey == "sm") = "DG-CA1";
    regionClean(regionKey == "dg_ca1") = "DG-CA1";

    regionClean(regionKey == "inner_dg") = "Hilus";
    regionClean(regionKey == "hilus") = "Hilus";

    regionClean(regionKey == "upper_ca1") = "SO";
    regionClean(regionKey == "so") = "SO";

    regionClean(regionKey == "dg") = "DG";
    regionClean(regionKey == "ca1") = "CA1";
    regionClean(regionKey == "ca3") = "CA3";
    regionClean(regionKey == "slm") = "SLM";
end


function colName = findFirstColumn(varNames, candidates)

    varNames = string(varNames);
    candidates = string(candidates);

    normVars = normalizeName(varNames);
    normCandidates = normalizeName(candidates);

    colName = "";

    for i = 1:numel(normCandidates)

        idx = find(normVars == normCandidates(i), 1);

        if ~isempty(idx)
            colName = varNames(idx);
            return;
        end
    end
end


function namesNorm = normalizeName(names)

    namesNorm = lower(string(names));
    namesNorm = regexprep(namesNorm, '[^a-z0-9]', '');
end


function x = toDoubleVector(v)

    if isnumeric(v)
        x = double(v);
    else
        x = str2double(string(v));
    end

    x = x(:);
end
