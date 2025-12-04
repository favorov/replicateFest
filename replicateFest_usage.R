# 2025-01-24
# Ludmila Danilova, Leslie Cope
#================================================================
# usage of replicateFest on mock data with replicates
#================================================================
library(tools)
library(replicateFest)
#=============
# run test data with replicates
inputDir = "./tests/testthat/testdata/with_replicates/"

# list paths to files with data
files = list.files(inputDir, full.names = T,
                   pattern = "tsv", recursive = TRUE)

filenames = file_path_sans_ext(basename(files))
sampAnnot = splitFileName(filenames)

# run all clones in a patient and time point and return the results
res = runExperiment(files,
                    peptides = sampAnnot$condition,
                    refSamp= "control",
                    fdrThr = 0.05,
                    nReads = 20,
                    percentThr = 0,
                    xrCond = NULL,
                    ntLevel = F,
#                   outputFile = "mock-data-replicates-output.xlsx",
                    saveToFile = F)

# save results to be used during package tests in testthat
#saveRDS(res, file = paste0(inputDir, "replicate_results.rds"))


