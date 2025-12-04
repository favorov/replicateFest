# 2025-12-01
# Ludmila Danilova
#================================================================
# run replicateFest on the previously published HIV data
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2020.00591/full
#================================================================
library(replicateFest)
#=============
# load previously saved RDA with all data (saved through shiny app)
load("Blankson_inputData.rda")

# extract conditions
sampAnnot = splitFileName(names(mergedData))

orThr = 1
fdrThr = .05
percentThr = 0
refSamp = "NoPep"

# get clones to test
clonesToTest = getClonesToTest(mergedData, nReads = 50)
# run the analysis for selected clones
res = fitModelSet(clonesToTest,
                  mergedData,
                  peptides = sampAnnot$condition,
                  excludeCond = NA,
                  refSamp)
rownames(res) = res$clone

posClones = getPositiveClonesReplicates(res,
                                        mergedData,
                                        refSamp,
                                        samp = sampForAnalysis,
                                        orThr,
                                        fdrThr)

resTable = createResTableReplicates(res,
                                    mergedData,
                                    refSamp,
                                    orThr,
                                    fdrThr)
resToExcel = list()
resToExcel$expandedTable = replicateFest:::getExpanded(res,mergedData, 
                  refSamp,
                  orThr,
                  fdrThr)

saveResults(resToExcel, "text.xlsx")

