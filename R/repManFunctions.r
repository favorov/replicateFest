
# Functions for analysis of FEST data in shiny app
# using input with replicates
# based on Leslie's development for the Cervical SPORE project

#' readMergeSave
#' Reads files, removes non-productive sequences, extracts counts,
#' creates all necessary objects for further analysis
#'
#' @importFrom dplyr bind_rows filter select group_by summarise %>%
#' @importFrom multcomp glht
#' @importFrom stats p.adjust setNames
#' @importFrom pheatmap pheatmap

#' @export
#' @param files a list of filenames with full paths
#' @param filenames a vector of filenames
#' @return a list of AA level counts, nucleotide level counts, and total number of reads
#
readMergeSave = function(files, filenames = NULL)
{
	if (is.null(files))
	{
		print('There is no files to read')
			return(NULL);
	}
  # TODO add check if files are exist

		require(tools)
    require(immunarch)

  # output objects
		mergedData = ntData = list()

		# read all files with immunarch functionality
		repertoire = repLoad(.path = files)

		# the names of files that were actually read in
		readFiles = c()
		for(i in names(repertoire$data))
		{
#			print(i)
		  dat = repertoire$data[[i]]
				#count reads of productive sequences only
				mergedData[[i]] = tapply(dat$Clones, dat$CDR3.aa, sum, na.rm = T)
				# nucleotide level data
				ntData[[i]] =tapply(dat$Clones, dat$CDR3.aa, sum, na.rm = T)
				readFiles = c(readFiles, i)
		}
		if (length(mergedData) == 0)
		{
			print(paste('There are no data to read'))
			return(NULL)
		}
		# if file names are not supplied, use internal shiny server file names (0,1,2,..)
		if (is.null(filenames))
		{
		  filenames = readFiles
		} else {
		  # it files names are supplied (actual file names that were loaded)
		  # take file names that were actually read and create an index
		  # shiny server saves files with 0,1,2,.. names
		  ind = as.numeric(names(repertoire$data))+1
		  filenames = sapply(unlist(filenames[ind]),file_path_sans_ext)
		}
		# assign file names as names to objects
#		browser()
		# if not all loaded files
		names(mergedData) = names(ntData) = filenames
		return(list(mergedData = mergedData,ntData = ntData))
}

#' cTabPR
#' function to make the count matrix for 1 clone, from merged data
#' @param clone a clone to get counts for
#' @param countData a list of per clone counts for all samples from readMergeSave
#' @param correct a parameter to add to all counts to avoid 0s
#' @return a vector with counts for a clone across all samples
#'
#'
#### requires 1 input:
###### 1) countData=merged data from readMergeSave
#### correct is a parameter to add to all counts to avoid 0s
cTabPR=function(clone,countData,correct=.5){
    # replace missing values with 0
     minna=function(x){  ### function to deal with missing values
        if(is.na(x)) x=0
        return(x)}
    # get counts for a clone across all samples
     cts=sapply(countData,function(x) return(x[clone]))
     # replace missing values with 0
     cts=sapply(cts,minna)
     # get total read count for each sample minus the count for the clone
     sms=sapply(countData,sum)-cts
     ans=cbind(cts,sms)+correct
     return(ans)
}

#### function to make the data frame for regression
#### requires 3 inputs:
###### 1) cTab = counts object created by cTabPR,
###### a matrix with two columns: counts of a clone and sums minus counts of the clone in all samples
###### 2) peps = vector of peptides corresponding to columns in merged data
###### 3) refSamp=name of refSamp peptide
cDfPR=function(cTab,peps,refSamp="NPA"){
  # if(!(refSamp%in%peps)){ break("refSamp peptide not found")}else{
  datPR=data.frame(as.vector(cTab), ## convert matrix to vector: counts of a clone (c+) come first, then sums minus counts of the clone (c-)
                                    ## then sums minus counts of the clone (c-)
                   factor(rep(c("c+","c-"),rep(nrow(cTab),2)),levels=c("c-","c+")),
                   factor(rep(peps,2),
                          levels=c(refSamp,sort(setdiff(peps,refSamp)))))
  colnames(datPR)=c("cts","clone","pep")
  return(datPR)
}

#### function to convert regression output to OR scale and format
#### requires 1 input:
###### 1) dt = row from mhcMod$coef matrix
######

prepStatsPR=function(dt){  ## dt=row from coefficent matrix

    names(dt)=c("Est","SE","Z","P")
    OR=round(exp(dt["Est"]),3)
    LCB=round(exp(dt["Est"]-1.96*dt["SE"]),3)
    UCB=round(exp(dt["Est"]+1.96*dt["SE"]),3)
    pval=dt["P"]
    ans=c(OR,LCB,UCB,pval)
    names(ans)=c("OR","LCB","UCB","pval")
    return(ans)
}

#### function to organize counts and regression results into 1 data.frame for export
#### messy, surely can be improved
#### requires 2 inputs:
###### 1) cts = counts object created by cTabPR
###### 2) cfs = matrix extracted from model results as
###### coefs[interact,] where
##### coefs=summary(mhcMod)$coef; interact=grep(":",rownames(coefs),fixed=T)
##### probably better to either operate straight on mhcMod or add a function like cTabPR to do so
dfResultPR=function(cts,cfs){
    #rownames(cfs)=gsub(":","|",rownames(cfs),fixed=T)
    ctsM=rbind(c("","counts","sums","",""),cbind(rownames(cts),cts,matrix(rep("",nrow(cts)*2),ncol=2)))
    cfsM=rbind(c("",colnames(cfs)),cbind(rownames(cfs),cfs))
    data=data.frame(rbind(ctsM,matrix(rep("",10),nrow=2),cfsM))
    return(data)
}


#####################################
#' fitModel
#' Perform analysis for one clone
#' @description  The function fits negative binomial regression model for a clone
#' @export
#' @param clone a clone to get counts for
#' @param countData a list of per clone counts for all samples from readMergeSave
#' @param peptides a vector of peptides corresponding to columns in merged data
#' @param refSamp name of reference condition
#' @param c.corr a parameter to add to all counts to avoid 0s
#' @return odds ratios and p-values of regression for comparisons to the reference condition
#' and comparison of the best to the second best condition
#'


fitModel = function(clone,countData,peptides,refSamp,
                    c.corr=1){

  #### peptides is a vector, equal in length to merged data indicating
  ####  which peptide is represented in each rep.
  #### could extend with a matrix of covariates, rows = length merged data

  #  print(clone)
  ### make the count matrix for a clone across all samples
  ctsPR=cTabPR(clone,countData,correct=c.corr)
  ### make the regression data
  datPR=cDfPR(ctsPR,peps=peptides, refSamp=refSamp)

  ### perform the regression
  #mhcMod <- glm(cts ~ clone*pep, data = datPR, family = poisson(link = "log"))

  # run negative binomial regression
  # if the model fails, return NULL
  tryCatch({
    mhcMod <- glm.nb(cts ~ clone*pep, data = datPR)
  }, error = function(e) {
    print(clone)
    print(e)
    return(NULL)
  })
  # if the model doesn't fail, extract coefficients
  if(!exists("mhcMod")){return(NULL)}

  # if the model fits fine, extract coefficients
  coefs=summary(mhcMod)$coef
  # find interactions to include in the output
  interact=grep(":",rownames(coefs),fixed=T)
  # subset coefficients for interactions that represents results of comparison to refSamp
  coefs_int = coefs[interact,]
  # add condition from which coefficients were extracted
  # and convert coefficients to OR
  condition = gsub("clonec+:pep","",rownames(coefs_int), fixed = T)
  OR = setNames(round(exp(coefs_int[,"Estimate"]),3),
                paste0("OR: ", condition, "_vs_",refSamp))
  pval = setNames(coefs_int[,"Pr(>|z|)"],
                paste0("pval: ", condition, "_vs_",refSamp))

#browser()
  # compare with the second best clone and output OR and p-value
  ### order conditions and find the best and second best interaction coefficients
  pepOrd=coefs_int[order(coefs_int[, "z value"],decreasing=T),]
  # names of the best and the second best conditions
  best=rownames(pepOrd)[1]
  scnd=rownames(pepOrd)[2]

  ### make a contrast matrix “cMat” with columns matching model coefficients,
  # and 1 row
  ### initialize with 0 for all irrelevant conditions,
  # and use the value 1 for the most expanded condition,
  # and -1 for the second most expanded condition
  cMat <- matrix(rep(0,length(mhcMod$coef)), 1) ## initialize contrast matrix
  # set colnames as names of coefficients
  colnames(cMat) = rownames(coefs)
  ### specify contrast top peptide to second
  cMat[1,c(best,scnd)]=c(1,-1)

  ### fit the contrast using glht function,
  # requires the original model,
  # “mhcMod” and the new contrast matrix
  # get an estimate
  pepComp<- glht(mhcMod, linfct=cMat)  ### fit

  ### extract statistics from the contrast result,
  # in this instance just p-value,
  # but OR between those conditions
  bestCond = gsub("clonec+:pep","",best, fixed = T)
  scndCond = gsub("clonec+:pep","",scnd, fixed = T)
  # convert an estimate to OR
  OR_scnd = setNames(round(exp(summary(pepComp)$test$coeff),3),
                "OR: best_vs_second")
  # p-value
  pval_scnd = setNames(summary(pepComp)$test$pvalues,
                  "pval: best_vs_second")
  # names of the best and the second best conditions
  best_vs_second = setNames(paste0(bestCond,"_vs_",scndCond), "second_comparison")

  # combine output and return it
  res = c(clone = clone,OR,pval, OR_scnd, pval_scnd, best_vs_second)

  return(res)

}


#===========
# wrapper for running the full analysis using negative binomial
# from reading files to output all results
#################
#' @export
#' @title runExperiment
#' @description Reads in files with TCR repertoires from
#' a FEST experiment with replicate samples per condition (stimulating peptide).
#' It fits negative binomial model to find expanded clones
#' comparing to a reference samples.
#' It also compares top conditions to find unique expansions.
#' The results are return and saved in an Excel file.
#' @param files a list of filenames with full paths
#' @param peptides a vector of peptides corresponding to columns in merged data
#' @param nReads minimal number of reads required to consider a clone
#' @param refSamp name of reference condition
#' @param orThr threshold for OR to consider a clone expanded
#' @param fdrThr threshold for FDR to consider a clone expanded
#' @param excludeCond a vector of conditions to exclude from the analysis
#' @param xrCond a vector of cross-reactive conditions
#' @param percentThr a threshold for percentage of reads in a sample to consider a clone expanded
#' @param outputFile name of the output file
#' @param saveToFile logical, if TRUE save results to a file
#' @param permute logical, if TRUE permute sample labels to run a permutation test
#' @return a list of all expanded clones, uniquely expanded clones
#' and parameters of the run
#'
#### returns a list of all expanded clones, uniquely expanded clones
#### and parameters of the run

#### requires the following inputs:
#### main arguments are filenames (full paths in current implementation)
#### and a vector of peptide ids for each file.
#### additional arguments include a minimal number of reads required
#### to consider a clone
#### the reference peptide ID,
#### a vector of conditions to exclude from the analysis
#### a threshold for OR and FDR to consider a clone expanded
#### a vector of cross-reactive conditions
#### added option for permuting labels to answer to Kellie's request. Need to remove for the app and package

#### 2025-01-29 added comparison to the second best to find unique expansions

runExperiment=function(files,
                       peptides,
                       nReads=50,
                       refSamp,
                       orThr=1,
                       fdrThr = 0.05,
                       excludeCond = NA,
                       xrCond = NULL,
                       percentThr = 0,
                       ntLevel = FALSE,
                       outputFile = "output.xlsx",
                       saveToFile = T, permute = FALSE)
{

  #### start algorithm read data
  #### start algorithm read data
  inputData = readMergeSave(files, filenames = NULL)
  # extract aa or nt level data for the downstream analysis
  if (ntLevel) mergedData = inputData$ntData else mergedData = inputData$mergedData

  # permute sample labels in mergedData to run permutation test
  if (permute){
    set.seed(123456)
    # create sampling
    s = sample(1:length(mergedData),size = length(mergedData), replace = F)
    # update sample names
    names(mergedData) = names(mergedData)[s]
    # update peptides
    peptides = peptides[s]
  }

  # get clones to test
  goodClones = getClonesToTest(mergedData, nReads = nReads)
  print(c("good clones #",length(goodClones)))

  if(length(goodClones) == 0)
  {
    print("There are not clones to analyze. Try to reduce the number of tempaltes.")
    return(NULL)
  }
#browser()
  # run the analysis for selected clones
  fitResults = fitModelSet(goodClones,
                           mergedData,
                           peptides,
                           excludeCond = excludeCond,
                           refSamp=refSamp,c.corr=1)
  rownames(fitResults) = fitResults$clone


  # get positive (uniquely expanded) clones
  posClones = getPositiveClonesReplicates(fitResults,
                                          mergedData,
                                          refSamp = refSamp,
                                          excludeCond = excludeCond,
                                          orThr = orThr,
                                          fdrThr = fdrThr,
                                          percentThr = percentThr)
  tablesToXls = list()
  # if there is no positive clones
  if (nrow(posClones)==0)
  {
    print('There are no positive clones. Try to adjust thresholds')
    tablesToXls$summary = data.frame('There are no positive clones', row.names = NULL, check.names = F)
  }else{
    # if there are positive clones, save them in to Excel file
  # create table with results
#    browser()
  tablesToXls = createPosClonesOutput(posClones,
                                      mergedData,
                                      refSamp,
                                      replicates = TRUE)

  }

  resTable = createResTableReplicates(fitResults,
                                      mergedData,
                                      percentThr,
                                      refSamp,
                                      orThr,
                                      fdrThr)
  # if there are clones
  if(nrow(resTable) > 0)
  {
    tablesToXls$ref_comparison_only = resTable
  }else{
    tablesToXls$ref_comparison_only =
      data.frame(res = 'There are no significant clones')
  }

  # add a spreadsheet with cross-reactive clones if specified
  if(!is.null(xrCond))
  {
    tablesToXls$cross_reactive = getXR(fitResults,
                                       peptides,
                                       refSamp,
                                       xrCond = xrCond,
                                       excludeCond = excludeCond,
                                       percentThr = percentThr,
                                       countData = mergedData,
                                       orThr = orThr,
                                       fdrThr = fdrThr)

  }


  # save parameters of analysis
  s = names(mergedData)
  productiveReadCounts = sapply(mergedData, sum)
  param = c("Data with replicates",
            'Reference condition',
            'Excluded conditions',
            'Cross-reactive conditions',
            'n template threshold','FDR threshold',
            'OR threshold','percent threshold',
            'Nucleotide level analysis',
            'n samples',
            paste(s, 'n templates',sep = '_'))
  value = c(TRUE,
            refSamp,
            paste(excludeCond, collapse = ', '),
            paste(xrCond, collapse = ', '),
            nReads,
            fdrThr,
            orThr,
            percentThr,
            ntLevel,
            length(s), productiveReadCounts[s])

  tablesToXls$parameters = data.frame(param, value)

  # save into Excel file
  if(saveToFile)
  {
    saveResults(tablesToXls, outputFile = outputFile)
  }
  return(tablesToXls)
}

#' @export
#' @title getAbundances
#' @description returns the read count for clones of interest in all samples
#' @param clones a vector of clones of interest
#' @param countData a list of counts for all samples
#' @return a data.frame of clone abundances across samples in `countData`
# input: a list of merged data, a vector of clones of interest
getAbundances = function(clones,countData)
{
  # create output matrices
  output_counts = data.frame(matrix(0,nrow = length(clones),
                                    ncol = length(countData)))
  rownames(output_counts) = clones
  colnames(output_counts) = names(countData)

  # get the read count for clones in each sample
  for (i in names(countData))
  {
    rows = intersect(names(countData[[i]]),rownames(output_counts))
    output_counts[rows,i] = countData[[i]][rows]
  }
  # update colnames to add "abundance"
  colnames(output_counts) = paste(names(countData),'abundance', sep = '_')
  return(output_counts)
}

#' @title getClonesToTest
#' @description function that returns clones of interest in all samples
#' @param countDat a list of merged data
#' @param nReads a minimal number of reads required to consider a clone
#' @return a vector of clones of interest
#' @export
#'
getClonesToTest = function(countDat, nReads = 50)
{
  # all clones
  clones=unlist(sapply(countDat,function(x)  return(names(x))))
  # the corresponding counts
  cts=unlist(sapply(countDat,function(x)  return(x)))
  #
  maxCt=tapply(cts,clones,max)

  # clones that have more than nReads reads to run the analysis
  goodClones=names(maxCt)[which(maxCt>nReads)]

  return(goodClones)
}

#' fit model for a set of clones
#' @param clones a list of clones to fit the model
#' @param countData a list of counts for all samples
#' @param peptides a vector of peptides corresponding to columns in merged data
#' @param excludeCond a vector of conditions to exclude from the analysis
#' @param ... additional parameters to pass to fitModel
#' @return a matrix with all ORs, p-values and FDRs
#' @export
# return a matrix with all ORs, p-values and FDRs
# input: a list of clones, merged data, peptides,

fitModelSet = function(clones, countData, peptides,
                       excludeCond = NA, ...)
{
    # run model for "good" clones
  if (length(excludeCond)>0)
  {
    # get indexes of conditions to include in the analysis
    incl = which(!(peptides %in% excludeCond))
    cat("Conditions to include:", peptides[incl],"\n")
    fitResults=lapply(clones,fitModel,countData=countData[incl],
                  peptides=peptides[incl],...)
  }else{
    fitResults=lapply(clones,fitModel,countData=countData,
                  peptides=peptides,...)
  }

  #browser()

  # remove enties with NULL
  fitResults = fitResults[!sapply(fitResults, is.null)]
  # convert list to a data.frame
  fitResults = as.data.frame(bind_rows(fitResults))

  # add FDR adjustment
  # find columns with p-values
  pvalCol = grep("pval",colnames(fitResults), value = T)
  fdrs = c()
  for(i in pvalCol)
  {
    fdrs = cbind(fdrs,p.adjust(as.numeric(fitResults[,i]), method = "BH"))
  }
  # add colnames for fdrs and add to the fitResults matrix
  colnames(fdrs) = paste0("FDR:",gsub("pval:","",pvalCol))
  fitResults = cbind(fitResults, fdrs)


  return(fitResults)
}

#' function that returns the expanded clones relative to reference
#' @param fitResults a data frame with ORs, p-values and FDRs
#' @param countData a list of counts for all samples
#' @param refSamp name of reference sample
#' @param orThr threshold for OR to consider a clone expanded
#' @param fdrThr threshold for FDR to consider a clone expanded
#' @return a data frame with expanded clones
#' @export
# input: a data frame with ORs, p-values and FDRs
# output: a data frame with expanded clones
getExpanded = function(fitResults, countData,
                       refSamp,
                       orThr = 1,
                       fdrThr = 0.05)
{
  # find colunms related to comparison to reference condition
  contCol = grep(refSamp,colnames(fitResults), value = T)
  # subset for those columns only
  fitResults = fitResults[,c("clone",contCol)]
  # find all significantly expanded clones
  # find comparison names
  comp = grep("OR",colnames(fitResults), value = T)
  comp = gsub("OR: ","",comp)

  # convert OR and FDR columns of fitResults into numeric values
  orCols = grep("OR:", colnames(fitResults), value = T)
  fdrCols = grep("FDR:", colnames(fitResults), value = T)

#   # remove columns corresponding to best_vs_second
#   orCols = orCols[!grepl("best_vs_second", orCols)]
#   fdrCols = fdrCols[!grepl("best_vs_second", fdrCols)]
#   comp = comp[!grepl("best_vs_second", comp)]

  # convert those columns to numeric values
  fitResults[,orCols] = sapply(fitResults[,orCols],
                                            as.numeric)
  fitResults[,fdrCols] = sapply(fitResults[,fdrCols],
                               as.numeric)
  # find expanded clones
  expandedClones = c()
  for (i in comp){
    #print(i)
    expandedClones = union(expandedClones,
                           rownames(fitResults)[which(as.numeric(fitResults[,paste0("OR: ",i)]) >= orThr & # expanded
                                                    as.numeric(fitResults[,paste0("FDR: ",i)]) < fdrThr)])# significant
  }

  # get the results for expanded clones only
  res_exp = fitResults[expandedClones,]


  # #
  # # # a table with significant clones only
  # res_exp = fitResults[which(fitResults[,orCols]>= orThr &
  #                              fitResults[,fdrCols]< fdrThr),]

  # add columns that indicates how many and what comparisons were significant
  # get T/F matrix for significant comparisons

 # browser()
  sig = (res_exp[,fdrCols] < fdrThr &
           res_exp[,orCols] >= orThr)
  # list significant comparisons
  sigComp = apply(sig, 1, function(x){
    # select significant comparisons and get first condition before "vs"
    s = sapply(strsplit(comp[x], split = "_vs_"), getElement, 1)
    # list significant comparisons using comma
    paste(s, collapse = ",")
    })
  # add the number and the list to the results
  res_exp = cbind(clone = res_exp[,"clone"],
                  n_significant_comparisons = rowSums(sig, na.rm = T),
                  significant_condition = sigComp,
                  res_exp[,setdiff(colnames(res_exp),c("clone"))])
   # add abundance and percentage of the top clones in each condition

  #================
  # TODO replace with getCountsPercent
  # get abundance for the top clones
  abundance = getAbundances(rownames(res_exp), countData)
  # get total read count for each sample
  totalReadCountPerSample = sapply(countData, sum)
  # calculate the percentage of each clone in each sample
  percentage = round(sweep(abundance, 2, totalReadCountPerSample, "/")*100,3)
  colnames(percentage) = paste(names(countData),'percent', sep = '_')

  res_exp = cbind(res_exp,
                  abundance[rownames(res_exp),],
                  percentage[rownames(res_exp),])

  return(res_exp)

}

# function that returns the cross-reactive clones
# input: a data frame with expanded clones
# and a vector of cross-reactive conditions
# output: a data frame with cross-reactive clones

getXR = function(res, conditions, refSamp, xrCond,
                 excludeCond = NULL,percentThr = 0, ...)
{
  # get all expanded clones relative to the refSamp
   res_exp = getExpanded(res,refSamp,
                        ...)

  #=================
  # keep clones with maximum percentage across
  # all analyzed samples higher that a specified threshold
  #=================
  # grep columns with percentage
  percCol = grep("percent",colnames(res_exp), value = T)
  #browser()
  # exclude columns with excludeCond
  if(!is.null(excludeCond))
    percCol = percCol[!grepl(paste(excludeCond,collapse = "|"),
                             percCol)]
  # get clones with maximum percentage higher than specified threshold
  res_exp_filtered = res_exp[apply(res_exp[,percCol],1,max) > percentThr,]

  if(nrow(res_exp_filtered) == 0) return(NULL)

    # a vector of conditions that shouldn't be cross-reactive
  # find all and take difference
  allCond = res_exp_filtered %>% filter(n_significant_comparisons == 1) %>%
    dplyr::select(significant_condition) %>% unlist() %>% unique()

  # conditions that shouldn't be cross-reactive
  excludeCond = setdiff(conditions, xrCond)

  # find clones that are reactive in more than 1 condition (cross-reactive clones)
  res_xr = res_exp_filtered %>% filter(n_significant_comparisons > 1)
  # find clone that are cross-reactive for conditions in xrCond
  inclXR = res_xr %>% filter(grepl(paste(xrCond,collapse = "|"),significant_condition))

  # exclude conditions that are not in xrCond
  excl = inclXR %>% filter(!grepl(paste(excludeCond,collapse = "|"),significant_condition))
  return(excl)
}

#' saveResults
#' @export
#' @title saveResults
#' @description Saves results to an excel file
#' @param results a list of data frames with results
#' @param outputFile name of the output file

saveResults = function(results, outputFile = "output.xlsx")
{
#  library(openxlsx)
  # create a workbook
  wb = createWorkbook()

  for(i in names(results))
  {
    # add sheets
    addWorksheet(wb, i)
    # add data to sheets
    writeData(wb, i, results[[i]])
  }
  # save the workbook
  saveWorkbook(wb, outputFile, overwrite = TRUE)
}

#' splitFileName
#' function that extracts condition and replicate information
#' from the file names. The condition and replicate should
#' be separated by "_" and be the last two elements.
#' @param filenames a vector of file names
#' @param sep a separator to split the file names, default is "_"
#' @return a data frame with condition and replicate information
#' @export

splitFileName = function(filenames, sep = "_")
{
  # split the file names by "_"
  splitNames = strsplit(filenames, split = sep)
  # create a data frame with condition and replicate information
  condRep = data.frame(matrix(nrow = length(filenames), ncol = 2,
                       dimnames = list(filenames,c("condition","replicate"))))
  # get the last two elements to fill in the data frame
  for( i in 1:(length(filenames)))
  {
    l = length(splitNames[[i]])
    if(l>1) condRep[i,] = splitNames[[i]][(l-1):l]
  }
  # add file name
  condRep = cbind(file = filenames, condRep)
  rownames(condRep) = NULL
  return(condRep)
}

#' getPositiveClonesReplicates
#' returns positive clones for version with replicates
#' @param analysisRes a data frame with results of analysis
#' @param mergedData a list of data frames with read counts for each sample
#' @param refSamp a name of reference condition to compare to
#' @param samp a vector of sample names to include in the analysis
#' @param excludeCond a vector of conditions to exclude from the analysis
#' @param orThr a threshold for odds ratio
#' @param fdrThr a threshold for FDR
#' @param percentThr a threshold for percentage of reads in a sample to consider a clone expanded
#' @return a data frame with positive clones, significant condition,
#' @export

getPositiveClonesReplicates = function(analysisRes,
                                       mergedData,
                                       refSamp,
                                       samp = names(mergedData),
                                       excludeCond = NA,
                                       orThr = 1,
                                       fdrThr = 0.05,
                                       percentThr = 0)
{
  # get all expanded clones relative to the refSamp
  res_exp = getExpanded(analysisRes, mergedData, refSamp,
                        orThr = orThr, fdrThr = fdrThr)
#browser()
  #=================
  # keep clones with maximum percentage across
  # all analyzed samples higher that a specified threshold
  #=================
  # grep columns with percentage
  percCol = grep("percent",colnames(res_exp), value = T)
  #browser()
  # exclude columns with excludeCond
  if(!is.null(excludeCond)) percCol = percCol[!grepl(paste(excludeCond,collapse = "|"),percCol)]
  # get clones with maximum percentage higher than specified threshold
  res_exp_filtered = res_exp[apply(res_exp[,percCol],1,max) > percentThr,]

  #======
  # find uniquely expanded clones by checking the second best clone
  # save the second best comparison results
  # these are the rest of the columns that are not comparison to reference
  screen_scndBest = analysisRes[rownames(res_exp_filtered),
                                grepl("second",colnames(analysisRes))]
  # check for uniqueness. it should be expanded and significant in comparison to the second best as well
  unique_exp = (screen_scndBest[,grep("OR", colnames(screen_scndBest))]>1 &
                  screen_scndBest[,grep("FDR", colnames(screen_scndBest))]<fdrThr)
  # merge the results
  # add columns that indicates how many and what comparisons were significant
  # results for the second best comparison
  # the rest of info
  sigComCol = c("clone","n_significant_comparisons")
  res = cbind(res_exp_filtered[unique_exp,sigComCol],
              significant_comparison = res_exp_filtered[unique_exp,"significant_condition"],
              res_exp_filtered[unique_exp,setdiff(colnames(res_exp),sigComCol)],
              screen_scndBest[unique_exp,])

  # extract conditions from the file names
  sampAnnot = splitFileName(names(mergedData))
  # extend output by replacing significant_condition with sample names,
  # so for every clone, there will several replicates instead of one condition
  resExt = merge(res[,c("clone","significant_comparison")],
                 sampAnnot[,c("file","condition")],
                 by.x = "significant_comparison",
                 by.y = "condition")
  # fix colnames to match input for createPosClonesOutput
  colnames(resExt) = c("significant_comparison","clone","significant_condition")

  #=====
  return(resExt[,c("clone","significant_condition","significant_comparison")])
}

#' createResTableReplicates
#' @description Creates a table with significantly expanded clones
#' relative to the reference condition and all corresponding data,
#' such as OR, FDR, abundances, and percentages
#' @param res a table with results of analysis
#' @param mergedData a list of data frames with read counts for each sample
#' @param percentThr a threshold for percentage of reads in a sample to consider a clone expanded
#' @param ... additional parameters to pass to getExpanded function
#' @return a data frame with significant clones and the corresponding
#'  OR, FDR, counts, and percentages for all conditions
#' @export
# write output
# OR, p-value, FDR, abundance, percent
createResTableReplicates = function(res,mergedData,
                          percentThr = 0,
                          ...)
{

  # get all expanded clones
  #browser()
  tab = getExpanded(res,mergedData, ...)

  # grep columns with percentage
  percCol = grep("percent",colnames(tab), value = T)
  # get clones with maximum percentage higher than specified threshold
  tab = tab[apply(tab[,percCol],1,max) > percentThr,]


  # add check if there are any rows in tab
  if (nrow(tab) == 0)
  {
    m = 'There is no significant clones after applying percent and condition thresholds'
    print(m)
    return (NULL)
  }


  return(tab)
}
