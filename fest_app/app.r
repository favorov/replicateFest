# 2017-11-14 (v6)
# add fields for the number of cells and probability to calculate baseline frequency
# to filter out clones with low frequency
# this version uses only clones that appear in baseline above calculated frequency
# previously, we used 10 reads as a threshold and do not restrict to baseline clones

# 2018-01-23 (v7)
# Added a threshold for the number of reads into interface
# Allowed user to ignore baseline threshold
# Make only one heatmap with all positive clones
# Output all significant clones relative to the reference even there are no positive clones

# 2018-03-07 (v7.1)
# Trying to solve out of memory issue when saving large table to Excel
# Solved by using WriteXLS and saving all results at once
# fixed problem with using AA level data when selecting clones to analyze with Use nucleotide level option selected

# 2018-04-18 (v8)
# Added a tab with user's manual to UI
# Added baseline percent and the corresponding number of templates in baseline to output parameter sheet

# 2018-06-07 (v9)
# switch to the parser from the tcR package to eliminate sequencing error on 5' end of the nt sequence

# 2018-07-20 (v10)
# fixed reading function. It doesn't save mergedData and ntData to the disc;
# it returns a list with two elements - mergedData and ntData
# fixed creating Excel file with no significant clones

# 2019-06-04
# fixed output parameters issue
# add the OR threshold when compare to all other conditions (previously, it was only the FDR threshold)

# 2019-06-14 (v12)
# improved an input process. Now a user can select what file format to upload: Adaptive, VDJtools, or the previously saved R object
# add excluded samples to output
# remove nRead requirement when compare with other wells

# 2020-05-03 (v13)
# switch to immunarch: can't install on shiny server
# make comparison with reference as an option (checkbox)
# if off, then each condition will be compared with top two other conditions

# 2023-12-01
# switched to immunarch

# 2025-02-05 (v14)
# Split loading the data and analysis to the different tabs
# added percentage threshold to the Fisher's test version
# fixed output to heatmap when there are no or one clone
# added all functions to the replicateFest package

# 2025-05-06 (v15)
# - made load data, run analysis, and save the results to the different tabs
# - added analysis with replicates
# - increased the default for the number of reads
# - removed the baseline option
# - switched to openxlsx for saving results
# - made a list with analysis results and parameters that were used for the analysis
# - made custom separator for file names in the splitFileName function to separate replicates
# - switched to pheatmap
# - added filter for the percent of non-zero wells
# - added threshold for frequency

# TODO:
# - add custom separator into interface
#   Couldn't easily add, because the getPerSampleSummary function
#   also requires separator, but it's not straightforward to pass it
#   Probably, need to save that as another element in analysisRes object
#  added a custom separatro to the splitFileName function for now

# - what do you think about removing the “compare to reference* button
# in favor of  testing reference==’none’?
# (Maybe want to reinitialize when new data is loaded
# to avoid errors)

# - A message saying “Analysis in progress” would help,
# to be replaced with “Analysis is done” , maybe also a
# note about expected wait time

#==============================
# User interface
#==============================
ui <- fluidPage(
  #  title = "FEST data analysis",
  headerPanel("FEST data analysis"),
  # add description of the app and steps
  tags$div(
    tags$p("The FEST (Functional Expansion of Specific T-cells)
           application for analysis of TCR sequencing data of short-term,
           peptide-stimulated cultures to identify antigen-specific
           clonotypic amplifications. The app loads TCR sequencing data
           in the tab-delimited format, performs the analysis, and
           visualizes and saves results. For analysis, we use only
           productive clones and summarize template counts for
           nucleotide sequences that translated into the same amino acid
           sequence."),
    tags$h3("Steps"),
    tags$ol(
      tags$li(HTML("<b>Load Data</b>: upload FEST files or a previously saved R
              object with data")),
      tags$li(HTML("After the data is successfully loaded,
              select the <b>Run Analysis</b>
                   tab to run the analysis.
                   <br><b>Important!</b> Please make sure to specify
                   whether the input data has replicates or not.
                   If the input data has replicates, the conditions
                   will be extracted from the input file names.
                   To be able to correctly extract the conditions, file names should
                   follow the format 'sampleID_condition_replicate.ext'.
                   E.g. 'sample1_HIV_1.tsv'.
                   <br>If the input data has replicates,
                   the negative binomial model will be used in the analysis;
                   if the input data does not have replicates,
                   the Fisher's exact test will be used.")),
      tags$li(HTML("After the analysis is done,
              select the <b>Save Results</b> tab
                   to save the results.")),
    ),
    tags$p(HTML("For more details, please refer to the
                <b>User Manual</b> tab"))
  ),

  # layout with tabs
  tabsetPanel(
    #============================
    # tab for loading the data
    tabPanel("Load Data",
             sidebarLayout(
               # the left side panel
               sidebarPanel(
                 # select the format of input files
                 selectInput('inputFiles', 'Select an input format',
                             choices = c('FEST files','an R object with data'),
                             selected = NULL, multiple = FALSE,selectize = TRUE,
                             width = NULL, size = NULL),

                 # if a previously saved R object with the input data will be uploaded
                 conditionalPanel(
                   condition = "input.inputFiles == 'an R object with data'",
                   fileInput('inputObj', 'Upload an R object with data (.rda)',
                             multiple = FALSE,
                             accept=c('rda', '.rda'))
                 ),
                 # show the fileInput control depending on the selected values of inputFiles
                 # if files with raw data will be uploaded
                 conditionalPanel(
                   condition = "input.inputFiles == 'FEST files'",
                   fileInput('sourceFiles', 'Upload FEST files',
                             multiple = TRUE,
                             accept=c('text/tsv', 'text/comma-separated-values,text/plain', '.tsv')),
                   downloadButton('saveInputObj', 'Save Input Object')
                 ), # end conditionalPanel
              ), # end sidebarPanel

              # the main panel
              mainPanel(
                htmlOutput("message_load")
              )

          ) # end sidebarLayout
    ), # end tabPanel with loading data

    #==============================
    # tab for analysis
    #==============================
  	tabPanel("Run Analysis",
  	 sidebarLayout(
  	   # the left side panel
  		sidebarPanel(
  		  # check box that specifies the input data format
    		checkboxInput('replicates','Analyze with replicates',
    		              value = FALSE),

    		# if analysis with replicates, then specify separator
    		# conditionalPanel(
    		#   condition = "replicates == TRUE",
    		#   tags$span("Specify a separator:"),
    		#   textInput('separator', NULL,
    		#             value = "_", width = "3em",
    		#             placeholder = NULL),
    		# ),



    		# specify if analysis should be performed on nucleotide level data
    	checkboxInput('nucleotideFlag','Use nucleotide level data'),    		# check box that controls the type of analysis -
    		# if the comparison with a reference should be performed or not
    		shinyjs::useShinyjs(),
    		checkboxInput('compareToRef','Compare to reference', value = TRUE),

    		# drop down menu to select a reference
    		selectInput('refSamp', 'Select a reference', choices = list('None'), selected = NULL, multiple = FALSE,
    			selectize = TRUE, width = NULL, size = NULL),


    		selectInput('excludeSamp', 'Select conditions to exclude', choices = list('None'), selected = NULL, multiple = TRUE,
    			selectize = TRUE, width = NULL, size = NULL),
    		textInput('nReads', 'Specify the minimal number of templates (increase in case of large files)',
    		          value = "50", width = NULL, placeholder = NULL),

   		  actionButton('runAnalysis', 'Run Analysis', width = '60%'),

  		),

  		# the main panel
  	  mainPanel(
  		 # textOutput('contents'),
  		  htmlOutput("message_analysis")
  		)
  	 )
	 ),# end tabPanel Analysis

    #=============================
    # tab for saving results of the analysis
    #=============================
    tabPanel("Save Results",
         sidebarLayout(
           # the left side panel
           sidebarPanel(
             tags$p(HTML("<b>Specify thresholds:</b>")),
             numericInput('fdrThr', 'for FDR ', value = 0.05,
                          max = 1, step= 0.01),
             textInput('orThr', 'for odds ratio', value = "5",
                       width = NULL, placeholder = NULL),
             # keep clones that have maximal abundance above that values
             textInput('percentThr', 'for clone abundance (in percent)',
                       value = "0"),
             # a threshold for percent of conditions where a clone has non-zero abundance
             shinyjs::useShinyjs(),
             textInput('condsThr', 'for conditions with non-zero abundance (in percent)',
                       value = "0"),
             downloadButton('saveResults', 'Download Results'),
             downloadButton('saveHeatmaps', 'Save Heatmap')
           ), # end sidebarPanel

           # the main panel
           mainPanel(
             htmlOutput("save_results")
           )

         ) # end sidebarLayout
    ), # end tabPanel with saving data


    #=============================
    # the tab with manual
    #=============================
	 tabPanel("User Manual",
	          tags$iframe(
	            src = "userManual.html",
	            width = "100%",
	            height = "800px",
	            frameborder = 0
	          )
	 )
 ) # end tabsetPanel
)# end of the UI


#====================
# Server
#====================
server <- function(input, output,session) {
  #rm(list = ls())
  # increase file upload limit to 100M
  options(shiny.maxRequestSize=100*1024^2, java.parameters = "-Xmx8000m")
  library(shiny)
  library(tools)
  if (!require(replicateFest)) devtools::install_github("OncologyQS/replicateFest")
  library(replicateFest)
  library("Matrix")
  library(DT)
  library(kableExtra)
  library(dplyr)
  library(shinyjs)
  if (!require(openxlsx)) install.packages("openxlsx")
  if (!require(immunarch)) install.packages("immunarch")
  if (!require(devtools)) install.packages("devtools")
  if (!require(contrast)) install.packages("contrast")
  if (!require(multcomp)) install.packages("multcomp")

#  source("R/manafest_shiny_functions.r")
#  source("R/repManFunctions.r")

  # -----------------------
  # Reactive storage instead of global env
  # -----------------------
  rv <- reactiveValues(
    mergedData = NULL,
    ntData = NULL,
    analysisRes = NULL
  )

  # -----------------------
  # Helper functions to avoid duplicated code between saveResults and saveHeatmaps
  # -----------------------
  get_save_params <- function() {
    sp <- reactiveValuesToList(input)
    sp$orThr     <- suppressWarnings(as.numeric(sp$orThr))
    sp$fdrThr    <- suppressWarnings(as.numeric(sp$fdrThr))
    sp$percentThr<- suppressWarnings(as.numeric(sp$percentThr))
    sp$condsThr  <- suppressWarnings(as.numeric(sp$condsThr))
    sp
  }

  prepare_obj_and_samp <- function(analysisRes) {
    obj <- if (isTRUE(analysisRes$params$nucleotideFlag)) rv$ntData else rv$mergedData
    sampForAnalysis <- if (!is.null(obj)) setdiff(names(obj), c(analysisRes$params$excludeSamp, analysisRes$params$refSamp)) else character()
    list(obj = obj, sampForAnalysis = sampForAnalysis)
  }

  get_positive_clones <- function(analysisRes, obj, sampForAnalysis, saveParams) {
    if (isTRUE(analysisRes$params$compareToRef)) {
      if (isTRUE(analysisRes$params$replicates)) {
        getPositiveClonesReplicates(
          analysisRes$res, obj,
          refSamp = analysisRes$params$refSamp,
          samp = sampForAnalysis,
          excludeCond = analysisRes$params$excludeSamp,
          orThr = saveParams$orThr,
          fdrThr = saveParams$fdrThr,
          percentThr = saveParams$percentThr
        )
      } else {
        getPositiveClones(
          analysisRes$res, obj,
          samp = sampForAnalysis,
          orThr = saveParams$orThr,
          fdrThr = saveParams$fdrThr,
          percentThr = saveParams$percentThr,
          condsThr = saveParams$condsThr
        )
      }
    } else {
      getPositiveClonesFromTopConditions(
        analysisRes$res,
        orThr = saveParams$orThr,
        fdrThr = saveParams$fdrThr,
        percentThr = saveParams$percentThr,
        condsThr = saveParams$condsThr,
        mergedData = obj,
        samp = sampForAnalysis
      )
    }
  }

  #===================
  # read input files
  observeEvent(input$sourceFiles,{
    # remove all previously loaded data and analysis
#    rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    # reset reactive storage (do not touch global env)
    rv$mergedData <- NULL
    rv$ntData <- NULL
    rv$analysisRes <- NULL

    output$message_load <- renderUI({
      req(input$sourceFiles)
      if (length(input$sourceFiles$name) < 2) {
        return(HTML('There should be at least two files to analyze'))
      }

      # Use progress indicator and safe read
      tryCatch({
        withProgress(message = "Reading files...", value = 0.1, {
          files_dir <- dirname(input$sourceFiles$datapath[1])
          incProgress(0.2)
          res <- readMergeSave(files = files_dir, filenames = unlist(input$sourceFiles$name))
          incProgress(0.7)
        })

        if (is.null(res)) stop("readMergeSave returned NULL")

        rv$mergedData <- res$mergedData
        rv$ntData <- res$ntData

        updateSelectInput(session, "refSamp", choices = c('None', names(rv$mergedData)))
        updateSelectInput(session, "excludeSamp", choices = names(rv$mergedData))

        HTML(c('Uploaded files:<br/>', paste(names(rv$mergedData), collapse = '<br/>')))
      }, error = function(e) {
        message("Error reading files: ", e$message)
        HTML(paste("Error in reading files:", e$message))
      })
    })
  })


  #===================
  # load RDA file with the input data
  observeEvent(input$inputObj,{
    # remove all previously loaded data and analysis
 #   rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    # clear reactive storage and load the uploaded .rda into rv
    rv$mergedData <- NULL;
    rv$ntData <- NULL;
    rv$analysisRes <- NULL

    output$message_load <- renderUI({
      # check if there is a file to analyze
      if (length(input$inputObj) == 0)
      {
        HTML('There are no files to analyze. Please upload files')
      }else if(length(input$inputObj$name) > 0 )
      {
        if (tolower(file_ext(input$inputObj$name)) != 'rda') {
          return(HTML('This is not an .rda file. Please check the input file'))
        }
        tmpenv <- new.env(parent = emptyenv())
        load(input$inputObj$datapath, envir = tmpenv)
        if (!exists("mergedData", envir = tmpenv)) {
          return(HTML('There are no input data to analyze. Please check the input file'))
        }
        rv$mergedData <- get("mergedData", envir = tmpenv)
        rv$ntData <- get0("ntData", envir = tmpenv)
        updateSelectInput(session, "refSamp", choices=c('None',names(rv$mergedData)))
        updateSelectInput(session, "excludeSamp", choices=names(rv$mergedData))
        # generate message about loaded samples
        message = c('Uploaded samples:<br/>',paste(names(rv$mergedData), collapse = '<br/>'))
        # return the message as HTML
        HTML(message)
      }
    })
  })

  # enable the Save Input Object button only if data is loaded
  observe({
    shinyjs::toggleState("saveInputObj", !is.null(rv$mergedData))
  })
  # enable the Run Analysis button only if data is loaded
  observe({
    shinyjs::toggleState("runAnalysis", !is.null(rv$mergedData))
  })
  # enable the Download Results and Save Heat Map buttons
  # only if analysis is done
  observe({
    shinyjs::toggleState("saveResults", !is.null(rv$analysisRes))
  })
  observe({
    shinyjs::toggleState("saveHeatmaps", !is.null(rv$analysisRes))
  })

  #===================
  # save the inputData object with uploaded files for further analysis when the Save Input Object button is clicked
  output$saveInputObj <- downloadHandler(
    filename='inputData.rda',
    content=function(file){
      #		load('inputData.rda')
      if (!is.null(rv$mergedData)) {
        #
        mergedData = rv$mergedData
        ntData = rv$ntData
        save(mergedData,ntData,file=file)
      }else{
        print('No data to save')
        output$message_load = renderText('No data to save')
        emptyObj = NULL
        save(emptyObj,file=file)
      }
    })

   #==================
  # an observerEvent for checking the The input with replicates checkbox
  # if it checked, disable the compare to reference and condition threshold fields
  observeEvent(input$replicates, {
    # if the checkbox is selected
    if (input$replicates == TRUE){
      # make compareToRef checked before disabling
      updateCheckboxInput(session, "compareToRef", value = TRUE)
    # disable the Compare to reference option and condition threshold
      shinyjs::disable('compareToRef')
      # and enable
      shinyjs::disable('condsThr')
    }else{ # and enable otherwise
      shinyjs::enable('compareToRef')
      shinyjs::enable('condsThr')
    }
    # if data is loaded
    if(!is.null(rv$mergedData))
    {
      mergedData = rv$mergedData

      if (input$replicates == TRUE){
        # if the input data has replicates, then the reference sample should be selected
        # update conditions
        # check if there are objects to run analysis
         # extract conditions from the file names
          sampAnnot = splitFileName(names(mergedData))
        # check if there are conditions
        if (all(is.na(sampAnnot$condition)))
        {
          output$message_analysis = renderText('There are no conditions to analyze.
                                                 Please check the input file names')
        }else{
          # output a table with sample annotations
          output$message_analysis = renderTable({
            # Create a table using kable
            kable(sampAnnot, format = "html") %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
          }, sanitize.text.function = function(x) x)

          # update inputs in the dropdowns
          updateSelectInput(session, "refSamp",
                            choices=c('None',sampAnnot$condition))
          updateSelectInput(session, "excludeSamp", choices=sampAnnot$condition)

        }
        # make compareToRef checked before disabling
        updateCheckboxInput(session, "compareToRef", value = TRUE)
        # and disable
        shinyjs::disable('compareToRef')

      # if The input with replicates is unchecked
      }else{
        # if the input data does not have replicates, then the reference sample should be selected
        # update conditions
        updateSelectInput(session, "refSamp", choices=c('None',names(mergedData)))
        updateSelectInput(session, "excludeSamp", choices=names(mergedData))
        shinyjs::enable('compareToRef')
        # output sample names in the main panel
        output$message_analysis = renderText(
          HTML(c('Uploaded samples:<br/>',paste(names(mergedData),
                                                collapse = '<br/>'))))

      }
    }
  } )

  #================
  # run analysis with the Run Analysis button is clicked
  observeEvent(input$runAnalysis,{
    #browser()
    output$message_analysis = renderText('Analysis is running...')

    # check if there are objects to run analysis
    if(!is.null(rv$mergedData) && !is.null(rv$ntData))
    {
      # Check if comparison to reference is enabled but no reference sample is selected
      if (input$compareToRef && input$refSamp == 'None') {
        output$message_analysis <- renderText('There is no reference sample. Please select a reference sample.')
        return() # Exit the observer to prevent further execution
      }
      # run analysis on aa level data
        # load object with input data
        if (!input$nucleotideFlag){
          obj = rv$mergedData
        }else{ # run analysis on nucleotide level data
          #if 'Use nucleotide level data' checkbox is selected
          obj = rv$ntData
        }
      # get samples to analyze
      sampNames = names(obj)

        # create an object with results
        # it's a list with two elements:
        # one is the analysis output
        # and the other one is the parameters of the analysis
        analysisRes = list(res = NULL,
                           params = reactiveValuesToList(input))
  #browser()
        #================================
        ## analysis without replicates
        #================================
        if (input$replicates == FALSE)
        {
          # list samples to analyze that excludes reference and other
          # samples that should be excluded from the analysis
          sampForAnalysis = setdiff(sampNames,
                                    c(input$excludeSamp,input$refSamp))

          # select clones to test
          clonesToTest = NULL
#    browser()
          # if the comparison to reference should be included in the analysis
          if(input$compareToRef)
          {
             # create comparing pairs (to refSamp)
              compPairs = cbind(sampForAnalysis,rep(input$refSamp,length(sampForAnalysis)))

              # run Fisher's test
              if(nrow(compPairs) == 1)
              {
                analysisRes$res = list(runFisher(compPairs,obj, clones = clonesToTest, nReadFilter = c(as.numeric(input$nReads),0)))
              }else{
                analysisRes$res = apply(compPairs,1,runFisher,obj, clones = clonesToTest, nReadFilter = c(as.numeric(input$nReads),0))
              }
              #			browser()
              if (!is.null(analysisRes$res))
              {
                names(analysisRes$res) = apply(compPairs,1,paste,collapse = '_vs_')
             }

          }else{
            # if there is no comparison with the reference, then compare only within conditions
  #          if(input$ignoreBaseline) clonesToTest = NULL
            # take only clones that have the number of reads more then nReads and compare with top second and top third conditions
            analysisRes$res = compareWithOtherTopConditions(obj,
                                                      sampForAnalysis,
                                                      nReads = as.numeric(input$nReads),
                                                      clones = clonesToTest)
         }
        } else {# end of analysis without replicates


        #================================
        ## analysis with replicates
        #================================
          # extract conditions from the file names
          sampAnnot = splitFileName(names(obj))
          # add that to the analysisRes object to be used in generating output
          analysisRes$sampAnnot = sampAnnot

          # check if there are conditions
          if (all(is.na(sampAnnot$condition)))
          {
            output$message_analysis = renderText('There are no conditions to analyze.
                                                 Please check the input file names')
          }else{
            # output a table with sample annotations
            output$message_analysis = renderTable({
              # Create a table using kable
              kable(sampAnnot, format = "html") %>%
                kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
            }, sanitize.text.function = function(x) x)

              # update inputs in the dropdowns
              updateSelectInput(session, "refSamp", choices=c('None',sampAnnot$condition))
  #            updateSelectInput(session, "baselineSamp", choices=c('None',sampAnnot$condition))
              updateSelectInput(session, "excludeSamp", choices=sampAnnot$condition)

              # get clones to test
              clonesToTest = getClonesToTest(obj, nReads = as.numeric(input$nReads))

              # check if there are enough clones to analyze
              if (length(clonesToTest)<1)
              {
                output$message_analysis = renderText('There are no clones to analyze. Try to reduce confidence or the number of templates 1')
                return()
              }
              # run the analysis for selected clones
              analysisRes$res = fitModelSet(clonesToTest, obj,
                                       sampAnnot$condition,
                                       excludeCond = input$excludeSamp,
                                       refSamp=input$refSamp,
                                       c.corr=1)
              rownames(analysisRes$res) = analysisRes$res$clone

        }
      }# end of analysis with replicates
       # check if there are any results to save
      if (!is.null(analysisRes$res))
      {
          output$message_analysis = renderText('Analysis is done. Go to the Save results tab')
          rv$analysisRes = analysisRes
#          assign('analysisRes',analysisRes, envir = .GlobalEnv)
      }	else{
          output$message_analysis = renderText('There are no clones to analyze. Try to reduce the number of templates')
      }


    } else{
      output$message_analysis = renderText('There are no data to analyze. Please load files')
    }
  })

  # save results with selected thresholds when the Download Results button is clicked
  output$saveResults <- downloadHandler(
    filename = function() {
      paste0('analysisRes_', Sys.Date(), '.xlsx')
    },
    contentType = "text/xlsx",
    content = function(file) {
      if (is.null(rv$analysisRes)) {
        output$save_results <- renderText('Please click the Run Analysis button')
        return()
      }

      analysisRes <- rv$analysisRes
      saveParams <- get_save_params()
      prep <- prepare_obj_and_samp(analysisRes)
      obj <- prep$obj
      sampForAnalysis <- prep$sampForAnalysis

      posClones <- tryCatch(
        get_positive_clones(analysisRes, obj, sampForAnalysis, saveParams),
        error = function(e) {
          message("Error getting positive clones: ", e$message)
          NULL
        }
      )

      tablesToXls <- list()

      if (is.null(posClones) || nrow(posClones) == 0) {
        tablesToXls$summary <- "There are no positive clone with the specified thresholds"
        output$save_results <- renderText('There are no positive clones. Try to adjust thresholds')
      } else {
        tablesToXls <- createPosClonesOutput(
          posClones,
          obj,
          analysisRes$params$refSamp,
          replicates = analysisRes$params$replicates
        )
      }

      if (isTRUE(analysisRes$params$compareToRef)) {
        resTable <- if (isTRUE(analysisRes$params$replicates)) {
          createResTableReplicates(
            analysisRes$res,
            obj,
            percentThr = saveParams$percentThr,
            refSamp = analysisRes$params$refSamp,
            orThr = saveParams$orThr,
            fdrThr = saveParams$fdrThr
          )
        } else {
          createResTable(
            analysisRes$res, obj,
            orThr = saveParams$orThr,
            fdrThr = saveParams$fdrThr,
            percentThr = saveParams$percentThr,
            condsThr = saveParams$condsThr,
            saveCI = FALSE
          )
        }

        if (is.null(resTable)) {
          tablesToXls$ref_comparison_only <- data.frame(res = 'There are no significant clones', stringsAsFactors = FALSE)
        } else if (nrow(resTable) > 0) {
          tablesToXls$ref_comparison_only <- resTable
        }
      }

      # parameters sheet
      s <- if (!is.null(obj)) names(obj) else character()
      productiveReadCounts <- if (!is.null(obj) && length(s) > 0) sapply(obj, sum) else numeric(0)

      param <- c(
        "Data with replicates",
        "Reference sample",
        "Excluded samples",
        "Compare to reference",
        "n template threshold",
        "FDR threshold",
        "OR threshold",
        "percent threshold",
        "percent non-zero conditions",
        "Nucleotide level analysis",
        "n samples",
        paste(s, "n templates", sep = "_")
      )

      value <- c(
        analysisRes$params$replicates,
        analysisRes$params$refSamp,
        paste(analysisRes$params$excludeSamp, collapse = ", "),
        analysisRes$params$compareToRef,
        analysisRes$params$nReads,
        saveParams$fdrThr,
        saveParams$orThr,
        saveParams$percentThr,
        ifelse(is.null(saveParams$condsThr), "", saveParams$condsThr),
        analysisRes$params$nucleotideFlag,
        length(s),
        if (length(productiveReadCounts) > 0) productiveReadCounts[s] else numeric(0)
      )

      tablesToXls$parameters <- data.frame(param = param, value = value, stringsAsFactors = FALSE)

      # write out xlsx
      saveResults(tablesToXls, file)
      output$save_results <- renderText('The results are saved')
    }
  )

  # save heatmaps with selected thresholds when the Download Heatmaps button is clicked
  output$saveHeatmaps <- downloadHandler(
    filename = function() {
      paste0('heatmap_', Sys.Date(), '.pdf')
    },
    contentType = "image/pdf",
    content = function(file) {
      if (is.null(rv$analysisRes)) {
        output$save_results <- renderText('Please click the Run Analysis button')
        return()
      }

      analysisRes <- rv$analysisRes
      saveParams <- get_save_params()
      prep <- prepare_obj_and_samp(analysisRes)
      obj <- prep$obj
      sampForAnalysis <- prep$sampForAnalysis

      posClones <- tryCatch(
        get_positive_clones(analysisRes, obj, sampForAnalysis, saveParams),
        error = function(e) {
          message("Error getting positive clones for heatmap: ", e$message)
          NULL
        }
      )

      if (is.null(posClones) || nrow(posClones) <= 1) {
        m <- "There are not enough positive clones to plot. Try to adjust thresholds"
        output$save_results <- renderText(m)
        pdf(file)
        plot.new()
        text(0.5, 0.5, m)
        dev.off()
        return()
      }

      ref <- if (isTRUE(analysisRes$params$replicates)) NULL else analysisRes$params$refSamp

      tryCatch({
        makeHeatmaps(unique(posClones$clone), obj,
                     refSamp = ref,
                     fileName = file, size = 7)
        output$save_results <- renderText('The heat map is saved')
      }, error = function(e) {
        message("Error creating heatmap: ", e$message)
        output$save_results <- renderText(paste("Heatmap generation failed:", e$message))
      })
    }
  )
  # Register cleanup code when session ends
  session$onSessionEnded(function() {
    gc()
#    message("Global environment cleaned up.")
  })

}# end the sever function

shinyApp(ui = ui, server = server)
