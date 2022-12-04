#' The toolbox for analysing single-cell protein intensity from fluorescence flow cytometry (FFC).
#'
#' @description This function is an end-to-end toolbox for analysing the data from FFC in FCS format, including data normalisation and cluster analysis.
#'
#' @param FCSpath Path to the input directory where input FCS files are stored (one file per well). Will look for FCS files recursively in that directory.
#' @param Outpath Path to the output directory where intermediate results and final results will be stored
#' @param batch_key Data from plate-based experiments (batch_key="plate") or not (batch_key="batch")
#' @param bkb.v A vector of the names of the backbone markers (MUST match to the names in the FCS file).
#' @param bkb.upper.quntile The cut-off (default = 0.9) for selecting cells used for estimating the parameter of signal.
#' @param bkb.lower.quntile The cut-off (default = 0.1) for selecting cells used for estimating the parameters of noise.
#' @param bkb.min.quntile The cut-off (default = 0.01) for omitting the cells with the smallest values to minimise the impact of outliers.
#' @param trans.dat Method for transforming raw intensities (default: Logicle).
#'
#' @author Hsiao-Chi Liao
#'
#' @export
#' @return Normalised backbone measurements. Cluster analysis for normalised backbones. Graphs will be provided if set visualisation = TRUE.
#'
#' @usage
#' pipBackbone(FCSpath="/PathToFCSfiles/",
#' Outpath="/PathToOutputFolder/Output/",
#' batch_key="plate",
#' bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
#' "F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"),
#' bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01,
#' trans.dat="cent.lgc")
#'
#'
pipBackbone <-
function(FCSpath,
                              Outpath,
         batch_key="plate",
                              bkb.v,
                              bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01,
                              trans.dat="cent.lgc")
  {

    ##/!\ Potentially add a check here to make sure parameters are consistent with FCS files

    settings <- initialize(
      path_to_fcs=FCSpath,
      path_to_output=Outpath,
      verbose=TRUE
    )
    paths <- settings$paths

    ###

    ##1.fcs to rds
    fcs_to_rds(paths, batch_key=batch_key)

    ##2-1.bkc bkb
    bkc_bkb(paths, bkb.v,
            bkb.upper.quntile=0.9, #cells used for estimating parameter of signal
            bkb.lower.quntile=0.1, #cells used for estimating parameters of noise
            bkb.min.quntile=0.01) #the lowest 1% of values will not be used to minimise the impact of outliers on sig

    ##3-0.initM
    initM(paths, bkb.v,
          trans.dat="cent.lgc") #may give "cent.log.bkc" in the future

    ##3-1.rmWellEffect
    rmWellEffect(paths,
                 visualisation = TRUE)

    ##7.cluster analysis bkbOnly
    cluster.analysis.bkbOnly(paths,
                     bkb.v)

    message("\tCompleted!")
  }
