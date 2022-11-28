#' Converting FCS files to RDS files
#' 
#' This function has been designed to convert the raw FCS files to data matrix and export to RDS files.
#' 
#' @param paths Character vector of paths to store intput, intermediary results, outputs...
#' @param batch_key Data from plate-based experiments (batch_key="plate") or not (batch_key="batch")
#' 
#' @author Hsiao-Chi Liao
#' 
#' @importFrom flowCore read.flowSet keyword parameters
#' @importFrom Biobase pData exprs
#' @importFrom stringr str_extract
#' 
#' @return Generating fcs_metadata_df.rds and fcs_rawInten_mt.rds files.
#' 
fcs_to_rds <-
function(paths, batch_key=batch_key) #inpath with a set of fcs files
  {
    # FCSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/sample_data/"
    # RDSpath="/Users/chelsea/Desktop/PhD/ProteinExpression/InFlow/Programming/infinity/R_CODE/72.mk_pkg/data/rds/"
    
    ## re-make large mt and consider SSC and FSC
    # import fcs.impu.raw files ###
    fs <- read.flowSet(path = file.path(paths["input"]))
    name.desc.fs <- pData(parameters(fs[[1]]))
    name.desc.fs$desc[which(is.na(name.desc.fs$desc))] <- name.desc.fs$name[which(is.na(name.desc.fs$desc))]
    
    ####
    ##check the keywords
    # sink(paste0(outpath, "keyword.txt"))
    # options(max.print=1000000)
    # i=1
    # print(keyword((fs[[i]])))
    # sink()
    
    #### data - raw ####
    fcs.raw.mt.list = fcs.raw.meta.list = list()
    for(i in 1:length(fs)){
      fcs.raw.mt.list[[i]] <- exprs(fs[[i]]) #22 columns
      fcs.raw.meta.list[[i]] <- cbind(NO.in.well=(1:dim(fcs.raw.mt.list[[i]])[1]),
                                      Filenam=unlist(keyword(fs[[i]], "GUID")))
    }
    fcs.raw.mt <- do.call(rbind, fcs.raw.mt.list) #266000*2
    colnames(fcs.raw.mt) <- name.desc.fs$desc[match(colnames(fcs.raw.mt), name.desc.fs$name)]
    ####
    
    # making metadata ###
    fcs.raw.meta.df <- as.data.frame(do.call(rbind, fcs.raw.meta.list)) #266000*2
    fcs.raw.meta.df$NO.in.well <- as.integer(fcs.raw.meta.df$NO.in.well)
    #use lapply for multiple changing class: lapply(mydf[,2:3], as.factor)
    
    if(batch_key == "plate"){
      #plate
      fcs.raw.meta.df$Plate <- str_extract(fcs.raw.meta.df$Filenam, "Plate[0-9]")
      
      #well
      fcs.raw.meta.df$Well <- str_extract(fcs.raw.meta.df$Filenam, "[A-H]\\d{1,2}") #[:digit:], or the shorthand \\d
      orig.lab <- c(paste0(rep(c("A","B","C","D","E","F","G","H"), each = 9), c(1:9)),
                    paste0(rep(c("A","B","C","D","E","F","G","H"), each = 9), c(10:12)))
      new.lab <- c(paste0(rep(c("A","B","C","D","E","F","G","H"), each = 9), 0, c(1:9)),
                   paste0(rep(c("A","B","C","D","E","F","G","H"), each = 9), c(10:12)))
      map <- setNames(new.lab, orig.lab)
      fcs.raw.meta.df$Well <- map[fcs.raw.meta.df$Well]
      
      #column
      fcs.raw.meta.df$Column <- paste0("Col.", sub("[A-H]", "", fcs.raw.meta.df$Well))
      # fcs.raw.meta.df$Column <- sub("[A-G]", "", fcs.raw.meta.df$Well)
      # orig.lab <- c("1","2","3","4","5","6","7","8","9","10","11","12")
      # new.lab <- paste0("Col.", 
      #                   c("01","02","03","04","05","06","07","08","09","10","11","12"))
      # map <- setNames(new.lab, orig.lab)
      # fcs.raw.meta.df$Column <- map[fcs.raw.meta.df$Column]
      
      #row
      fcs.raw.meta.df$Row <- sub("\\d{1,2}", "", fcs.raw.meta.df$Well)
      orig.lab <- c("A","B","C","D","E","F","G","H")
      new.lab <- paste0("Row.", 
                        c("01","02","03","04","05","06","07","08"))
      map <- setNames(new.lab, orig.lab)
      fcs.raw.meta.df$Row <- map[fcs.raw.meta.df$Row]
      
      #well.lab
      fcs.raw.meta.df$Well.lab <- paste0(sub("late", "", fcs.raw.meta.df$Plate),
                                         "_",
                                         fcs.raw.meta.df$Well)
      
      #ordering before giving NO.in.all
      # unique(fcs.raw.meta.df$Well.lab)
      ord.fcs.raw.meta.df <- fcs.raw.meta.df[order(fcs.raw.meta.df$Well.lab),]
      # unique(ord.fcs.raw.meta.df$Well.lab)
      ord.fcs.raw.mt <- fcs.raw.mt[order(fcs.raw.meta.df$Well.lab),]
    }else{
      # fcs.raw.meta.df$batch 
    }

    
    #global ID for each cell
    ord.fcs.raw.meta.df$NO.in.all <- 1:nrow(ord.fcs.raw.mt)
    
    # tail(ord.fcs.raw.meta.df); unique(ord.fcs.raw.meta.df$Well.lab)
    
    #out rds
    ord.fcs.raw.meta.df.out <- ord.fcs.raw.meta.df[,c(8,1,3,5,6,4,7)]
    saveRDS(ord.fcs.raw.meta.df.out, file = file.path(paths["intermediary"], "/fcs_metadata_df.rds"))
    saveRDS(ord.fcs.raw.mt, file = file.path(paths["intermediary"], "/fcs_rawInten_mt.rds"))
  }
