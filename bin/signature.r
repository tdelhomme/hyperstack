args <- commandArgs(TRUE)
parseArgs <- function(x) {
  res = strsplit(sub("^--", "", x), "=")
  if(length(unlist(res))==1) res[[1]][2]=""
  return(res)
}

argsL <- as.list(do.call("cbind", parseArgs(args))[c(F,T)])
names(argsL) <- as.list(do.call("cbind", parseArgs(args))[c(T,F)])
args <- argsL;rm(argsL)

if(! is.null(args$help)) {
  cat("
      ---- Estimation of mutational signature from hyperstack method on single cell data ----
      Mandatory arguments:
      NONE
      
      Optional arguments:
      --ITC_folder                - path to the output folder ITC from our hyperstack method
                                    containing the rlm object (list of 96 intercept estimations)
                                    (default: .)
      --output_folder             - path to the output folder (default: signatures)
      --sigbin_folder             - path to the signature bin folder, i.e. containing functions to
                                    load/script to source to perform sig estimation (default: .)
      --help                      - print this text
      example: Rscript signature.r --ITC_output_folder=/path/to/itc/folder \n\n")
  
  q(save="no")
}

if(is.null(args$ITC_folder)) {ITC_folder="."} else {ITC_folder=as.character(args$ITC_folder)}
if(is.null(args$output_folder)) {output_folder="signatures"} else {output_folder=as.character(args$output_folder)}
if(is.null(args$sigbin_folder)) {sigbin_folder="."} else {sigbin_folder=as.character(args$sigbin_folder)}

library(Biostrings) # added this to reverse DNA string
library(magrittr)

for( f in list.files(path = ITC_folder, pattern = "rlm_BC+")){
  load(paste(ITC_folder,"/",f,sep=""))
  bc = unlist(strsplit(f, "_"))[2]
  r = get(ls(pattern = bc))
  #pow10 = (10^(1:10)) [ which((min(unlist(r[which(r >0)])) * 10^(1:10) > 1))[1] ]
  #r2 = t(as.matrix(unlist(lapply(r, function(x) round(x * pow10))))) # in number of mut
  r2 = t(as.matrix(unlist(r))) # in prop of mut
  if(exists("triSig.dataM")){
    r2 = r2[,colnames(triSig.dataM)] # ordering, just in case
    triSig.dataM = rbind(triSig.dataM, r2)
    rownames(triSig.dataM)[nrow(triSig.dataM)] = bc
  } else {
    triSig.dataM = r2
    rownames(triSig.dataM) = bc
  }
}
subs96 = colnames(triSig.dataM)
subs96 = subs96 %>% gsub("\\[","", .) %>% gsub("\\]","",.) %>% gsub(">","",.)
subs96 = paste(substr(subs96,1,2), substr(subs96,4,4), substr(subs96,3,3), sep = "")     
subs96[49:96] = unlist(lapply(subs96[49:96], function(x) as.character(Biostrings::complement(DNAString(x)))))
subs96 = paste(substr(subs96,1,3), substr(subs96,4,4), sep=">")
colnames(triSig.dataM) = subs96

setwd(sigbin_folder)
source("NMF_signatureDiscovery_hyperstack.R")

#plot signatures
#pdf(paste(input_folder,"/signatures_plot.pdf", sep=""))
plot_grid(plotlist = plots, ncol = 1)
#dev.off()

