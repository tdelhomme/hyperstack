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
      ---- Application of a trained Random Forest to assign FDR to needlestack variant calls ----

      Mandatory arguments:
      --vcf                       - VCF to apply the RF model (should be annotated with annovar and in bgzip format/indexed with tabix)
      --model                     - a trained random forest model (in .Rdata format)
      --bin_path                  - path to bin folder containing functions to load
      --help                      - print this text

      Optional arguments:
      --output_folder             - path to the output folder (default=vcfname_RF_output)
      --genome                    - genome version (default=hg18)
      --minQVAL                   - filtering calls with QVAL<minQVAL (default=20)
      --features                  - vcf features used to train the RF model
                                    (default=QVAL,AO,AF,DP,ERR,QUAL,RVSB,FS)
      --mappability_file          - file containing mappability scores for each variant in input vcf (line format: )

      example: Rscript FDR_RF_apply.r --vcf=myvcf.bgz --model=myRF.Rdata \n\n")

  q(save="no")
}

set.seed(98)

if(is.null(args$vcf)) {stop("Option --vcf should be provided")} else{vcf=args$vcf}
if(is.null(args$bin_path)) {stop("Option --bin_path should be provided")} else{bin_path=args$bin_path}
if(is.null(args$model)) {stop("Option --model should be provided")} else{model=args$model}
load(model) # the model whould be an object named "rf"
if(is.null(args$genome)) {genome="hg18"} else {genome=args$genome}
if(is.null(args$output_folder)) { output_folder="."} else {output_folder = args$output_folder}
system(paste("mkdir -p",output_folder,sep=" "))
out_vcf = paste(output_folder, "/", paste( sub(".vcf.gz", "", sub('.vcf.bgz', '', basename(vcf))), "RF_needlestack.vcf", sep="_"), sep="")
if(is.null(args$features)) {features=c("QVAL","AO","AF","DP","ERR","QUAL","RVSB","FS")} else {features=as.character(unlist(strsplit(args$features,",")))}
if(is.null(args$minQVAL)) {minQVAL=20} else {minQVAL=args$minQVAL}
if(is.null(args$mappability_file)) {mappability = FALSE} else {
  print("INFO: mappability scores have been provided")
  mappability = TRUE
  map_dat = read.table(args$mappability_file, header = T, stringsAsFactors = F)
  map_vect = map_dat$MAPPABILITY
  names(map_vect) = paste(map_dat$CONTIG, map_dat$START, sep="-")
}

suppressMessages(library(VariantAnnotation))
suppressMessages(library(randomForest))
source(paste(bin_path,"/functions.r",sep=""))

vcf = open(VcfFile(vcf, yieldSize = 500000))
all_calls = readVcf(vcf, genome)

while(dim(all_calls)[1] != 0) {
  print(paste("Starting a new chunk at:", date(), sep=" "))
  all_calls = all_calls[which(apply(geno(all_calls, "QVAL"), 1, max) >= minQVAL), ] # here can't use QUAL if the vcf was separated into different pieces (QVAL is not recalculated)

  ### populate the table of all mutations with features and ethnicities ###

  n_samples = length(samples(header(all_calls)))

  kept_variants = which(as.vector(t(geno(all_calls)$QVAL)) >= minQVAL)
  all_mut_table = data.frame(matrix(NA, nrow = length(kept_variants), ncol = length(features)))
  colnames(all_mut_table) = features
  rownames(all_mut_table) = paste(rep(rownames(all_calls), each=n_samples), rep(samples(header(all_calls)),n_samples), sep="\\")[kept_variants]

  # assign features
  for (f in features){
    if( f %in% names(geno(all_calls))){ # start with genotype because a variable can have same name in both geno and info (prioritize geno)
      all_mut_table[,f] = as.vector(t(geno(all_calls)[[f]]))[kept_variants] # do this in each if, in order to do not keep na values for features not found
    } else { # only do this if not found in geno, otherwise vector will be replace by the info one
      if(f %in% colnames(info(all_calls))){
        all_mut_table[,f] = rep(info(all_calls)[,f], each = n_samples)[kept_variants] # to get same structure than as.vector(t(geno))
      } else {
        all_mut_table[,f] = rep(as.data.frame(rowRanges(all_calls)[,f])[,f], each = n_samples)[kept_variants]
      }
    }
  }
  if(mappability){
    id_sub = paste(as.character(seqnames(rowRanges(all_calls,"seqnames"))), start(ranges(rowRanges(all_calls,"seqnames"))), sep="-")
    map_scores = as.numeric(map_vect[id_sub])
    all_mut_table[,"CRGmap"] = rep(map_scores, each=n_samples)[kept_variants]
  }
  if("RVSB" %in% features) all_mut_table[which(all_mut_table$RVSB <0.5),"RVSB"]=0.5

  # apply the random forest model
  all_mut_table$FPprob = predict(rf, all_mut_table, type="prob")[,"FP"] # extract probability to be false to have the FPR
  FPRF = rep(NA, length(rep(rownames(all_calls), each=n_samples)))
  FPRF[kept_variants] = all_mut_table$FPprob

  # annotate the VCF with FPRF statistic
  geno(header(all_calls))["FPRF",]=list("A","Float","False Probability of being a variant for a trained Random-Forest model")
  geno(all_calls)$FPRF = matrix(data=FPRF, nrow=nrow(all_calls), byrow = T)

  # write out the annotated VCF file
  con = file(out_vcf, open = "a")
  writeVcf(all_calls, con)
  close(con)
  all_calls = readVcf(vcf, genome)
}
