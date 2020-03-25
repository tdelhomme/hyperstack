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
      ---- INTERCEPT METHOD ----

      Mandatory arguments:
      --vcf                       - input VCF with FPRF geno stat (should be in bgzip format and indexed with tabix)
      --bin_path                  - path to bin folder containing functions to load
      --help                      - print this text

      Optional arguments:
      --mutation_type             - type of mutations that would be used(all(default), exonic, intergenic or intronic)
      --no_plots                  - to do not output plots
      --fdr_range                 - fdr range in format: fdr_min,fdr_max (default is 0,1)
      --output_folder             - path to the output folder (default=vcf_name_intercept_output)
      --genome                    - genome version (default=hg18)
      --sm                        - sample identifier (default=input vcf file name)
      --germline_mutations        - txt file containing germline mutations (chr pos ref alt sm)
      --minQVAL                   - filtering calls with QVAL<minQVAL (default=10)

      example: Rscript intercept_method.r --vcf=file.bgz \n\n")
  
  q(save="no")
}

if(is.null(args$vcf)) {stop("Option --vcf should be provided")} else{vcf=args$vcf}
if(is.null(args$bin_path)) {stop("Option --bin_path should be provided")} else{bin_path=args$bin_path}
if(is.null(args$mutation_type)) {args$mutation_type="all"}
if(! args$mutation_type %in% c("all","exonic","intergenic","intronic")) {stop("Option --mutation_type should be: all, exonic, intergenic or intronic")} else {mutation_type=args$mutation_type}
if(is.null(args$no_plots)) {no_plots=FALSE} else {no_plots=TRUE}
if(is.null(args$fdr_range)) {fdr_range=c(0,1)} else {fdr_range=as.numeric(unlist(strsplit(args$fdr_range,",")))}
if(is.null(args$output_folder)) { output_folder="." } else {output_folder = args$output_folder}
system(paste("mkdir -p",output_folder,sep=" "))
if(is.null(args$genome)) {genome="hg18"} else {genome=args$genome}
if(is.null(args$sm)) {sm=substr(basename(vcf), 1, 21)} else {sm=args$sm}
if(is.null(args$germline_mutations)) {germline_mutations=NULL} else {germline_mutations=as.character(read.table(args$germline_mutations,h=F)[,1])}
if(is.null(args$minQVAL)) {minQVAL=10} else {minQVAL=args$minQVAL}

# loading required libraries
suppressMessages(library(VariantAnnotation))
suppressMessages(library(Hmisc))
suppressMessages(library(robustbase))
suppressMessages(library(MutationalPatterns))
suppressMessages(library(BSgenome))
suppressMessages(library(RCurl))

source(paste(bin_path,"/objects.r",sep=""))

# set the genome
ref_genome = available.genomes()[which(grepl(genome,available.genomes()) & ! grepl("masked",available.genomes()))]
if (! is.element(ref_genome, installed.packages()[,1])) BiocManager::install(ref_genome)
library(ref_genome, character.only = TRUE)

# compute fdr ranges to fit the linear regression
min_fdr = fdr_range[1] ; max_fdr = fdr_range[2]

# format the input vcf into variantannotation vcf file and keep interesting chromosomes
vcf = open(VcfFile(vcf,  yieldSize=500000))
calls = readVcf(vcf, genome)

while(dim(calls)[1] != 0){
  print(paste("Starting a new chunk at:", date(), sep=" "))
  print(paste("Number of mutations in the chunk:", nrow(calls), sep=" "))
  chrs = as.character(seqnames(rowRanges(calls,"seqnames")))
  human_chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22","chrX", "chrY")
  calls = calls[which(chrs %in% human_chrs),]
  seqlevels(calls) = human_chrs
  calls = calls[which(apply(geno(calls, "QVAL"), 1, max) >= minQVAL), ] # here can't use QUAL if the vcf was separated into different pieces (QVAL is not recalculated)
  
  # remove germline mutations if in input
  if(!is.null(germline_mutations)) calls = calls[which(! rownames(calls) %in% germline_mutations),]
  
  # filter the calls on the requested mutation type
  #refgene = unlist(info(calls)$Func.refGene)
  #if(mutation_type != "all") calls = calls[which(refgene == mutation_type)]
  
  # get the 3nucleotides context
  nut3_context = type_context(rowRanges(calls), ref_genome)
  all_mut1 = paste(substr(nut3_context$context,1,1), "[", nut3_context$types, "]", substr(nut3_context$context,3,3), sep = "")
  all_mut2 = rep(all_mut1, each = length(samples(header(calls))))
  # compute the FDRs
  all_fdr1 = as.vector(unlist(t(geno(calls)[["FPRF"]]))) # use the statistic given by needlestack
  
  # get only non na mutations i.e. mutations with a FPRF value so for which we have apply the model (typically QVAL>20)
  all_fdr_tot = all_fdr1[which(!is.na(all_fdr1))]
  all_muts_tot = all_mut2[which(!is.na(all_fdr1))]
  
  if(! exists("all_fdr")) {all_fdr=all_fdr_tot} else {all_fdr = c(all_fdr, all_fdr_tot)}
  if(! exists("all_muts")) {all_muts=all_muts_tot} else {all_muts = c(all_muts, all_muts_tot)}
  
  calls = readVcf(vcf, genome)
}

# initiate the results
by_par = 0.15
fdr_ranges = seq(min_fdr, max_fdr + by_par, by=by_par)
fdrs = fdr_ranges[1:(length(fdr_ranges)-1)] + by_par/2

dat_counts = dat_total = as.data.frame(matrix(0, nrow = length(fdr_ranges)-1, ncol = length(TRIPLETS_96)))
colnames(dat_counts) = colnames(dat_total) = TRIPLETS_96

# compute the proportion of subsitutions (in column) by fdr bin (in row)
for(sub in TRIPLETS_96){
  # local FDRs
  dat_counts[,sub] = unlist(lapply(1:(length(fdr_ranges)-1), function(i){
    length(which(all_fdr >= fdr_ranges[i] & all_fdr <= fdr_ranges[i+1] & all_muts == sub))
  }))
  dat_total[,sub] = unlist(lapply(1:(length(fdr_ranges)-1), function(i){
    length(which(all_fdr >= fdr_ranges[i] & all_fdr <= fdr_ranges[i+1] ))
  }))
}

# plot the results and output the estimated robust linear regression
list_mod = list()
list_r2 = list()
if(! no_plots) pdf(paste(output_folder, "/local_fdrs_",sm,"_",mutation_type,".pdf", sep=""), 12, 10)
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
for(sub in TRIPLETS_96){
  assign(paste("wald",sub,sep="_"), lapply(1:nrow(dat_counts), function(i) binconf(dat_counts[i,sub], dat_total[i,sub], method="asymptotic")))
  pcts = dat_counts[,sub] / dat_total[,sub]
  pcts[which(is.na(pcts))] = 0
  plot(fdrs, pcts, pch=19, main=paste(sub," n=",table(all_muts)[sub],sep=""), xlab="local fdr bin", ylab="percent of all mutations with x local fdr", ylim=c(0,2*max(pcts)))
  lines(fdrs, pcts, lwd=2)
  if(max(pcts)!= 0){
    msg.UCV <- "XXX" # to resolve a bug, see https://stackoverflow.com/questions/56334077
    mod <- robust::lmRob(Pc ~ Fd, data = data.frame(Pc=pcts, Fd=fdrs))
    list_mod[[sub]] = ifelse(as.numeric(mod$coefficients["(Intercept)"]) < 0 , 0, as.numeric(mod$coefficients["(Intercept)"]))
    list_r2[[sub]] = as.numeric(cor.test(pcts, as.numeric(mod$coefficients["Fd"]) * fdrs + as.numeric(mod$coefficients["(Intercept)"]))$estimate)
    abline(b = as.numeric(mod$coefficients["Fd"]), a = as.numeric(mod$coefficients["(Intercept)"]), col="red", lwd=2)
    for(i in 1:(nrow(dat_counts)-1)){
      wald1 = get(paste("wald",sub,sep="_"))[[i]]
      wald2 = get(paste("wald",sub,sep="_"))[[i+1]]
      segments(x0=fdr_ranges[i] + by_par/2, x1=fdr_ranges[i+1] + by_par/2, y0=wald1[3], y1=wald2[3], lty=2)
      segments(x0=fdr_ranges[i] + by_par/2, x1=fdr_ranges[i+1] + by_par/2, y0=wald1[2], y1=wald2[2], lty=2)
    }
  } else {list_mod[[sub]] = 0 ; list_r2[[sub]] = NA}
} ; mtext(paste(sm,mutation_type,sep=" - "), outer = TRUE, cex = 1.5)

# return the list of estimated linear models for each substitution
assign(paste("rlm_ITC",sm,mutation_type,sep="_"), list_mod)
save(list = ls(pattern="rlm_ITC_"), file=paste(output_folder,"/","rlm","_",sm,"_",mutation_type,"_ITC.Rdata",sep=""))
assign(paste("r2_ITC",sm,mutation_type,sep="_"), list_r2)
save(list = ls(pattern="r2_ITC_"), file=paste(output_folder,"/","r2","_",sm,"_",mutation_type,"_ITC.Rdata",sep=""))
if(! no_plots) dev.off()
