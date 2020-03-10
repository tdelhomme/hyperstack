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
      ---- Random Forest algorithm to assign FDR to needlestack variant calls ----

      Mandatory arguments:
      --vcf                       - VCF to train/apply the RF model (should be annotated with annovar and in bgzip format/indexed with tabix)
      --help                      - print this text

      Optional arguments:
      --output_folder             - path to the output folder (default=vcfname_RF_output)
      --genome                    - genome version (default=hg18)
      --minQVAL                   - filtering calls with QVAL<minQVAL (default=20)
      --features                  - vcf features to train the model
                                    (default=QVAL,AO,AF,DP,ERR,QUAL,RVSB,FS)
      --ethnicity                 - to use ethnicity to train (use --ethnicity=TRUE)
      --mappability_file          - file containing mappability scores for each variant in input vcf (line format: )
      --normal_id_pattern         - pattern to identify normal samples (Default: BN, as in Li et al 2012)
      --min_normal_cells          - minimum nb of normal samples to set a call as TP (Default: 3, can be 1 for bulk)
      
      example: Rscript FDR_RF_train.r --vcf=myvcf.bgz \n\n")

  q(save="no")
}

set.seed(98)

if(is.null(args$vcf)) {stop("Option --vcf should be provided")} else{vcf=args$vcf}
if(is.null(args$genome)) {genome="hg18"} else {genome=args$genome}
if(is.null(args$output_folder)) { output_folder="."} else {output_folder = args$output_folder}
system(paste("mkdir -p",output_folder,sep=" "))
out_vcf = paste(output_folder, "/", paste( sub(".vcf.gz", "", sub('.vcf.bgz', '', basename(vcf))), "RF_needlestack.vcf", sep="_"), sep="")
if(is.null(args$minQVAL)) {minQVAL=20} else {minQVAL=args$minQVAL}
if(is.null(args$features)) {features=c("QVAL","AO","AF","DP","ERR","QUAL","RVSB","FS")} else {features=as.character(unlist(strsplit(args$features,",")))}
if(is.null(args$ethnicity)) {ethnicity = FALSE} else {ethnicity = TRUE}
if(is.null(args$normal_id_pattern)) {normal_id_pattern = "BN"} else {normal_id_pattern = a.character(args$normal_id_pattern)}
if(is.null(args$min_normal_cells)) {min_normal_cells=3} else {min_normal_cells=args$min_normal_cells}
if(is.null(args$mappability_file)) {mappability = FALSE} else {
  print("INFO: mappability scores have been provided")
  mappability = TRUE
  map_dat = read.table(args$mappability_file, header = T, stringsAsFactors = F)
  map_vect = map_dat$MAPPABILITY
  names(map_vect) = paste(map_dat$CONTIG, map_dat$START, sep="-")
}

suppressMessages(library(VariantAnnotation))
suppressMessages(library(randomForest))
suppressMessages(library(caret))
suppressMessages(library(ROCR))
suppressMessages(library(parallel))

vcf = open(VcfFile(vcf, yieldSize=500000))
all_calls = readVcf(vcf, genome)

while(dim(all_calls)[1] != 0) {
  print("new chunk")
  
  # in case minQVAL is lower than --min_qval used for the calling
  # don't use QUAL if the vcf was separated into different pieces (QVAL is not recalculated)
  # all_calls = all_calls[which(apply(geno(all_calls, "QVAL"), 1, max) >= minQVAL), ]

  ### populate the table of all mutations with features and ethnicities ###

  n_samples = length(samples(header(all_calls)))

  kept_variants = which(as.vector(t(geno(all_calls)$QVAL)) >= minQVAL)
  all_mut_table = data.frame(matrix(NA, nrow = length(kept_variants), ncol = length(features)))
  colnames(all_mut_table) = features
  rownames(all_mut_table) = paste(rep(rownames(all_calls), each=n_samples), rep(samples(header(all_calls)),n_samples), sep="\\")[kept_variants]

  all_pop = c("EAS","AMR","FIN","OTH","SAS","NFE","AFR")
  for(pop in all_pop){ # compute for each pop the relative proportion of the variants (0.5=50% of samples with that SNP were from this pop)
    assign(paste("rp",pop,sep="_"), info(all_calls)[paste("ExAC",pop,sep="_")][,1] /
             rowSums(data.frame(info(all_calls)$ExAC_AMR, info(all_calls)$ExAC_EAS, info(all_calls)$ExAC_FIN, info(all_calls)$ExAC_OTH, info(all_calls)$ExAC_SAS,
                 info(all_calls)$ExAC_NFE, info(all_calls)$ExAC_AFR), na.rm=T)
             )
  }
  # compute the number of variant that are unique for the population
  all_rp = unlist(lapply(all_pop, function(p) {
    rp = rep(get(paste("rp", p, sep="_")), each=n_samples)[kept_variants]
    length(which(rp>=0.99))
  }))
  ethn = all_pop[which.max(all_rp)]

  sm_ethn = rep(get(paste("rp", ethn, sep="_")), each=n_samples)[kept_variants]

  pdat = data.frame(get(paste("rp", all_pop[which(all_pop != ethn)][1], sep="_")))
  for(p in all_pop[which(all_pop != ethn)][2:length(all_pop[which(all_pop != ethn)])]){
    pdat = cbind(pdat, get(paste("rp", p, sep="_")))
  }
  other_ethn = rep(apply(pdat, 1, function(r){if(sum(!is.na(r))==0) {NA} else {max(r, na.rm=T)}}), # compute max relative proportion of all other ethnicities
                   each=n_samples)[kept_variants]

  exac_all = rep(info(all_calls)$ExAC_ALL, each=n_samples)[kept_variants]

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

  if( ! ethnicity ){
    print("INFO: using recurrence of mutation to build the sets of TP/FP")
    normal_calls = all_calls[,which(grepl(normal_id_pattern, colnames(all_calls)))]
    nbq = apply(geno(normal_calls)[["QVAL"]], 1, function(r) length(which(r>=minQVAL)))
    all_nbq = rep(nbq, each = n_samples)[kept_variants]
    # assign status
    all_mut_table$status = NA # status as NA is for variants not used in the training
    all_mut_table[which(all_nbq>=min_normal_cells),"status"] = "TP" # if found in at least 3 normal cells
    all_mut_table[which(all_nbq==1 & (exac_all == 0 | is.na(exac_all))),"status"] = "FP" # if found in only 1 normal cell
  }
    
  if( ethnicity ){
    print("INFO: using ethnicity of the sample to build the sets of TP/FP")
    # assign status
    all_mut_table$status = NA # status as NA is for variants not used in the training
    all_mut_table[which(sm_ethn>=0.99),"status"] = "TP" # here add & exac_all>=0.001 if want to filter rare variants
    all_mut_table[which(other_ethn >=0.99),"status"] = "FP" # here too
  }

  # correct rvsb feature
  if("RVSB" %in% features) all_mut_table[which(all_mut_table$RVSB <0.5),"RVSB"]=0.5

  # return the table used for training
  train_table_chunk = all_mut_table[which(!is.na(all_mut_table$status)),]
  if(! exists("train_table")) {train_table=train_table_chunk} else {train_table = rbind(train_table, train_table_chunk)}

  all_calls = readVcf(vcf, genome)
}

# re-equilibrate the classes, otherwise the random forest with give uncalibrated probabilities
# i.e. the the point that maximize the sens/spec will not be at 0.5 but at x=proportion of FP
train_table = train_table[c(which(train_table$status=="TP"), 
                            sample(which(train_table$status=="FP"), length(which(train_table$status=="TP")))),]

#plots
pdf(paste(output_folder,"/RF_plots.pdf",sep=""),8,8)
par(mfrow=c(2,2))
for(f in features) {
  hist(train_table[which(train_table$status=="TP"),f], main=f, xlab=f, col=adjustcolor("darkgreen",0.7))
  hist(train_table[which(train_table$status=="FP"),f],main=f, xlab=f, col=adjustcolor("red",0.7))
}

# training of the random forest
rf = randomForest(as.factor(status) ~ .,
                  data = train_table,
                  importance = TRUE, # to allow us to inspect variable importance
                  ntree = 5000, sampsize = as.numeric(table(train_table$status)["TP"]), probs=T, mtry = (length(colnames(train_table)) -1))
save(rf, file=paste(output_folder,"/RF_needlestack.Rdata",sep=""))

# compute a ROC curve from k-fold cross validation
folds <- createFolds(train_table$status, 10)
all_pred_rf = c()
all_status_rf = c()

for(i in 1:10){
  #print(paste("fold: ",i,sep=""))
  test = train_table[folds[[i]],]
  train = train_table[-folds[[i]],]

  train$old_status=train$status; test$old_status=test$status
  train$status = 0 ; train[which(train$old_status=="FP"),"status"] = 1 # FP are status 1
  test$status = 0 ; test[which(test$old_status=="FP"),"status"] = 1 # FP are status 1
  rf_fold = randomForest(as.factor(status) ~ ., data = train[,c(features,"status")], importance = TRUE,
                         ntree = 500, sampsize = min(as.numeric(table(train$status))), probs=T, mtry = (length(colnames(train)) -2))
  prediction = predict(rf_fold, test, type="prob")[,"1"] # should be sure that we take the column corresponding to the false positives
  all_pred_rf = c(all_pred_rf, prediction)
  all_status_rf = c(all_status_rf, test$status)
}

perf = performance( prediction( all_pred_rf, all_status_rf ), "rec" ,"spec") #sens in y and tdr in x
plot(perf, colorize=T, lwd=3, xlab="specificity", ylab="sensitivity", xaxt='n')
# roc corresponds to the accuracy to recognize the false positive, not true ones (so it is not an accuracy of the calling but of the classifier)
auc = performance( prediction( all_pred_rf, all_status_rf ), "auc" )@y.values[[1]]
text(0.98, 1, paste("auc=",round(auc,3)))
garbage <- dev.off()
