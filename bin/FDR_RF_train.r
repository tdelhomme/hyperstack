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
      --minQVAL                   - filtering calls with QVAL<minQVAL (default=50)
      --features                  - vcf features to train the model 
                                    (default=QVAL,AO,AF,DP,ERR,QUAL,RVSB,FS)

      example: Rscript FDR_RF_train.r --vcf=myvcf.bgz \n\n")
  
  q(save="no")
}

if(is.null(args$vcf)) {stop("Option --vcf should be provided")} else{vcf=args$vcf}
if(is.null(args$genome)) {genome="hg18"} else {genome=args$genome}
if(is.null(args$output_folder)) { output_folder="."} else {output_folder = args$output_folder}
system(paste("mkdir -p",output_folder,sep=" "))
out_vcf = paste(output_folder, "/", paste( sub(".vcf.gz", "", sub('.vcf.bgz', '', basename(vcf))), "RF_needlestack.vcf", sep="_"), sep="")
if(is.null(args$minQVAL)) {minQVAL=50} else {minQVAL=args$minQVAL}
if(is.null(args$features)) {features=c("QVAL","AO","AF","DP","ERR","QUAL","RVSB","FS")} else {features=as.character(unlist(strsplit(args$features,",")))}

suppressMessages(library(VariantAnnotation))
suppressMessages(library(randomForest))
suppressMessages(library(caret))
suppressMessages(library(ROCR))

all_calls = readVcf( open(VcfFile(vcf, yieldSize=1000000)), genome)
all_calls = all_calls[which(apply(geno(all_calls, "QVAL"), 1, max) >= minQVAL), ] # here can't use QUAL if the vcf was separated into different pieces (QVAL is not recalculated)

### populate the table of all mutations with features and ethnicities ###

n_samples = length(samples(header(all_calls)))

kept_variants = which(as.vector(t(geno(all_calls)$QVAL)) >= minQVAL)
all_mut_table = data.frame(matrix(NA, nrow = length(kept_variants), ncol = length(features)))
colnames(all_mut_table) = features
rownames(all_mut_table) = paste(rep(rownames(all_calls), each=n_samples), rep(samples(header(all_calls)),n_samples), sep="\\")[kept_variants]

sm_ethn = rep(info(all_calls)$ExAC_NFE, each=n_samples)[kept_variants]

other_ethn = rep ( apply(data.frame(info(all_calls)$ExAC_AMR, info(all_calls)$ExAC_EAS, info(all_calls)$ExAC_FIN, 
                                    info(all_calls)$ExAC_OTH, info(all_calls)$ExAC_SAS), 1, function(r){if(sum(!is.na(r))==0) {NA} else {max(r, na.rm=T)}}),
                   each=n_samples)[kept_variants]

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

# assign status
all_mut_table$status = NA # status as NA is for variants not used in the training
all_mut_table[which(all_mut_table$QVAL>=50 & all_mut_table$AO>=10 & all_mut_table$DP>=10 & sm_ethn>=0.2),"status"] = "TP"
all_mut_table[which(all_mut_table$QVAL>=50 & all_mut_table$AO>=10 & all_mut_table$DP>=10 & 
                      ((sm_ethn<0.2 | is.na(sm_ethn)) & other_ethn >=0.2)),"status"] = "FP"

# correct rvsb feature
if("RVSB" %in% features) all_mut_table[which(all_mut_table$RVSB <0.5),"RVSB"]=0.5

train_table = all_mut_table[which(!is.na(all_mut_table$status)),]

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
                  ntree = 5000, sampsize = as.numeric(table(train_table$status)["TP"]), probs=T)
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
  train$status = 0 ; train[which(train$old_status=="TP"),"status"] = 1
  test$status = 0 ; test[which(test$old_status=="TP"),"status"] = 1
  rf_fold = randomForest(as.factor(status) ~ ., data = train[,c(features,"status")], importance = TRUE,
                         ntree = 500, sampsize = as.numeric(table(train$status)["1"]), probs=T) 
  prediction = predict(rf_fold, test, type="prob")[,2]
  all_pred_rf = c(all_pred_rf, prediction)
  all_status_rf = c(all_status_rf, test$status)
}

perf = performance( prediction( all_pred_rf, all_status_rf ), "rec" ,"ppv") #sens in y and tdr in x
plot(perf, colorize=T, lwd=3, xlab="1-FDR", ylab="sensitivity", xaxt='n')
auc = performance( prediction( all_pred_rf, all_status_rf ), "auc" )@y.values[[1]]
text(0.98, 1, paste("auc=",round(auc,3)))
garbage <- dev.off()

# apply the RF to all the mutations
all_mut_table$FPprob = predict(rf, all_mut_table, type="prob")[,2]
FPRF = rep(NA, length(rep(rownames(all_calls), each=n_samples)))
FPRF[kept_variants] = all_mut_table$FPprob

# annotate the VCF with FPRF statistic
geno(header(all_calls))["FPRF",]=list("A","Float","False Probability of being a variant for a trained Random-Forest model")
geno(all_calls)$FPRF = matrix(data=FPRF, nrow=nrow(all_calls), byrow = T)

# write out the annotated VCF file
con = file(out_vcf, open = "a")
writeVcf(all_calls, con)
close(con)
