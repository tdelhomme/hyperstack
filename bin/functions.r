## Functions to load

min_distance <-function(calls, n_samples, af_ratio=10, minQVAL=50){
  
  qvals = as.vector(t(geno(calls)$QVAL))
  chrs = rep(as.character(seqnames(rowRanges(calls,"seqnames"))), each = n_samples) 
  starts = rep(start(ranges(rowRanges(calls,"seqnames"))), each = n_samples)
  refs = rep(as.character(ref(calls)), each = n_samples)
  alts = rep(unlist(lapply(alt(calls), as.character)), each = n_samples)
  samples = rep(colnames(calls), nrow(calls))

  cur_data = all_my_data[-row_number_variant,] #remove current variant from distances computations
  cur_data = cur_data[which(cur_data$Chr==all_my_data[row_number_variant,"Chr"] & 
                              cur_data$old_SM==all_my_data[row_number_variant,"old_SM"] &
                              cur_data$AF>=af_ratio*all_my_data[row_number_variant,"AF"]) ,]
  indels = which(as.numeric(cur_data$End) > as.numeric(cur_data$Start))
  min( unlist(lapply(all_my_data[row_number_variant,"Start"]:all_my_data[row_number_variant,"End"], function(var_pos) {
    #create a vector of all position where we observed a variant
    all_pos_sort = sort(as.numeric(unique(c(cur_data$Start, unlist(lapply(indels, function(i) cur_data[i,"Start"]:cur_data[i,"End"]))))))
    min(var_pos - all_pos_sort[which(all_pos_sort<=var_pos)[length(which(all_pos_sort<=var_pos))]] ,
        all_pos_sort[which(all_pos_sort>=var_pos)[1]] - var_pos, na.rm = T) #na if only one variant compared
  })) ) #min to minimum distance over all bases of the variant 
}
