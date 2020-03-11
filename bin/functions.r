## Functions to load

to_NbVarWin <-function(all_my_data, row_number_variant, wind_length=100){
  cur_data = all_my_data[-row_number_variant,] #remove current variant from distances computations
  row_chr = unlist(strsplit(rownames(all_my_data[row_number_variant,]),":"))[1]
  row_start =  as.numeric(unlist(strsplit(unlist(strsplit(rownames(all_my_data[row_number_variant,]),":"))[2], "_"))[1])
  chrs = unlist(lapply(rownames(cur_data), function(r) unlist(strsplit(r,":"))[1]))
  starts = as.numeric( unlist(lapply( unlist(lapply(rownames(cur_data), function(r) unlist(strsplit(r,":"))[2])), function(rr) unlist(strsplit(rr,"_"))[1])) )
  
  target_data = cur_data[which(chrs == row_chr & 
                                 starts <= (row_start +  wind_length/2) &
                                 starts >= (row_start -  wind_length/2)),]
  nrow(target_data)
}