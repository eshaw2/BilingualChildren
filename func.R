find_max_corr = function(data_df, poly_corr=NA) {
  comb_flag = F
  if (is.na(poly_corr)){
    poly_corr = as.data.frame(matrix(0, 
                                     ncol = ncol(data_df), 
                                     nrow = ncol(data_df)))
  }
  
  pearson_corr = cor(data_df, method = "pearson")
  spearman_corr = cor(data_df, method = "spearman")
  
  if(ncol(data_df)!=ncol(poly_corr)){
    pearson_corr = pearson_corr[rownames(poly_corr),
                                colnames(poly_corr)]
    spearman_corr = spearman_corr[rownames(poly_corr),
                                  colnames(poly_corr)]
  }
  
  max_corr = which_corr = as.data.frame(matrix(0, 
                                               ncol = ncol(pearson_corr), 
                                               nrow = nrow(pearson_corr)))
  for (i in 1:nrow(pearson_corr)){
    for (j in 1:ncol(pearson_corr)) {
      corrs = c(pearson_corr[i,j], spearman_corr[i,j], poly_corr[i,j])
      max_ind = which.max(abs(corrs))
      which_corr[i,j] = max_ind
      max_corr[i,j] = max(abs(corrs))
    }
  }
  return(list(corr=max_corr, type_ind=which_corr))
}

data_max_corr = function(data_df) {
  ordinal_data = data_df %>% 
                  na.omit() %>%
                  select(c(contains(c("brief")),"pvt_number_of_lapses",
                           "pvt_count_falsestarts"))
  
  continuous_data = data_df %>% 
                      na.omit() %>%
                      select(c(contains(c("flanker")),"pvt_mean_rt",
                               "pvt_mean_lapse_rt"))
  
  ordinal_psych_output = psych::polychoric(ordinal_data)
  ordinal_max_corr = find_max_corr(ordinal_data, ordinal_psych_output$rho)
  colnames(ordinal_max_corr$corr) = colnames(ordinal_data)
  rownames(ordinal_max_corr$corr) = colnames(ordinal_data)
  colnames(ordinal_max_corr$type_ind) = colnames(ordinal_data)
  rownames(ordinal_max_corr$type_ind) = colnames(ordinal_data)
  
  cont_max_corr = find_max_corr(continuous_data)
  rownames(cont_max_corr$corr) = colnames(continuous_data)
  colnames(cont_max_corr$corr) = colnames(continuous_data)
  rownames(cont_max_corr$type_ind) = colnames(continuous_data)
  colnames(cont_max_corr$type_ind) = colnames(continuous_data)
  
  comb_psych_output = psych::polyserial(continuous_data,ordinal_data)
  cross_max_corr = find_max_corr(EF_numeric_cln %>% na.omit(), comb_psych_output)
  rownames(cross_max_corr$corr) = rownames(comb_psych_output)
  colnames(cross_max_corr$corr) = colnames(comb_psych_output)
  rownames(cross_max_corr$type_ind) = rownames(comb_psych_output)
  colnames(cross_max_corr$type_ind) = colnames(comb_psych_output)
  
  mix_cont_corr = bind_rows(cont_max_corr$corr, cross_max_corr$corr)
  mix_ord_corr = bind_cols(ordinal_max_corr$corr,cross_max_corr$corr)
  all_EF_corr = bind_rows(as.data.frame(t(mix_cont_corr)),mix_ord_corr)
  
  mix_cont_ind = bind_rows(cont_max_corr$type_ind, cross_max_corr$type_ind)
  mix_ord_ind = bind_cols(ordinal_max_corr$type_ind,cross_max_corr$type_ind)
  all_EF_ind = bind_rows(as.data.frame(t(mix_cont_ind)),mix_ord_ind)
  
  return(list(corr=all_EF_corr, type_ind=all_EF_ind))
}

data_cons_corr = function(data_df) {
  nnorm_df = data_df %>% 
              select("flanker_percenterrors_congruent","pvt_mean_rt",
                      "pvt_mean_lapse_rt","pvt_count_falsestarts","brief_raw_initiate")
  continuous_df = data_df %>%
              select(c(contains(c("brief","flanker")),"pvt_number_of_lapses")) %>%
              select(-c("brief_raw_initiate","flanker_percenterrors_congruent"))
  
  spear_corr = cor(bind_cols(nnorm_df, continuous_df), method = "spearman")
  cont_corr = cor(continuous_df, method = "pearson")
  all_corr = spear_corr
  all_corr[rownames(cont_corr),colnames(cont_corr)] = cont_corr
  return(all_corr)
}