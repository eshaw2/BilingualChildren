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
                      "pvt_mean_lapse_rt","brief_raw_initiate")
  continuous_df = data_df %>%
              select(c(contains(c("brief","flanker")),"pvt_number_of_lapses",
                       "pvt_count_falsestarts")) %>%
              select(-c("brief_raw_initiate","flanker_percenterrors_congruent"))
  
  spear_corr = cor(bind_cols(nnorm_df, continuous_df), method = "spearman")
  cont_corr = cor(continuous_df, method = "pearson")
  all_corr = spear_corr
  all_corr[rownames(cont_corr),colnames(cont_corr)] = cont_corr
  return(all_corr)
}

calc_projection = function(data_df, eigen_matrix, cg_process=NULL) {
  eigen_ordered = eigen_matrix[order(rownames(eigen_matrix)),]
  data_ordered = data_df[,order(colnames(data_df))]
  
  proj_data = as.matrix(data_ordered) %*% as.matrix(eigen_ordered)
  if (cg_process == "EF"){
    colnames(proj_data) = c("EF_Comp1","EF_Comp2","EF_Comp3")
  }else if (cg_process == "lang"){
    colnames(proj_data) = c("lang_Comp")
  }
  
  return(proj_data)
}

project_data = function(data_df,EF_Evec,lang_Evec,demo_df){
  EF_proj = calc_projection(data_df %>% select(contains(c("brief","flanker","pvt"))),
                            EF_Evec,"EF")
  lang_Evec = as.matrix(lang_Evec)
  lang_proj = calc_projection(data_df %>% select(contains(c('vocab','bpvs'))), lang_Evec, "lang")
  proj_data = bind_cols(demo_df,
                         as_tibble(EF_proj) %>% mutate_all(funs(scale)),
                         as_tibble(lang_proj) %>% mutate_all(funs(scale))
              )
  return(proj_data)
}

calc_power = function(lm_result,alpha) {
  n = 89
  r2 = summary(lm_result)$r.squared
  return(pwr.f2.test(u = length(lm_result$coefficients)-1, 
                     v = n-(length(lm_result$coefficients)-1), 
                     f2 = r2/(1-r2), 
                     sig.level = alpha))
}

rename_ef_cols = function(ef_df) {
  return( ef_df %>% mutate( BIH = brief_raw_inhibit,
                            BSM = brief_raw_self.monitor,
                            BS = brief_raw_shift,
                            BEC = brief_raw_emotional_control,
                            BIT = brief_raw_initiate,
                            BWM = brief_raw_working_memory,
                            BPO = brief_raw_plan_organise,
                            BTM = brief_raw_task_monitor,
                            BOM = brief_raw_organisation_of_materials,
                            FCE = flanker_percenterrors_congruent,
                            FIE = flanker_percenterrors_incongruent,
                            FCR = flanker_mean_rt_congruent,
                            FIR = flanker_mean_rt_incongruent,
                            PBR = pvt_mean_rt,
                            PLC = pvt_number_of_lapses,
                            PLR = pvt_mean_lapse_rt,
                            PFC = pvt_count_falsestarts)
    
  )
}

project_data_report = function(data_df,EF_Evec,lang_Evec,demo_df){
  EF_data = data_df %>%
              select(BPO,BS,BWM, BIT, BSM, BIH, BTM)
  EF_proj = calc_projection_report(EF_data, EF_Evec, "EF")
  lang_Evec = as.matrix(lang_Evec)
  lang_proj = calc_projection_report(data_df %>% select(contains(c('vocab','bpvs'))), lang_Evec, "lang")
  proj_data = bind_cols(demo_df,
                        as_tibble(EF_proj) %>% mutate_all(funs(scale)),
                        as_tibble(lang_proj) %>% mutate_all(funs(scale))
                          )
  return(proj_data)
}

calc_projection_report = function(data_df, eigen_matrix, col_name) {
  proj_data = as.matrix(data_df) %*% as.matrix(eigen_matrix)
  colnames(proj_data) = c(col_name)
  
  return(proj_data)
}
calc_power_report = function(n,params,r2,alpha) {
  return(pwr.f2.test(u = params, 
                     v = n-params-1, 
                     f2 = r2/(1-r2), 
                     sig.level = alpha))
}
