# science paper functions
calc_specificity = function(mat) {
  t(apply(mat, 1, function(x) x / mean(x)))
}

specificity_correlate = function(mat_a, mat_b, method="spearman") {
  
  mat_a_spec = calc_specificity(mat_a)
  mat_b_spec = na.omit(calc_specificity(mat_b))
  #mat_a_spec = t(apply(mat_a, 1, function(x)  calc_gene_specificity(x)))
  #mat_b_spec = na.omit(t(apply(mat_b, 1, function(x) calc_gene_specificity(x))))
  #mat_a_spec = t(apply(mat_a, 1, scale, scale=F))
  #mat_b_spec = na.omit(t(apply(mat_b, 1, scale, scale=F)))
  #colnames(mat_a_spec) = colnames(mat_a)
  #colnames(mat_b_spec) = colnames(mat_b)
  #mat_a_spec = log(mat_a_spec)
  #mat_b_spec = log(mat_b_spec)
  mat_a_spec = mat_a_spec[match(rownames(mat_b_spec), rownames(mat_a_spec)),]
  obj_cor = cor(mat_a_spec, mat_b_spec, method=method)
  obj_cor
}