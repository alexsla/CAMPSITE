MNND <- function(x){
  a <- x[!is.na(x)]
  if(length(a) > 1){
    a <- as.matrix(dist(a))
    diag(a) <- NA
    return(mean(apply(a, 2, min, na.rm = TRUE)))
  }
  else(return(0))
}

VNND <- function(x){
  a <- x[!is.na(x)]
  if(length(a) > 1){
    a <- as.matrix(dist(a))
    diag(a) <- NA
    return(var(apply(a, 2, min, na.rm = TRUE)))
  }
  else(return(0))
}

sumModels <- function(model){
  result1 <- model$all$lin_mat
  mod.length <- round(max(model$all$lin_mat[,4])*100, 0)
  trait_df <- list()
  for (k in 1:length(model$all$trait_mat)){
    trait_df[[k]] <- model$all$trait_mat[[k]]
    trait_df[[k]] <- trait_df[[k]][5:length(trait_df[[k]])]
    length(trait_df[[k]]) <- mod.length
    names(trait_df[[k]]) <- c(1:mod.length)
  }
  trait_df <- matrix(unlist(trait_df), ncol = mod.length, byrow = TRUE)
  result2 <- apply(trait_df, 2, var, na.rm = TRUE)
  result3 <- apply(trait_df, 2, mean, na.rm = TRUE)
  result4 <- apply(trait_df, 2, MNND)
  result5 <- apply(trait_df, 2, VNND)
  result6 <- model$gsp_extant$tips
  result7 <- model$gsp_extant$tree$Nnode
  result8 <- model$gsp_extant$tree
  result9 <- model$gsp_fossil$tree
  me <- list(
    lineages = result1,
    traits_var = result2,
    traits_mean = result3,
    MNND = result4,
    VNND = result5,
    tip_traits = result6,
    Nnode_extant = result7,
    tree_extant = result8,
    tree_fossil = result9
  )
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

tipTraitPlot <- function(tip_traits, comp, selec){
  n.obs <- sapply(tip_traits, length)
  seq.max <- seq_len(max(n.obs))
  tip_traits_df <- as_tibble(sapply(tip_traits, "[", i = seq.max))
  colnames(tip_traits_df) <- c(1:100)
  tip_traits_df <- add_column(tip_traits_df,
                              competition = comp,
                              selection = selec)
  tip_traits_df
}

sumBL <- function(x){
  max.bl <- max(distRoot(x, method = "patristic"))
  bl <- numeric()
  for (i in 0:(max.bl-1)){
    if(i < round(max.bl-1)) bl[[i+1]] <- sum(distRoot(timeSliceTree(x, i, plot = F), method = "patristic")) -
        sum(distRoot(timeSliceTree(x, i+1, plot = F), method = "patristic"))
    else bl[[i+1]] <- sum(distRoot(timeSliceTree(x, i, plot = F), method = "patristic")) + 1*(ceiling(max.bl)-max.bl)
  }
  length(bl) <- 50
  bl[is.na(bl)] <- 1
  return(data.frame(c(1:50),
                    rev(bl)))
}

extractSignal <- function(model, method = c("K", "lambda")){
  tree <- model$tree_extant
  tips <- model$tip_traits
  tips <- tips[tree$tip.label]
  if (method == "K")
  {K <- phylosig(tree, tips, "K") 
  return(as.numeric(K))} else
    if (method == "lambda")
    {lambda <- phylosig(tree, tips, "lambda")
    return(lambda$lambda)}
}

fitTraitModels <- function(model){
  models <- numeric()
  tree <- model$tree_extant
  tree$edge.length <- sapply(tree$edge.length, function(x) x + rnorm(1, sd = 0.001))
  tips <- as.matrix(model$tip_traits)
  models[[1]] <- transformPhylo.ML(tips, tree, model = "BM")$AICc
  models[[2]] <- transformPhylo.ML(tips, tree, model = "OU")$AICc
  models[[3]] <- transformPhylo.ML(tips, tree, model = "ACDC")$AICc
  models[[4]] <- fit_t_comp(tree, tips[,1], model = "MC")$aicc
  models[[5]] <- fit_t_comp(tree, tips[,1], model = "DDlin")$aicc
  models[[6]] <- fit_t_comp(tree, tips[,1], model = "DDexp")$aicc
  names(models) <- c("BM", "OU", "ACDC", "MC", "DDlin", "DDexp")
  if (names(sort(models))[[1]] == "ACDC"){
    if (transformPhylo.ML(tips, tree, model = "ACDC")$ACDC[,3] < 0) {names(models)[[3]] <- "ACDCdec"} 
    else {names(models)[[3]] <- "ACDCinc"}
  }
  return(names(sort(models))[[1]])
}

fitDivModels <- function(model){
  models <- numeric()
  tree <- model$tree_extant
  f.null <- function(t,y){0}
  f.cst <- function(t,y){y[1]}
  f.lin <- function(t,y){abs(y[1] + y[2] * t)}
  f.exp <- function(t,y){y[1] * exp(y[2] * t)}
  models[[1]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.cst, f.mu = f.null, lamb_par = c(0.09), mu_par = c(), f = 1, cst.lamb = T, fix.mu = T, cond = "stem")$aicc
  models[[2]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.lin, f.mu = f.null, lamb_par = c(0.09, 0.001), mu_par = c(), f = 1, fix.mu = T, cond = "stem")$aicc
  models[[3]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.exp, f.mu = f.null, lamb_par = c(0.05, 0.01), mu_par = c(), f = 1, expo.lamb = T, fix.mu = T, cond = "stem")$aicc
  models[[4]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.cst, f.mu = f.cst, lamb_par = c(0.09), mu_par = c(0.005), f = 1, cst.lamb = T, cst.mu = T, cond = "stem")$aicc
  models[[5]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.lin, f.mu = f.cst, lamb_par = c(0.09, 0.001), mu_par = c(0.005), f = 1, cst.mu = T, cond = "stem")$aicc
  models[[6]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.exp, f.mu = f.cst, lamb_par = c(0.05, 0.01), mu_par = c(0.005), f = 1, expo.lamb = T, cst.mu = T, cond = "stem")$aicc
  models[[7]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.cst, f.mu = f.lin, lamb_par = c(0.09), mu_par = c(0.005, 0.0001), f = 1, cst.lamb = T, cond = "stem")$aicc
  models[[8]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.lin, f.mu = f.lin, lamb_par = c(0.09, 0.001), mu_par = c(0.005, 0.0001), f = 1, cond = "stem")$aicc
  models[[9]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.exp, f.mu = f.lin, lamb_par = c(0.05, 0.01), mu_par = c(0.005, 0.0001), f = 1, expo.lamb = T, cond = "stem")$aicc
  models[[10]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.cst, f.mu = f.exp, lamb_par = c(0.09), mu_par = c(0.0035, 0.001), f = 1, cst.lamb = T, expo.mu = T, cond = "stem")$aicc
  models[[11]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.lin, f.mu = f.exp, lamb_par = c(0.09, 0.001), mu_par = c(0.0035, 0.001), f = 1, expo.mu = T, cond = "stem")$aicc
  models[[12]] <- fit_bd(phylo = tree, tot_time = max(nodeHeights(tree)), f.lamb = f.exp, f.mu = f.exp, lamb_par = c(0.05, 0.01), mu_par = c(0.0035, 0.001), f = 1, expo.lamb = T, expo.mu = T, cond = "stem")$aicc
  names(models) <- c("Bcst", "Blin", "Bexp", "BcstDcst", "BlinDcst", "BexpDcst", "BcstDlin", "BlinDlin", "BexpDlin", "BcstDexp", "BlinDexp", "BexpDexp")
  return(names(sort(models))[[1]])
}