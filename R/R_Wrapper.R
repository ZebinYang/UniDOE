
DesignEval<-function(x = matrix(0), crit="CD2")
{
  x = as.matrix(x)
  nlevel = max(x) - min(x) + 1
  return(CritEval(x,nlevel,crit))
}

GenUD <- function(n,s,q,init="rand",initX=matrix(0),crit="CD2",
                  maxiter=10000,hits_ratio = 0.1,levelpermt=FALSE,rand_seed=0,vis=FALSE)
{
  #check the arguments
  if( (n != round(n)) || (s != round(s)) || (q != round(q))){ stop("Wrong types of n,s,q.")}
  else if(n%%q != 0){stop("n should be multiple of q.")}
  else if(n > q^s){ stop(("n should not be larger than q^s")) }
  else if(s<1 || n<2 || q <2){ stop(("Invalid design table.")) }
  else if(init=="input"){
    initX = as.matrix(initX)
    bflag_init = FALSE
    for (i in 1:nrow(xp)) { bflag_init = bflag_init || (max(table(initX[,i]))>n/q)}
    if ((n != nrow(initX)) || (s != ncol(initX))) {
      stop("The size of the input design matrix does not match the given n,s.")
     } else if ((1 != min(initX)) | (q != max(initX)) || (!all(round(initX) == initX))){
      stop("The values of the input design matrix should be integers within: 1,2,3...,q.")
     } else if (bflag_init){
      stop("initX does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in initX.")
     }

  }

  list <- SATA_UD(n,s,q,init,initX,crit,maxiter,hits_ratio,levelpermt,rand_seed)
  names(list) = c("initial_design","final_design","initial_criterion",
                  "final_criterion","time_consumed","criterion_history")
  if(vis == TRUE){
    plot(list$criterion_history,type="l")
    bst_score = round(list$final_criterion,5)
    min_index = which.min(list$criterion_history)[1]
    abline(v = min_index,col=2)
    abline(h =list$criterion_history[min_index],col=4 )
    title(main = paste("Best value = ", bst_score, " in ", round(list$time_consumed,3), " sec",sep = ""))
  }
  return (list)
}

GenAUD <- function (xp,n,s,q,init="rand",initX=matrix(0),crit="CD2",
                    maxiter=10000, hits_ratio = 0.1,levelpermt=FALSE,rand_seed=0,vis=FALSE)
{
  #check the arguments
  xp = as.matrix(xp)
  bflag_xp = FALSE
  for (i in 1:ncol(xp)) { bflag_xp = bflag_xp || (max(table(xp[,i]))>n/q)}
  if( (n != round(n)) || (s != round(s)) || (q != round(q))){ stop("Wrong types of n,s,q.")}
  else if(n%%q != 0){stop("n should be multiple of q.")}
  else if(n > q^s){ stop(("n should not be larger than q^s")) }
  else if(s<1 || n<2 || q <2){ stop(("Invalid design table.")) }
  else if(is.matrix(xp)==FALSE){stop("Please input X0 to do the augmented searching. End of program.") }
  else if ((n <= nrow(xp)) | (s != ncol(xp)) ){
    stop("The size of the existing design matrix xp does not match the given n,s.")}
  else if ((!all(round(xp) == xp)) | (1 > min(xp)) | (q < max(xp))) {
        stop("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")}
  else if (bflag_xp) {
        stop("xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp.")}
  else if (init=="input"){
        initX = as.matrix(initX)
        bflag_init = FALSE
        for (i in 1:nrow(xp)) { bflag_init = bflag_init || (max(table(initX[,i]))>n/q)}
        if ((n <= nrow(initX)) || (s != ncol(initX))) {
          stop("The size of the input design matrix does not match the given n,s.")
         } else if ((1 > min(initX)) | (q < max(initX)) || (!all(round(initX) == initX))){
          stop("The values of the input design matrix should be integers within: 1,2,3...,q.")
         } else if (bflag_init){
          stop("initX does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in initX.")
         }
  }

  np = nrow(xp)
  nnew = n - np

  list= SATA_AUD(xp, nnew, s, q, init, initX, crit, maxiter, hits_ratio,levelpermt,rand_seed)
  names(list) = c("initial_design","final_design","initial_criterion",
                  "final_criterion","time_consumed","criterion_history")
  if(vis == TRUE){
    plot(list$criterion_history,type="l")
    bst_score = round(list$final_criterion,5)
    min_index = which.min(list$criterion_history)[1]
    abline(v = min_index,col=2)
    abline(h =list$criterion_history[min_index],col=4 )
    title(main = paste("Best value = ", bst_score, " in ", round(list$time_consumed,3), " sec",sep = ""))
  }
  return(list)
}

GenAUD_COL <- function (xp,n,s,q,init="rand",initX=matrix(0),crit="CD2",
                        maxiter=10000, hits_ratio = 0.1, levelpermt=FALSE,rand_seed,vis=FALSE)
{
  #check the arguments
  xp = as.matrix(xp)
  bflag_xp = FALSE
  for (i in 1:ncol(xp)) { bflag_xp = bflag_xp || (max(table(xp[,i]))>n/q)}
  if( (n != round(n)) || (s != round(s)) || (q != round(q))){ stop("Wrong types of n,s,q.")}
  else if(n%%q != 0){stop("n should be multiple of q.")}
  else if(n > q^s){ stop(("n should not be larger than q^s")) }
  else if(s<1 || n<2 || q <2){ stop(("Invalid design table.")) }
  else if(is.matrix(xp)==FALSE){stop("Please input X0 to do the augmented searching. End of program.") }
  else if ((n != nrow(xp)) | (s <= ncol(xp)) ){
    stop("The size of the existing design matrix xp does not match the given n,s.") }
  else if ((!all(round(xp) == xp)) | (1 != min(xp)) | (q != max(xp))) {
        stop("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")}
  else if (bflag_xp) {
        stop("xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp.")}
  else if (init=="input"){
        initX = as.matrix(initX)
        bflag_init = FALSE
        for (i in 1:nrow(xp)) { bflag_init = bflag_init || (max(table(initX[,i]))>n/q)}
        if ((n != nrow(initX)) || (s <= ncol(initX))) {
          stop("The size of the input design matrix does not match the given n,s.")
         } else if ((1 != min(initX)) | (q != max(initX)) || (!all(round(initX) == initX))){
          stop("The values of the input design matrix should be integers within: 1,2,3...,q.")
         } else if (bflag_init){
          stop("initX does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in initX.")
         }
  }

  nvp = ncol(xp)
  list= SATA_AUD_COL(xp, s - nvp, q, init, initX, crit, maxiter, hits_ratio, levelpermt,rand_seed)
  names(list) = c("initial_design","final_design","initial_criterion",
                  "final_criterion","time_consumed","criterion_history")
  if(vis == TRUE){
    plot(list$riterion_history,type="l")
    bst_score = round(list$final_criterion,5)
    min_index = which.min(list$criterion_history)[1]
    abline(v = min_index,col=2)
    abline(h =list$criterion_history[min_index],col=4 )
    title(main = paste("Best value = ", bst_score, " in ", round(list$time_consumed,3), " sec",sep = ""))
  }
  return(list)
}
