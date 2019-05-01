load("./data/UD_CD2.rda")
load("./data/UD_MD2.rda")

DesignQuery <- function(n,s,q,crit="CD2", ShowCrit = TRUE)
{
  if( (n != round(n)) || (s != round(s)) || (q != round(q))){ stop("Wrong types of n,s,q.")}
  else if(n%%q != 0){stop("n should be multiple of q.")}
  else if(s<=1 || n<=2 || q <=1){ stop(("The size of design should be larger than 2*2.")) }

  if (crit=="CD2") DataX = UD_CD2
  if (crit=="MD2") DataX = UD_MD2

  idx = which(DataX$n == n & DataX$s == s & DataX$q == q)
  if (length(idx) == 0) {
    warning("No design found.")
    return(NULL)
  } else{
    D = DataX$Design[idx][[1]]
    if(ShowCrit) cat("CD2 =", DesignEval(D, "CD2"),
                     "MD2 =", DesignEval(D, "MD2"),
                     "Maximin =", DesignEval(D, "maximin"), "\n")
    return (as.data.frame(D))
  }
}


GenUD_MS <- function(n, s, q, crit="CD2", maxiter=30, nshoot = 5, vis=FALSE)
{
  if( (n != round(n)) || (s != round(s)) || (q != round(q))){ stop("Wrong types of n,s,q.")}
  else if(n%%q != 0){stop("n should be multiple of q.")}
  else if(s<=1 || n<=2 || q <=1){ stop(("The size of design should be larger than 2*2.")) }

  crit_list = c()
  shoot_idx = c()
  time_list = c()
  bestcrit = 1e10
  for (i in 1:nshoot)
  {
    list0 = GenUD(n=n,s=s,q=q, crit=crit, maxiter = maxiter)
    crit_list = c(crit_list, list0$criterion_history)
    shoot_idx = c(shoot_idx, length(list0$criterion_history))
    time_list = c(time_list, list0$time_consumed)
    tmp = DesignEval(list0$final_design, crit = crit)
    if (tmp < bestcrit) {
      bestcrit = tmp
      bestdesign = list0$final_design
    }
  }
  if(vis) {
    par(mar=rep(2,4))
    plot(crit_list, type = "l", xlab = "", ylab = "")
    abline(v = which.min(crit_list), col=2, lty = 2)
    if (nshoot>1)  abline(v = cumsum(shoot_idx)[1:(nshoot-1)], col=4, lty=2)
    title(main = paste("Best value = ", bestcrit, " in ", round(sum(time_list),3), " sec",sep = ""))
  }
  return(as.data.frame(bestdesign))
}

GenAUD_MS <- function(xp, n, s, q, crit="CD2", maxiter=30, nshoot = 5, vis=FALSE)
{
  bflag_xp = FALSE
  for (i in 1:nrow(xp)) { bflag_xp = bflag_xp || (max(table(xp[,i]))>n/q)}
  if( (n != round(n)) || (s != round(s)) || (q != round(q))){ stop("Wrong types of n,s,q.")}
  else if(n%%q != 0){stop("n should be multiple of q.")}
  else if(s<=1 || n<=2 || q <=1){ stop(("The size of design should be larger than 2*2.")) }
  else if(is.matrix(xp)==FALSE){stop("Please input X0 to do the augmented searching. End of program.") }
  else if ((n <= nrow(xp)) | (s != ncol(xp)) ){
    stop("The size of the existing design matrix xp does not match the given n,s.")}
  else if ((!all(round(xp) == xp)) | (1 > min(xp)) | (q < max(xp))) {
        stop("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")}
  else if (bflag_xp) {
        stop("xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp.")}

  crit_list = c()
  shoot_idx = c()
  time_list = c()
  bestcrit = 1e10
  for (i in 1:nshoot)
  {
    list0 = GenAUD(xp=xp, n=n, s=s, q=q, crit=crit, maxiter = maxiter)
    crit_list = c(crit_list, list0$criterion_history)
    shoot_idx = c(shoot_idx, length(list0$criterion_history))
    time_list = c(time_list, list0$time_consumed)
    tmp = DesignEval(list0$final_design, crit = crit)
    if (tmp < bestcrit) {
      bestcrit = tmp
      bestdesign = list0$final_design
    }
  }
  if(vis) {
    par(mar=rep(2,4))
    plot(crit_list, type = "l", xlab = "", ylab = "")
    abline(v = which.min(crit_list), col=2, lty = 2)
    if (nshoot>1) abline(v = cumsum(shoot_idx)[1:(nshoot-1)], col=4, lty=2)
    title(main = paste("Best value = ", bestcrit, " in ", round(sum(time_list),3), " sec",sep = ""))
  }
  return(as.data.frame(bestdesign))
}

GenAUD_COL_MS <- function(xp, n, s, q, crit="CD2", maxiter=30, nshoot = 5, vis=FALSE)
{
  bflag_xp = FALSE
  for (i in 1:nrow(xp)) { bflag_xp = bflag_xp || (max(table(xp[,i]))>n/q)}
  if( (n != round(n)) || (s != round(s)) || (q != round(q))){ stop("Wrong types of n,s,q.")}
  else if(n%%q != 0){stop("n should be multiple of q.")}
  else if(s<=1 || n<=2 || q <=1){ stop(("\n The size of design should be larger than 2*2.")) }
  else if(is.matrix(xp)==FALSE){stop("Please input X0 to do the augmented searching. End of program.") }
  else if ((n != nrow(xp)) | (s <= ncol(xp)) ){
    stop("The size of the existing design matrix xp does not match the given n,s.") }
  else if ((!all(round(xp) == xp)) | (1 != min(xp)) | (q != max(xp))) {
        stop("The values of the existing design matrix x0 should be integers within: 1,2,3...,q.")}
  else if (bflag_xp) {
        stop("xp does not follow a balanced design, please increase the number of n or remove duplicated elements (per column) in xp.")}

  crit_list = c()
  shoot_idx = c()
  time_list = c()
  bestcrit = 1e10
  for (i in 1:nshoot)
  {
    list0 = GenAUD_COL(xp=xp, n=n, s=s, q=q, crit=crit, maxiter = maxiter)
    crit_list = c(crit_list, list0$criterion_history)
    shoot_idx = c(shoot_idx, length(list0$criterion_history))
    time_list = c(time_list, list0$time_consumed)
    tmp = DesignEval(list0$final_design, crit = crit)
    if (tmp < bestcrit) {
      bestcrit = tmp
      bestdesign = list0$final_design
    }
  }
  if(vis) {
    par(mar=rep(2,4))
    plot(crit_list, type = "l", xlab = "", ylab = "")
    abline(v = which.min(crit_list), col=2, lty = 2)
    if (nshoot>1) abline(v = cumsum(shoot_idx)[1:(nshoot-1)], col=4, lty=2)
    title(main = paste("Best value = ", bestcrit, " in ", round(sum(time_list),3), " sec",sep = ""))

  }
  return(as.data.frame(bestdesign))
}

DesignPairPlot <- function(D, Diag=FALSE) {
  if (Diag==TRUE) { pairs(D, lower.panel = panel.scatter, upper.panel = panel.heatmap, diag.panel = panel.bar, cex.labels=1.2) }
  else {
    pairs(D, lower.panel = panel.scatter, upper.panel = panel.heatmap, cex.labels=1.2)
  }
}

panel.bar <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  qx = table(x)
  par(usr = c(0, length(qx), 0, max(qx)*1.5))
  barplot(qx, width=1, space = 0, col=5, axes=FALSE, add = TRUE)
}

panel.scatter <- function(x, y, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  tmp = as.data.frame(table(x,y))
  points(tmp$x, tmp$y, cex=tmp$Freq/max(tmp$Freq), pch=19, col=4, xpd=TRUE)
  grid()
}

panel.heatmap <- function(x, y, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  tmp = table(x,y)
  adj = 1/2/(dim(tmp)-1)
  par(usr = c(-adj[1], 1+adj[1], -adj[2], 1+adj[2]))
  colmap = adjustcolor(terrain.colors(length(unique(tmp))), alpha.f=0.8)
  if (length(colmap)==1) colmap = adjustcolor(terrain.colors(3), alpha.f=0.8)[2]
  image(tmp, col=rev(colmap), add=TRUE)
}
