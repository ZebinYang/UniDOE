\name{Internal functions}
\title{UniDOE internal functions}
\alias{SATA_AUD}
\alias{SATA_UD}
\alias{SATA_AUD_COL}
\alias{CritEval}
\alias{panel.bar}
\alias{panel.scatter}
\alias{panel.heatmap}
\description{Internal core functions for constructing uniform design of experiments}
\usage{
GenUD(n,s,q,init,initX,crit,maxiter,hits_ratio,levelpermt,vis)
GenAUD(xp,n,s,q,initX,crit,maxiter,hits_ratio,levelpermt,vis)
GenAUD_COL(xp,n,s,q,init,initX,crit,maxiter,hits_ratio,levelpermt,vis)
CritEval(x,q,crit)
panel.bar(x, ...)
panel.scatter(x, y, ...)
panel.heatmap(x, y, ...)
}
\author{
  Aijun Zhang, Zebin Yang, Haoyu Li, Shijie Quan
}
\details{These functions are internal functions, not directly recalled by R users.}
\keyword{internal}
