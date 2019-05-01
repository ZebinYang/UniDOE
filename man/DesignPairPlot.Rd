\name{DesignPairPlot}
\alias{DesignPairPlot}
\title{Draw pair plot for design of experiments}
\usage{
DesignPairPlot(D,Diag)
}
\description{
  This function takes a design D and a boolean value Diag to draw pair plot.
}

\arguments{

\item{D}{a matrix object. Design of Experiment.}

\item{Diag}{a boolean R object.}

}

\value{
  A pair plot
}

\author{
   Aijun Zhang, Haoyu Li, Shijie Quan, Zebin Yang
}
\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}


\examples{
##e.g.1
n=12 #(must be multiples of q)
s=3
q=4
crit = "MD2"#(Mixture L2 criteria)
D = DesignQuery(n=n,s=s,q=q,crit="MD2")
DesignPairPlot(D)

##e.g.2
n=12 #(must be multiples of q)
s=3
q=3
crit = "MD2"#(Mixture L2 criteria)
D = DesignQuery(n=n,s=s,q=q,crit="MD2")
DesignPairPlot(D,Diag=TRUE)

}
