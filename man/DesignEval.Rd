\name{DesignEval}
\alias{DesignEval}
\title{Evaluate design in terms of criteria}
\usage{
DesignEval(x,crit)
}
\description{
  This function takes matrix X and crit to output the criterion value.
}
\arguments{

  \item{x}{an integer matrix object. Representing the design to be evaluated.}

  \item{crit}{a character R object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.
    }
}

\value{criterion value.}

\author{
   Aijun Zhang, Haoyu Li, Shijie Quan, Zebin Yang
}
\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}

\examples{

  x = matrix(c(1,1,1,2,2,2,3,3,3),nrow=3,byrow=TRUE)
  crit="MD2"
  value = DesignEval(x,crit)

}
