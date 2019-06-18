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

             "CD2" -- Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2" -- Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimum Coherence ;

             "A2" -- Mean Squared Correlation.
    }
}

\value{criterion value.}

\author{
   Aijun Zhang, Zebin Yang, Haoyu Li, Shijie Quan
}


\examples{

  x = matrix(c(1,1,1,2,2,2,3,3,3),nrow=3,byrow=TRUE)
  crit="MD2"
  value = DesignEval(x,crit)

}
