\name{DesignQuery}
\alias{DesignQuery}
\title{Evaluate design in terms of criteria}
\usage{
DesignQuery(n,s,q,crit,ShowCrit)
}
\description{
  This function takes size of desired design,criterion crit.
  If the required design exists in database, then return the design, else return NULL.
}
\arguments{

  \item{n}{an integer R object. Run of Experiment}

  \item{s}{an integer R object. Factor of Experiment.}

  \item{q}{an integer R object. Level of Experiment.}

  \item{crit}{a character R object. Type of criterion to use:

             "CD2" -- Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2" -- Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimum Coherence ;

             "A2" -- Mean Squared Correlation.
    }

\item{ShowCrit}{If TRUE, print CD2,MD2,maximin value of required design.

                Default: TRUE}

}

\value{
 Desired design
}

\author{
   Aijun Zhang, Zebin Yang, Haoyu Li, Shijie Quan
}

\examples{
  n = 9
  s = 3
  q = 3
  crit="MD2"
  D = DesignQuery(n,s,q,crit)
  D
}
