\name{GenAUD_COL_MS}
\alias{GenAUD_COL_MS}
\title{Generating sequential Uniform Design of Experiments (Augmenting Factors) using diffrent initial designs}
\usage{
  GenAUD_COL_MS(xp, n, s, q, crit, maxiter, nshoot, rand_seed, vis=FALSE)
}
\description{
  This function takes n,s,q and other arguments to output a list(described below).
}
\arguments{

  \item{xp}{an integer matrix R object. Representing the previous existing design matrix.}

  \item{n}{an integer R object.}

  \item{s}{an integer R object. Factor of Experiment (including the number of factors in previous design xp).}

  \item{q}{an integer R object. Level of Experiment.}

  \item{crit}{a character R object. Type of criterion to use:

             "CD2" -- Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2" -- Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimum Coherence ;

             "A2" -- Mean Squared Correlation.
    }

  \item{maxiter}{a positive integer  R  object. Maximum iteration number in outer while loop of SATA algorithm in each shoot.}

  \item{nshoot}{Total counts to try different initial designs.}

  \item{rand_seed}{a integer R object. It controls the random seed used for optimization.}

  \item{vis}{a boolean R object. If true, plot the criterion value sequence for all shoots.}
}

\value{
  Best design over all shoots.
}

\author{
   Aijun Zhang, Zebin Yang, Haoyu Li, Shijie Quan
}

\examples{
  D1 = DesignQuery(n=18, s=3, q=3, crit="MD2")
  D2 = GenAUD_COL_MS(xp=D1, n=18, s=5, q=3, crit="MD2",
                 maxiter=100, nshoot=5)
}
