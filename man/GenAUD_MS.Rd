\name{GenAUD_MS}
\alias{GenAUD_MS}
\title{Generating sequential Uniform Design of Experiments (Augmenting Runs) using diffrent initial designs}
\usage{
GenAUD_MS(xp, n, s, q, crit, maxiter, nshoot, rand_seed, vis=FALSE)
}
\description{
  This function takes n,s,q and other arguments to output a list(described below).
}
\arguments{

  \item{xp}{an integer matrix R object. Representing the previous existing design matrix.}

  \item{n}{an integer R object. Run of Experiment (including the previous design xp).}


  \item{s}{an integer R object. Factor of Experiment.}

  \item{q}{an integer R object. Level of Experiment.}

  \item{crit}{a character R object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.
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
   Aijun Zhang, Haoyu Li, Shijie Quan, Zebin Yang
}

\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}


\examples{
  D1 = DesignQuery(n=18, s=7, q=3, crit="MD2")
  D2 = GenAUD_MS(xp=D1, n=21, s=7, q=3, crit="MD2",
                 maxiter=100, nshoot=5)
}
