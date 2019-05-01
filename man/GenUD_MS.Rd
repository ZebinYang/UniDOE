\name{GenUD_MS}
\alias{GenUD_MS}
\title{Generating Uniform Design of Experiments using diffrent initial designs}
\usage{
GenUD_MS(n, s, q, crit, maxiter, nshoot, vis)
}
\description{
This function takes n,s,q,crit and nshoot to return a design. nshoot number of random initial designs are used in each shoot. The design returned is the best design over all shoots.
}
\arguments{

  \item{n}{an integer R object. Run of Experiment}

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
D = GenUD_MS(36, 4, 6, crit="CD2",
             maxiter=50, nshoot = 6,
             vis=TRUE)
}
