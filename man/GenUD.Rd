\name{GenUD}
\alias{GenUD}
\title{Generating Uniform Design of Experiments}
\usage{
GenUD(n,s,q,init,initX,crit,maxiter,hits_ratio,levelpermt,vis)
}
\description{
  This function takes n,s,q and other arguments to output a list(described below).
}
\arguments{

  \item{n}{an integer R object. Run of Experiment.}

  \item{s}{an integer R object. Factor of Experiment.}

  \item{q}{an integer R object. Level of Experiment.}

  \item{init}{a string vector R object. Initialization method for the design:

              "rand" --Randomly generate initial design (default);

              "input". --User specified.
  }

  \item{initX}{a user-defined integer matrix R object. This is the user-defined initial design matrix, and will be used when init="input".}

  \item{crit}{a character R object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.
    }

  \item{maxiter}{a positive integer  R  object. Maximum iteration number in outer while loop of SATA algorithm.}

  \item{levelpermt}{a boolean R  object. It controls whether to use level permutation.}

  \item{hits_ratio}{a float R object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop.}

  \item{vis}{a boolean R object. If true, plot the criterion value sequence.}
}

\value{
  A list that contains Initial design matrix(initial_design),optimal design matrix(final_design), initial criterion value(initial_criterion), final criterion value(final_criterion) and criterion list(criterion_history) in update process.
}

\author{
   Aijun Zhang, Haoyu Li, Shijie Quan, Zebin Yang
}
\references{
 Zhang, A. and Li, H. (2017). UniDOE: An R package for constructing uniform design of experiments via stochastic and adaptive threshold accepting algorithm. Technical Report.
}


\examples{
  ##Example 1
  n=12 #(must be multiples of q)
  s=3
  q=4
  crit = "CD2"#(Centered L2 criteria)
  res = GenUD(n,s,q,crit=crit,maxiter=100)

  ##Example 2
  n=10 #(must be multiples of q)
  s=3
  q=5
  init = "rand"
  crit = "MD2" #(Mixture L2 criteria)
  vis=TRUE
  res=GenUD(n,s,q,init=init,crit=crit,maxiter=100,vis=vis)

  ##Example 3
  #If init="input", algorithm will search for better a better design with same size as initX.
  init = "input"
  n=6
  s=2
  q=3
  initX = matrix(c(1,1,2,2,3,3,3,3,1,1,2,2),ncol=2)
  res = GenUD(n,s,q,init=init,initX = initX,maxiter=100)
}
