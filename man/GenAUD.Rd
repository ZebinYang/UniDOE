\name{GenAUD}
\alias{GenAUD}
\title{Generating sequential Uniform Design of Experiments (Augmenting Runs)}
\usage{
GenAUD(xp,n,s,q,initX,crit,maxiter,hits_ratio,levelpermt,vis)
}
\description{
    This function takes n,s,q; a unchanged initial design and other arguments to output a list (described below).}
\arguments{

  \item{xp}{an integer matrix R object. Representing the previous existing design matrix.}

  \item{n}{an integer R object. Run of Experiment (including the previous design xp).}

  \item{s}{an integer R object. Factor of Experiment.}

  \item{q}{an integer R object. Level of Experiment.}

  \item{init}{a string vector R object. Initialization method for the run-augmented design:

              "rand" --Randomly generate initial design (default);

              "input". --User specified.
  }

  \item{initX}{a user-defined integer matrix R object. This is the user-defined initial augmentation matrix, and will be used when init="input".}

  \item{crit}{a character R object. Type of criterion to use:

             "CD2"  --Centered L2 Discrepancy (default) ;

             "WD2" -- Wrap-around L2 Discrepancy ;

             "MD2"  --Mixture L2  Discrepancy ;

             "maximin" -- Maximin Discrepancy ;

             "MC" -- Minimize Coherence ;

             "A2" -- Minimize Average Chi-Square.
    }

  \item{maxiter}{a positive integer  R  object. Maximum iteration number in outer while loop of SOAT algorithm.}

  \item{hits_ratio}{a float R object. Default value is 0.1, which is the ratio to accept changes of design in inner for loop.}

  \item{levelpermt}{a boolean R  object. It controls whether to use level permutation.}

  \item{vis}{a boolean variable. If true, plot the trace of criterion values.}
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
  #Example.
  #Set a fixed initial matrix:
  n=12 #(must be multiples of q)
  s=3
  q=4
  mat0 = matrix(c(1,1,1,2,2,2,3,3,3),ncol=3,byrow=TRUE)# nb. of columns=s
  crit = "MD2" #(Mixture L2 criteria)
  res = GenAUD(xp=mat0,n,s,q,crit=crit,maxiter=100,vis = TRUE)
}
