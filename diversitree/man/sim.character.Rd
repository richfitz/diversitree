\name{sim.character}
\alias{sim.character}
\alias{make.sim.character}

\title{Simulate a Character Distribution on a Tree}

\description{Simulate a character distribution (state of each species)
  under some simple models of trait evolution.  Currently this does not
  return the full history (node states, and state changes) but this may
  be added in a future version.}

\usage{
sim.character(tree, pars, x0=0, model="bm", br=NULL)
make.sim.character(tree, pars, model="bm", br=NULL)
}

\arguments{
  \item{tree}{A bifurcating phylogenetic tree, in \code{ape}
    \dQuote{phylo} format.}

  \item{pars}{A set of model parameters (see Details below), as the
    order and interpretation depends on the model.}

  \item{x0}{Root state.  The default is zero, which may not make sense
    in all models.}

  \item{model}{Character string specifying which model to evolve the
    character under.  Possible values are \code{"bm"}, \code{"ou"},
    \code{"bbm"}, \code{"mk"} and \code{"meristic"}; see Details.}

  \item{br}{For cases where none of the models specifiable through the
    \code{model} argument are sufficient, you can provide your own
    function.  The function must have arguments \code{x0}, \code{t},
    which are the state at the base of a branch and the length of time
    to simulate over.  It must return a scalar state at the tip of the
    branch.  Future versions may change requirements of this function,
    especially to allow retention of character histories along
    branches.}
}

\details{This function duplicates functionality in other packages; see
  \code{sim.char} in \code{geiger} in particular.  The main difference
  here is that for continuous characters, this does not use the
  variance-covariance matrix, which can make it much faster for very
  large trees.  I believe that this approach is similar to \code{fastBM}
  in \code{phytools}.

  \itemize{
    \item \code{model="bm"}: Brownian Motion.  Takes a single
    parameter, representing the rate of diffusion (must be positive)

    \item \code{model="ou"}: Ornstein-Uhlenbeck process.  Takes a vector
    of three parameters, representing the rate of diffusion, strength of
    restoring force, and the "optimum" value.  The first two parameters
    must be non-negative, and the rate of diffusion must be positive.

    \item \code{model="bbm"}: Bounded Brownian Motion.  Takes a vector
    of three parameters (\code{s2}, \code{c}, \code{d}), representing
    the rate of diffusion, lower and upper bound, respectively.  The
    rate of diffusion must be positive.

    % TODO: What is the orientation of the matrix: \item
    \code{model="mk"}: Mk model (see \code{\link{make.mkn}}).  Takes a Q
    matrix as its argument.  The element \code{Q[i,j]} represents the
    rate of transition from state \code{i} to state \code{j}, and the
    diagonal elements must be such that \code{rowSums(Q)} is zero.

    \code{model="meristic"}: A special case of the Mk model, where the
    trait is meristic and character transitions are only possible
    between adjacent states.  There are three parameters (\code{k},
    \code{up}, \code{down}), representing the number of states, and rate
    of character change up (from state \code{i} to \code{i+1}) and down.
  }
}

% \examples{
% }

\author{Richard G. FitzJohn}
\keyword{models}
