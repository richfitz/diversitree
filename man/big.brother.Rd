\name{big.brother}
\alias{big.brother}

\title{Verbosely Track Function Evaluations}

\description{This is for making function evaluations more verbose;
  \code{big.brother} returns a function that prints its arguments and
  return value every iteration (or every several iterations if
  \code{interval} is specified)}

\usage{
  big.brother(f, interval=1, file="")
}

\arguments{
  \item{f}{A function}
  \item{interval}{Interval between printed evaluations}
  \item{file}{(Optional) file to log to.  By default, logs to the console.}
}

\section{Warning}{
  Using this directly on the likelihood functions may lead to unintended
  consequences as it affects the class attribute of the likelihood
  function and stop it causing the correct MLE function being
  dispatched.  Probably best to leave this alone.
}

\keyword{internal}
\author{Richard G. FitzJohn}
