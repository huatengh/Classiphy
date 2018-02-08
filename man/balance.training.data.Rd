\name{balance.training.data}
\alias{balance.training.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
balancing the training dataset, so that there are equal number of rows for two regimes
}
\description{
this is a short function for balancing the training data
}
\usage{
balance.training.data(x,repncolumnname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
trainning dataset as a data frame. should include a column called "model.id", and the repncolumname specified in the argument.
}
  \item{repncolumnname}{
a column in the input data frame indicating which rows (gene trees) belong to the same species tree, default to 'repn'
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For each species tree, equal number of rows are selected for each model in (model.id)
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
