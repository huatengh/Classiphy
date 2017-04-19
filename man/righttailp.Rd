\name{righttailp}
\alias{righttailp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculating customed Z score
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A function to calculate customed Z score using the peak of the density curve as mean, and the left half of the distribution for calculating variance.
}
\usage{
righttailp(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
a numeric vector
}
}
\details{
Used for distributions with a long right tail. Internal function used by function \code{\link{MDCcount}} and \code{\link{rfdistance}}.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  returns a list of customed z scores
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{ Huateng Huang
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MDCcount}},\code{\link{rfdistance}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
x<-rnorm(100, mean = 0, sd = 5)
righttailp(x)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
