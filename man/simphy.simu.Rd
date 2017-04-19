\name{simphy.simu}
\alias{simphy.simu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulating trees using Simphy%%  ~~function to do ... ~~
}
\description{
A wrapper function for using Simphy.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simphy.simu(simphy.executable = "simphy_lnx64", control.file, remove.individual.gtreefiles = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simphy.executable}{
  The file path to the executable simphy program
%%     ~~Describe \code{simphy.executable} here~~
}
  \item{control.file}{
  The file path to the control file for Simphy. For detailed information see \url{https://github.com/adamallo/SimPhy/wiki/Manual}
%%     ~~Describe \code{control.file} here~~
}
  \item{remove.individual.gtreefiles}{
  Simphy generate one file for each individual gene tree, this compiles all the files into one. Default is TRUE.
%%     ~~Describe \code{remove.individual.gtreefiles} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
returns a character vector, each is a folder path to one simulation replicate.
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
Huateng Huang
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
