\name{rfdistance}
\alias{rfdistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
RF distance between gene trees and species tree
%%  ~~function to do ... ~~
}
\description{
A function to calculate RF distance between a species tree and a set of gene trees
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rfdistance(repfolder = NULL, sfile = NULL, gfile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repfolder}{
  the path of the folder containing the s_tree.trees and g_trees.trees file
%%     ~~Describe \code{repfolder} here~~
}
  \item{sfile}{
  or the species tree file location
%%     ~~Describe \code{sfile} here~~
}
  \item{gfile}{
  with the gene tree file location
%%     ~~Describe \code{gfile} here~~
}
}
\details{
By default, this function assumes two files-- s_tree.trees and g_trees.trees--in the user specified folder path, if missing, using user supplied location for the two files. RF distances are turned in to Z scores using the \code{\link{righttailp}} function.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
return a data.frame with one column, column name rfzscore.
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
