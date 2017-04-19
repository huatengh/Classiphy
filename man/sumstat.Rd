\name{sumstat}
\alias{sumstat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculat summary statistics for a set of gene tree and a species tree
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sumstat(repfolder=NULL,sfile=NULL,gfile=NULL,phylonet.path,max.freq=0.1,c=T,tscore=T,mdcscore=T,subtree=T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repfolder}{
   The path to the folder containing s_tree.trees, l_trees.trees, and g_trees.trees. If missing, taking following two file paths
%%     ~~Describe \code{x} here~~
}
  \item{sfile}{species tree file}
  \item{gfile}{gene tree file, labels need to be the same as the species tree file.}
  \item{phylonet.path}{path to the PhyloNet program}
  \item{max.freq}{maximum subtree frequency. As HGT events generate rare clades on gene tree, the more common a clade is among gene trees, the less interesting it is. Not an important paramter, just to save computation time and space. Default to 10%}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a data frame ready to be classified using DAPC
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
