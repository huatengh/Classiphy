\name{MDCcount}
\alias{MDCcount}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimal deep coalescence count on the species tree
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MDCcount(repfolder = NULL, sfile = NULL, gfile = NULL, phylonet.path, save.middle.file = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{repfolder}{
%%     ~~Describe \code{repfolder} here~~
}
  \item{sfile}{
%%     ~~Describe \code{sfile} here~~
}
  \item{gfile}{
%%     ~~Describe \code{gfile} here~~
}
  \item{phylonet.path}{
%%     ~~Describe \code{phylonet.path} here~~
}
  \item{save.middle.file}{
%%     ~~Describe \code{save.middle.file} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
function (repfolder = NULL, sfile = NULL, gfile = NULL, phylonet.path,
    save.middle.file = T)
{
    require(ape)
    require(phangorn)
    require(phytools)
    require(stringr)
    outfile = paste(repfolder, "/phylonetinput.nexus", sep = "")
    resultfile = paste(repfolder, "/phylonetoutput.txt", sep = "")
    matrixfile = paste(repfolder, "/MDCcounts.txt", sep = "")
    if (is.na(repfolder)) {
        s <- read.tree(sfile)
        g <- read.tree(gfile)
        repfolder = sub(basename(gfile), "", gfile, perl = T)
    }
    else {
        s <- read.tree(paste(repfolder, "/s_tree.trees", sep = ""))
        s$edge.length <- NULL
        g <- read.tree(paste(repfolder, "/g_trees.trees", sep = ""))
        for (i in 1:length(g)) {
            g[[i]]$tip.label <- sub("_0_0$", "", g[[i]]$tip.label,
                perl = T)
            g[[i]]$edge.length <- NULL
            g[[i]] <- drop.tip(g[[i]], tip = "outgroup")
        }
    }
    subs <- cut.phy.branch(s)
    resultm <- matrix(0, ncol = length(subs[[1]]) + 1, nrow = length(g))
    for (gg in 1:length(g)) {
        sink(outfile)
        cat("#NEXUS\n\nBEGIN TREES;\n\n", sep = "")
        splist <- write.phylonet.subtrees(s, subs[[1]], "speciesTree")
        glist <- write.phylonet.subtrees(g[[gg]], subs[[1]],
            "geneTree")
        cat("END;\n\nBEGIN PHYLONET;\n\n", sep = "")
        for (i in 1:length(splist)) {
            cat("DeepCoalCount_tree {", splist[i], "} (", glist[i],
                ");\n", sep = "")
        }
        cat("END;\n", sep = "")
        sink()
        system(command = paste("java -jar ", phylonet.path, " ",
            outfile, " >", resultfile, sep = ""), wait = T)
        con = file(resultfile, "r")
        while (TRUE) {
            line = readLines(con, n = 1)
            if (length(line) == 0) {
                break
            }
            else if (length(grep("DeepCoalCount", line)) > 0) {
                line = sub(".*geneTree(\\d+).*", "\\1", line)
                x = as.integer(line)
            }
            else if (length(grep("Total number of extra lineages",
                line)) > 0) {
                resultm[gg, x + 1] = as.numeric(sub("\\D+(\\d+.*\\d*)",
                  "\\1", line))
            }
        }
        close(con)
    }
    resultm[, 2:dim(resultm)[2]] <- resultm[, 1] - resultm[,
        2:dim(resultm)[2]]
    bresult <- matrix(0, ncol = length(subs[[1]]) + 1, nrow = length(g))
    bresult[, 1] <- resultm[, 1]
    for (i in 1:length(subs[[2]])) {
        if (is.null(subs[[2]][[i]])) {
            bresult[, i + 1] <- resultm[, i + 1]
        }
        else {
            bresult[, i + 1] = resultm[, i + 1] - resultm[, subs[[2]][[i]][1] +
                1] - resultm[, subs[[2]][[i]][2] + 1]
        }
    }
    bresult[bresult < 0] <- 0
    unlink(outfile)
    unlink(resultfile)
    if (save.middle.file) {
        write.table(bresult, file = matrixfile, quote = F, sep = "\t",
            row.names = F, col.names = F)
    }
    xx <- bresult[, -1]
    Q95 <- apply(xx, MARGIN = 2, function(i) quantile(i, 0.95))
    Q05 <- apply(xx, MARGIN = 2, function(i) quantile(i, 0.05))
    IQR <- Q95 - Q05
    IQR[IQR == 0] <- 0.5
    k <- matrix(0, nrow = dim(xx)[1], ncol = dim(xx)[2])
    for (i in 1:dim(xx)[1]) {
        l <- xx[i, ] - Q95
        l[l < 0] = 0
        k[i, ] <- as.vector(as.numeric(l/IQR))
    }
    MDCbscore <- apply(k + 1, MARGIN = 1, prod)
    MDCzscore <- righttailp(bresult[, 1])
    MDCresult <- data.frame(MDCzscore = MDCzscore, MDCbscore = MDCbscore,
        stringsAsFactors = F)
    return(MDCresult)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
