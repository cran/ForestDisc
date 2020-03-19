\name{Select_cont_splits}
\alias{Select_cont_splits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function: Continuous cut points Selection
}
\description{
Build the optimal set of cut points C for discretization, based on moment matching.
The set of split points S extracted through Extract_cont_splits() function is mapped to a reduced set of cut points C.

}
\usage{
Select_cont_splits(cont_splits,max_splits,opt_meth)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cont_splits}{
Output of the function Extract_cont_splits().
}
  \item{max_splits}{
Maximum number of cut points allowed. Possible values range between 2 and 10.  Default value = 10.
}
  \item{opt_meth}{
The non-linear optimization algorithm to use in order to get the optimal set of cut points matching as well as possible the set of split points. The possible values are DIviding RECTangles algorithm "directL",  NelderMead Simplex method "NelderMead", Sequential Least-Squares Quadratic Programming "SLSQP". (more details about these non-linear optimization algorithms can be found in the documentation of the "NLopt" library).
}
}

\value{
List with 2 components:
\item{All_splits }{Data frame of solutions returned for k values ranging from 2 to 'max_splits'.}
\item{Selected_splits }{Data frame of the best solution returned.}
%% ...
}

\author{
Haddouchi Maïssae
}

\examples{
data(iris)
Mydata=iris
id_target=5
set.seed(1234)
X=Mydata[,1:(id_target-1)]
Y=Mydata[,id_target]
ntree=50
RFTREES=RF2Selectedtrees(X,Y,ntree)
RFCONTSPLITS=Extract_cont_splits(RFTREES)
RFSELECTCONTSPLITS=Select_cont_splits(cont_splits=RFCONTSPLITS,max_splits=10,opt_meth="NelderMead")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Optimization }% use one of  RShowDoc("KEYWORDS")
\keyword{ Discretization }% __ONLY ONE__ keyword per line
