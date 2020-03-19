\name{ForestDisc}
\alias{ForestDisc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate discretization for supervised learning using Random Forest and moment matching optimization
}
\description{
ForestDisc is a supervised, multivariate and non-parametric discretization algorithm based on tree ensembles learning and moment matching optimization. This version of the algorithm relies on random forest algorithm to learn a large set of split points that conserves the relationship between attributes and the target class, and on moment matching optimization to transform this set into a reduced number of cut points matching as well as possible statistical properties of the initial set of split points. For each attribute to be discretized, the set S of its related split points extracted through random forest is mapped to a reduced set C of cut points of size k.
}
\usage{
ForestDisc(data,id_target,ntree=50,max_splits=10,opt_meth="NelderMead")
}

\arguments{
  \item{data}{
Data frame to be discretized.
}
\item{id_target}{
Column id of the target class.
}
\item{ntree}{
Number of trees to grow using random forest algorithm in order to learn split points. The default value is 50.
}
\item{max_splits}{
Maximum number of cut points to be used for discretizing continuous attributes in the data. Possible values for 'max_splits' range between 2 and 10.  Default value = 10.
}
\item{opt_meth}{
The non-linear optimization algorithm to use in order to get the optimal set of cut points matching as well as possible the set of split points. The possible values are DIviding RECTangles algorithm "directL",  NelderMead Simplex method "NelderMead", Sequential Least-Squares Quadratic Programming "SLSQP". (more details about these non-linear optimization algorithms can be found in the documentation of the "NLopt" library). The default value used is "NelderMead".
}
}


\value{
List with components:
\item{Data_disc }{Discretized data.}
\item{cont_variables }{Continuous attributes column ids. }
\item{Listcutp }{List of cut points used to discretize continuous attributes.}
\item{cut_points }{Data frame summarizing the best solution returned.}
\item{opt_results }{Data frame summarizing all the solutions returned for different realizations. Each realization is determined by a size of the set of cut points, ranging between 2 and 'max_splits'.}
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
Mydata_Disc=ForestDisc(Mydata,id_target)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip }% use one of  RShowDoc("KEYWORDS")
\keyword{nonparametric }% __ONLY ONE__ keyword per line
\keyword{multivariate }% __ONLY ONE__ keyword per line
\keyword{supervised }% __ONLY ONE__ keyword per line
\keyword{trees }% __ONLY ONE__ keyword per line
\keyword{Optimization}
\keyword{Discretization}
\keyword{random forest}
\keyword{bining}
