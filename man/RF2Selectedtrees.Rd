\name{RF2Selectedtrees}
\alias{RF2Selectedtrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function: Trees extraction from Random Forest
}
\description{
Learn decision splits from random forest algorithm. The resulting model consists of a set of  trees where each tree is a collection of rules, and each rule is a combination of decision splits (pairs of variable/value(s)) defined from a root node to a terminal node.
}
\usage{
RF2Selectedtrees (X,Y,ntree,max_TreeRules = 'default',min_RuleSupport = 'default')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Descriptive attributes data frame.
}
  \item{Y}{
Target attribute (A response vector).
}
  \item{ntree}{
Number of trees to grow using Random Forest algorithm.
}
  \item{max_TreeRules}{
The maximum number of rules in each tree. It represents the maximum number of terminal nodes in each tree grown by random forest. The default value is the one set in random forest algorithm.
}
  \item{min_RuleSupport}{
The minimum support related to each rule (defined from a root node to a leef node). The support of a rule represents the size of its terminal node divided by the number of instances in the data. The default value used is the minimum size of terminal node set in random forest algorithm divided by the number of instances in the data.
}

}

\value{
List with components:
\item{ntree }{Number of trees.}
\item{list }{List of 'ntree' matrix where each one corresponds to a tree grown by random forest algorithm.  Each matrix consists of six columns and number of rows equal to the number of nodes in the tree. (more details can be found in the documentation of the function 'getTree' from “randomForest” package)}
\item{RF }{The original call to randomForest algorithm used.}
\item{xlevels }{vector of lists of size equal to the number of predictors. Each list corresponds to an attribute. In the case of categorical attribute, the categories are returned. In the case of continuous attribute, the distinct splits values performed by random Forest are returned.}
\item{continuous_var }{Vector of continuous predictors.}
\item{categorical_var }{Vector of categorical predictors.}
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ tree }% use one of  RShowDoc("KEYWORDS")

