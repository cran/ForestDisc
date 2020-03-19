\name{Extract_cont_splits}
\alias{Extract_cont_splits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function: Continuous split extraction from Random Forest
}
\description{
Extraction of the splits learned by random forest regarding continuous predictors.
}
\usage{
Extract_cont_splits(SelectedTREES)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SelectedTREES}{
The output of the function RF2Selectedtrees()
}
}

\value{
List with 2 components:
\item{continuous_var }{Vector of continuous predictors.}
\item{continuous_splits }{Data frame of splits learned by random forest algorithm regarding continuous predictors.}

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

