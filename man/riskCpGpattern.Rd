\name{riskCpGpattern}
\alias{riskCpGpattern}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Plot the detected association pattern with one phenotype using heatmap
}
\description{
    The detected association pattern is a way to visualize the p-values providied by the HIRE function.
}
\usage{
riskCpGpattern(pval_matr, main_title = "Detected association pattern", hc_row_ind = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval_matr}{
	the p-value matrix for one phenotype, where one row represents a CpG site and one column indicates one cell type. 
}
  \item{main_title}{
	the title name. The default is "Detected association pattern".
}
  \item{hc_row_ind}{
	whether we conduct hierarchical clustering in the row. The default is FALSE. 
}
}
\details{
	This function depends on the heatmap.2 function in the gplots R package.  
}

\value{
	return a heatmap 
}
%\references{
% HIRE manuscript
%}
\author{
xxx
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
\dontrun{
#case vs control
#Visualize the association pattern with the case/control status in the first 100 CpG sites
riskCpGpattern(ret_list$pvalues[1:100, c(2,1,3)], 
		main_title="Detected association pattern with disease status", hc_row_ind = FALSE)
#c(2,1,3) was used because of the label switching

#age
#Visualize the association pattern with the age in the first 100 CpG sites
riskCpGpattern(ret_list$pvalues[1:100, K+c(2,1,3)], 
		main_title="Detected association pattern with age", hc_row_ind = FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ HIRE }% use one of  RShowDoc("KEYWORDS")
\keyword{ riskCpGpattern }% __ONLY ONE__ keyword per line
