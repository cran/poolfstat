% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.graph.params.R
\docType{class}
\name{graph.params-class}
\alias{graph.params-class}
\alias{graph.params}
\title{S4 class to represent a population tree or admixture graph and its underlying parameter.}
\description{
S4 class to represent a population tree or admixture graph and its underlying parameter.
}
\details{
The graph is specified by a three column (character) matrix giving for each edge (whether admixed or not) to i) the child node; ii) the parent node; iii) the admixture proportion. For non-admixed edge, the third column must be blank. An admixed node should be referred two times as a child node with two different parent node and two different admixture proportions coded as alpha and (1-alpha) (parentheses are mandatory) if alpha is the name of the parameter for admixture proportion. The dot.graph element allows to plot the graph using grViz() from the DiagrammeR package or with the dot program after writing the files (e.g., dot -Tpng inputgraph.dot in terminal). Note that the dot file may be customized (e.g., to change leave color, parameter names...).
}
\section{Slots}{

\describe{
\item{\code{graph}}{The graph in 3 column format (see details)}

\item{\code{dot.graph}}{The graph in dot format}

\item{\code{is.admgraph}}{If FALSE the graph is binary tree (i.e., no admixture events), if TRUE the graph is an admixture graph}

\item{\code{n.leaves}}{Number of leaves of the graph}

\item{\code{leaves}}{Name of the leaves}

\item{\code{root.name}}{Name of the root}

\item{\code{n.nodes}}{Number of nodes (including root)}

\item{\code{nodes.names}}{Name of the nodes}

\item{\code{n.edges}}{Number of edges (including admixture edges)}

\item{\code{edges.names}}{Names of the edges (coded as "Parent node Name"<->"Child node Name")}

\item{\code{n.adm.nodes}}{Number of admixed nodes (=0 if is.admgraph=FALSE). This is also the number of admixed parameters since only two-ways admixture are assumed for a given node}

\item{\code{adm.params.names}}{Names of the admixed parameters}

\item{\code{graph.matrix}}{The graph incidence matrix consisting of n.leaves rows and n.edges columns. The elements of the matrix are the weights of each edge (in symbolic representation) for the different possible paths from the leaves to the graph root.}

\item{\code{root.edges.idx}}{Indexes of the graph.matrix columns associated to the (two) edges connected to the root}

\item{\code{f2.target}}{The (n.leaves-1) stats F2 involving popref (i.e., of the form F2(popref;pop))}

\item{\code{f2.target.pops}}{A matrix of (n.leaves-1) rows and 2 columns containing the names of populations of the F2 stats. The first column is by construction always popref. The order is the same as in f2.target}

\item{\code{f3.target}}{The (n.leaves-1)(n.leaves-2)/2 stats F3 involving popref as a target (i.e., of the form F3(popref;popA,popB))}

\item{\code{f3.target.pops}}{A matrix of (n.leaves-1)(n.leaves-2)/2 rows and 3 columns containing the name of popref in the first column and the names of the two populations involved in the F3 stats. The order is the same as in f3.target}

\item{\code{popref}}{The name of the reference population defining the fstats basis}

\item{\code{f.Qmat}}{A square matrix of rank n.leaves(n.leaves-1)/2 corresponding to the error covariance matrix of the F2 and F3 estimates}

\item{\code{Het}}{Estimated leave heterozygosities (if present in the fstats object)}
}}

\seealso{
To generate graph.params object, see \code{\link{generate.graph.params}}. The object may be used to estimate graph parameters with the function \code{\link{fit.graph}} or to generate files for the qpGraph software with \code{\link{graph.params2qpGraphFiles}}. See also \code{\link{graph.params2symbolic.fstats}} to obtain symbolic representation of Fstats from the matrix "Omega".
}
