#' S4 class to represent a population tree or admixture graph and its underlying parameter.
#'
#' @slot graph The graph in 3 column format (see details)
#' @slot dot.graph The graph in dot format
#' @slot is.admgraph If FALSE the graph is binary tree (i.e., no admixture events), if TRUE the graph is an admixture graph
#' @slot n.leaves Number of leaves of the graph
#' @slot leaves Name of the leaves
#' @slot root.name Name of the root
#' @slot n.nodes Number of nodes (including root)
#' @slot nodes.names Name of the nodes
#' @slot n.edges Number of edges (including admixture edges)
#' @slot edges.names Names of the edges (coded as "Parent node Name"<->"Child node Name")
#' @slot n.adm.nodes Number of admixed nodes (=0 if is.admgraph=FALSE). This is also the number of admixed parameters since only two-ways admixture are assumed for a given node
#' @slot adm.params.names Names of the admixed parameters
# #' @slot Omega A symbolic representation of the scaled covariance matrix of allele frequency with edges.names and adm.params.names
#' @slot graph.matrix The graph incidence matrix consisting of n.leaves rows and n.edges columns. The elements of the matrix are the weights of each edge (in symbolic representation) for the different possible paths from the leaves to the graph root.
#' @slot root.edges.idx Indexes of the graph.matrix columns associated to the (two) edges connected to the root
#' @slot f2.target The (n.leaves-1) stats F2 involving popref (i.e., of the form F2(popref;pop))
#' @slot f2.target.pops A matrix of (n.leaves-1) rows and 2 columns containing the names of populations of the F2 stats. The first column is by construction always popref. The order is the same as in f2.target
#' @slot f3.target The (n.leaves-1)(n.leaves-2)/2 stats F3 involving popref as a target (i.e., of the form F3(popref;popA,popB))
#' @slot f3.target.pops A matrix of (n.leaves-1)(n.leaves-2)/2 rows and 3 columns containing the name of popref in the first column and the names of the two populations involved in the F3 stats. The order is the same as in f3.target
#' @slot popref The name of the reference population defining the fstats basis
#' @slot f.Qmat A square matrix of rank n.leaves(n.leaves-1)/2 corresponding to the error covariance matrix of the F2 and F3 estimates
#' @slot Het Estimated leave heterozygosities (if present in the fstats object)
#' @details The graph is specified by a three column (character) matrix giving for each edge (whether admixed or not) to i) the child node; ii) the parent node; iii) the admixture proportion. For non-admixed edge, the third column must be blank. An admixed node should be referred two times as a child node with two different parent node and two different admixture proportions coded as alpha and (1-alpha) (parentheses are mandatory) if alpha is the name of the parameter for admixture proportion. The dot.graph element allows to plot the graph using grViz() from the DiagrammeR package or with the dot program after writing the files (e.g., dot -Tpng inputgraph.dot in terminal). Note that the dot file may be customized (e.g., to change leave color, parameter names...). 
#' @seealso To generate graph.params object, see \code{\link{generate.graph.params}}. The object may be used to estimate graph parameters with the function \code{\link{fit.graph}} or to generate files for the qpGraph software with \code{\link{graph.params2qpGraphFiles}}. See also \code{\link{graph.params2symbolic.fstats}} to obtain symbolic representation of Fstats from the matrix "Omega". 
graph.params<-setClass(Class = "graph.params",
                   representation(graph="matrix",dot.graph="character",is.admgraph="logical",n.leaves="numeric",leaves="character",root.name="character",n.nodes="numeric",nodes.names="character",n.edges="numeric",edges.names="character",n.adm.nodes="numeric",adm.params.names="character",graph.matrix="matrix",root.edges.idx="numeric",f2.target="numeric",f2.target.pops="matrix",f3.target="numeric",f3.target.pops="matrix",popref="character",f.Qmat="matrix",Het="numeric")
)


#' Check graph.params objects
#' @param x The name (or a path) of the graph.params objet
#' @export
is.graph.params <- function (x) 
{
  res <- (is(x,"graph.params") & validObject(x))
  return(res)
}

#' Show graph.params object
#' @param object Object of class graph.params
setMethod("show","graph.params",
          function ( object ){
            cat ( " * * * graph.params Object * * * \n" )
            cat("Example of useful functions are:\n  plot() to visualize the graph (interface for grViz() from the DiagrammeR package)\n  fit.graph() to estimate graph parameter values\n")
            cat ( " * * * * * * * * * * * * * * * * \n" )
          }
)

#' plot graph in graph.params object
#' @param x Object of class fitted.graph
#' @param y dummy argument
setMethod("plot","graph.params",
          function ( x,y ){
            grViz(x@dot.graph)
          }
)
