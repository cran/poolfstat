#' S4 class to represent a population tree or admixture graph and its underlying fitted parameter.
#'
#' @slot graph The graph in 3 column format originated from the fitted graph.params object
#' @slot dot.graph The fitted graph in dot format
#' @slot score the score of the model (squared Mahalanobis distance between the observed and fitted basis F-statistics vectors)
#' @slot bic The Bayesian Information Criterion associated to the model
#' @slot fitted.outstats a matrix containing the target values of the fstats, the fitted values and the Z-score measuring the deviation of the fitted values from the target values in units of standard errors (i.e., Z=(fitted.value-target.value)/se(target.value))
#' @slot edges.length a vector containing the estimated edges.length. Note finally, that the (two) edges coming from the roots are assumed of equal length (i.e., unrooted branch) as these are non-identifiable by the method.
#' @slot edges.length.scaled If drift.scaling=TRUE, the estimated edges.length in units of t/2N
#' @slot edges.length.ci A matrix with two columns (or four columns if drift scaled lengths are computed) containing for each edge length (in a row) the 95\% CI lower and higher bounds (columns 3 and 4 containing 95\% CI lower and higher bounds of drift scaled lengths, if any)
#' @slot admix.prop a vector containing the estimated admixture proportions (if any)
#' @slot admix.prop.ci a matrix with two columns containing for each admixture proportion (in a row) the 95\% CI lower and higher bounds
#' @slot nodes.het The estimated heterozygosities for all nodes (if available; see drift.scaling argument in fit.graph)
#' @slot fitted.f2.mat the matrix of all the fitted F2 statistics (obtained from fitted admixture graph parameter values) from which all the fitted fstats can be derived.
#' @slot optim.results list containing results of the optim call
#' @details The dot.graph element allows to plot the graph using grViz() from the DiagrammeR package or with the dot program after writing the files (e.g., dot -Tpng inputgraph.dot in terminal). Note that the dot file may be customized (e.g., to change leave color, parameter names...). 
#' @seealso To generate fitted.graph object, see \code{\link{fit.graph}}.
#' @aliases fitted.graph
gen.fitted.graph<-setClass(Class = "fitted.graph",
                   representation(graph="matrix",dot.graph="character",score="numeric",bic="numeric",fitted.outstats="matrix",edges.length="numeric",edges.length.scaled="numeric",edges.length.ci="matrix",admix.prop="numeric",admix.prop.ci="matrix",nodes.het="numeric",fitted.f2.mat="matrix",optim.results="list")
)


#' Check fitted.graph objects
#' @param x Object to be tested
#' @export
is.fitted.graph <- function (x) 
{
  res <- (is(x,"fitted.graph") & validObject(x))
  return(res)
}

#' Show fitted.graph object
#' @param object Object of class fitted.graph
setMethod("show","fitted.graph",
          function ( object ){
            cat ( " * * * fitted.graph Object * * * \n" )
            cat("Example of useful functions are compare.fitted.stats() to evaluate the fit and plot() to visualize the graph (implements an interface for grViz() from the DiagrammeR package)\n")
          }
)

#' plot pairwisefst object
#' @param x Object of class fitted.graph
#' @param y dummy argument
setMethod("plot","fitted.graph",
          function ( x, y ){
            grViz(x@dot.graph)
          }
)

