#' Online forum dataset
#'
#' A data set containing online forum posts from students at the University of California at Irvine, from 2004 (see Opsahl 2013). 
#' 
#' 
#' @name forum
#' @docType data
#'
#' @format A data set with a single array
#' \describe{
#'	\item{forum}{20 x 20 x 24 numeric matrix of weights. \code{NA} at \eqn{(i,j,t)} indicates that user \eqn{i} did not post to forum \eqn{j} in week \eqn{t}. }
#'}
#'
#' @details This data set contains online forum posts from students at the University of California at Irvine, from 2004 (see Opsahl 2013). The 20 most active users and the 20 forums to which these users posted the most are examined. The weights of the network are the number of characters posted to a given forum by a given user for each week. The 3-mode array \code{forum} contains the weights indexed by user, forum, and week, respectively. Data obtained June 8, 2018. See the link \url{http://opsahl.co.uk/tnet/datasets/OF_longitudinal_weightedchar.txt} for raw data. 
#'
#' @source \url{http://opsahl.co.uk/tnet/datasets/OF_longitudinal_weightedchar.txt}
#'
#' @references Opsahl, T. (2013). "Triadic closure in two-mode networks: Redefining the global and local clustering coefficients." Social Networks, 35(2), 159-167. <doi:10.1016/j.socnet.2011.07.001>
#'
#' @keywords datasets
#'
#' @examples
#' data("forum")
#' 
NULL