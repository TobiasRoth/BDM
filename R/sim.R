#'Calculate similarities for binary vegetation data
#'
#'This function is from the package \href{https://CRAN.R-project.org/package=simba}{simba.}
#'One of 56 (dis)similarity measures for binary data can be set to calculate
#'(dis)similarities.  The vegetational data can be in either database (list) or
#'matrix format. Same holds for the output. Simultaneous calculation of
#'geographical distances between plots and the virtual position of the
#'calculated similarity values between the parental units can be achieved if a
#'\code{data.frame} with coordinates is given.
#'
#'@param x Vegetation data, either as matrix with rows = plots and columns =
#'  species (similarities are calculated between rows!), or as
#'  \code{\link[base]{data.frame}} with first three columns representing plots,
#'  species and occurence information respectively.  All further columns are
#'  dumped before calculation.  Occurence is only considered as binary. If your
#'  list or matrix contains abundances or frequencies they are transformed
#'  automatically.
#'@param coord A \code{data.frame} with two columns containing the coordinate
#'  values of the sampling units. If given, it triggers the simultaneous
#'  calculation of the geographical distances between the sampling units, the
#'  coordinates of virtual centre-points between all possible pairs of plots,
#'  and the geographical distances in either x- or y-direction. If \code{coord}
#'  is given, output is always in \code{database} format (no matrix).
#'@param method Binary Similarity index (see \code{Details} for references and
#'  formulae), partial match to \code{"soerensen"}, \code{"jaccard"},
#'  \code{"ochiai"}, \code{"mountford"}, \code{"whittaker"}, \code{"lande"},
#'  \code{"wilsonshmida"}, \code{"cocogaston"}, \code{"magurran"},
#'  \code{"harrison"}, \code{"cody"}, \code{"williams"},    \code{"williams2"},
#'  \code{"harte"}, \code{"simpson"},    \code{"lennon"},   \code{"weiher"},
#'  \code{"ruggiero"},    \code{"lennon2"}, \code{"rout1ledge"},
#'  \code{"rout2ledge"},    \code{"rout3ledge"}, \code{"sokal1"}, \code{"dice"},
#'  \code{"kulcz1insky"}, \code{"kulcz2insky"}, \code{"mcconnagh"},
#'  \code{"manhattan"}, \code{"simplematching"}, \code{"margaleff"},
#'  \code{"pearson"}, \code{"roger"}, \code{"baroni"}, \code{"dennis"},
#'  \code{"fossum"}, \code{"gower"}, \code{"legendre"}, \code{"sokal2"},
#'  \code{"sokal3"}, \code{"sokal4"}, \code{"stiles"}, \code{"yule"},
#'  \code{"michael"}, \code{"hamann"}, \code{"forbes"}, \code{"chisquare"},
#'  \code{"peirce"}, \code{"eyraud"}, \code{"simpson2"}, \code{"legendre2"},
#'  \code{"fager"}, \code{"maarel"}, \code{"lamont"}, \code{"johnson"},
#'  \code{"sorgenfrei"}, \code{"johnson2"}. See details.
#'@param dn Neighbor definition. A geographic distance represented by a numeric
#'  or a two value vector defining a ring around each plot. Only takes effect
#'  when \code{coord} != NULL. If specified, the output does only contain
#'  similarities between neighboring plots. A plot is a neighbour to any given
#'  plot if it is within the range of the neighbor definition. See details.
#'@param normalize Logical value indicating whether the values for \code{a},
#'  \code{b} and \code{c} which are calculated in the process should be
#'  normalized to 100\% (per row, which means per plot comparison). If normalize
#'  = TRUE an asymmetric index must be chosen (see details).
#'@param listin if \code{x} is given in database (list) format this must be set
#'  to TRUE (there is no automatic detection of the format)
#'@param listout If output is wanted in database format rather than as a
#'  \code{\link[stats]{dist}}-object set this to TRUE. Output is automatically
#' @param ... 
#'  given in database-format, when \code{coord} is specified.
#'@return If listout = FALSE a distance matrix of class
#'  \code{\link[stats]{dist}} is returned. If listout = TRUE, a
#'  \code{\link[base]{data.frame}} is returned with 7 columns giving the names
#'  of the compared plots in the first two and the calculated similarity measure
#'  in the third column. The rest of the columns give the values for \code{a, b,
#'  c,} and \code{d} (in this order). Naming of the first three columns can be
#'  changed but defaults to \code{NBX} (one of the compared plots), \code{NBY}
#'  (the other one), \code{used index} (the values of the calculated index). If
#'  \code{coord} != NULL, the following columns are given in addition and the
#'  columns \code{a:d} shift to the end of the data.frame.
#' @export
#'
#' @examples
#' require(leaflet)
#' tmp <- makesq(coordID = c(645260), projection = "WGS")
#' leaflet() %>%
#'   addTiles() %>%
#'   addPolygons(data = tmp)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 


sim <- function(
    x, 
    coord=NULL, 
    method="soer", 
    dn=NULL, 
    normalize = FALSE, 
    listin = FALSE, 
    listout = FALSE, 
    ...)
{	
	if (!is.na(pmatch(method, "jaccard"))) 
        method <- "jaccard"
	METHODS <- c("soerensen", "jaccard", "ochiai", "mountford", "whittaker", "lande", "wilsonshmida", "cocogaston", "magurran", "harrison", "cody", "williams", "williams2", "harte", "simpson", "lennon", "weiher", "ruggiero", "lennon2", "rout1ledge", "rout2ledge", "rout3ledge", "sokal1", "dice", "kulcz1insky", "kulcz2insky", "mcconnagh", "manhattan", "simplematching", "margaleff", "pearson", "roger", "baroni", "dennis", "fossum", "gower", "legendre", "sokal2", "sokal3", "sokal4", "stiles", "yule", "michael", "hamann", "forbes", "chisquare", "peirce", "eyraud", "simpson2", "legendre2", "fager", "maarel", "lamont", "johnson", "sorgenfrei", "johnson2", "euclidean", "divergence")
	method <- pmatch(method, METHODS)
	if (is.na(method)){
		stop("invalid similarity method")
		}
	if (method == -1){
		stop("ambiguous similarity method")
		}
	if (listin) {
		x <- mama(x)
		x <- as.matrix(x)
		}
	x <- x > 0
	df <- as.matrix(x)
	zeina <- row.names(df)
    anz <- nrow(df)
	a <- df %*% t(df) 
	b <- df %*% (1 - t(df)) 
	c <- (1 - df) %*% t(df) 
	d <- ncol(df) - a - b - c
	if (normalize) {
		an <- a/(a+b+c)
		bn <- b/(a+b+c)
		cn <- c/(a+b+c)
		a <- an
		b <- bn
		c <- cn
		}
	if (method == 1) {
        dis <- (2*a)/((2*a) + b + c)
    }
    else if (method == 2) {
        dis <- a / (a + b + c)
    }
    else if (method == 3) {
        dis <- a / sqrt((a + b) * (a + c))
    }
    else if (method == 4) {
        dis <- (2 * a) / (a * (b + c) + (2 * b * c))
    }
    else if (method == 5) {
        dis <- ((a + b + c) / ((2 * a + b + c)/2))-1
    }
    else if (method == 6) {
        dis <- (b + c)/2
    }
    else if (method == 7) {
        dis <- (b+c)/((2*a)+b+c) 
    }
    else if (method == 8) {
        dis <- (b+c)/(a+b+c)
    }
    else if (method == 9) {
        dis <- ((2*a)+b+c)*(1-(a/(a+b+c)))
    }
    else if (method == 10) {
        dis <- pmin(b,c) / (pmax(b,c) + a)
    }
    else if (method == 11) {
        dis <- 1-((a*((2*a)+b+c))/(2*(a+b)*(a+c)))
    }
    else if (method == 12) {
        dis <- pmin(b,c) / (a+b+c) 
    }
    else if (method == 13) {
        dis <- ((b*c)+1) / ((((a+b+c)^2) - (a+b+c)) / 2)
    }
    else if (method == 14) {
        dis <- 1-((2*a) / ((2*a) + b + c)) 
    }
    else if (method == 15) {
        dis <- pmin(b,c) / (pmin(b,c) + a) 
    }
    else if (method == 16) {
        dis <- (2 * abs(b-c)) / ((2*a) + b +c)
    }
    else if (method == 17) {
        dis <- b + c
    }
    else if (method == 18) {
        dis <- a / (a+c)
    }
    else if (method == 19) {
        dis <- 1 - (log((2*a + b + c)/(a + b + c)) / log(2))
    }
    else if (method == 20) {
        dis <- (((a+b+c)^2)/(((a+b+c)^2)-(2*b*c)))-1
    }
    else if (method == 21) {
        dis <- log(2*a+b+c)-((1/(2*a+b+c))*2*a*log(2))-((1/(2*a+b+c))*((a+b)*log(a+b)+(a+c)*log(a+c))) 
    }
    else if (method == 22) {
        dis <- log(2*a+b+c)-((1/(2*a+b+c))*2*a*log(2))-((1/(2*a+b+c))*((a+b)*log(a+b)+(a+c)*log(a+c)))
    	dis <- exp(dis)-1
    }
    else if (method == 23) {
        dis <- a / (a + 2*(b + c))
    }
    else if (method == 24) {
        dis <- a / (pmin((b+a),(c+a))) 
    }
    else if (method == 25) {
        dis <- a / (b+c)
    }
    else if (method == 26) {
        dis <- ((a/2) * ((2*a) + b +c)) / ((a+b)*(a+c)) 
    }
    else if (method == 27) {
        dis <- ((a^2)-(b*c)) / ((a+b)*(a+c))
    }
    else if (method == 28) {
        dis <- (b+c) / (a+b+c+d) 
    }
    else if (method == 29) {
        dis <- (a+d) / (a+b+c+d) 
    }
    else if (method == 30) {
        dis <- (a * (a+b+c+d)) / ((a+b) * (a+c))
    }
    else if (method == 31) {
        dis <- ((a*d) - (b*c)) / sqrt((a + b)*(a + c)*(d + b)*(d + c))
    }
    else if (method == 32) {
        dis <- (a + d) / (a + 2*(b + c) +d)
    }
    else if (method == 33) {
        dis <- ((sqrt(a*d))+a) / ((sqrt(a*d))+b+c+a)
    }
    else if (method == 34) {
        dis <- ((a*d) - (b*c)) / (sqrt((a+b+c+d)*(a+b)*(a+c)))
    }
    else if (method == 35) {
        dis <- ((a+b+c+d) * (-1 * ((a/2)^2))) / ((a+b)*(a+c))
    }
    else if (method == 36) {
        dis <- (a - (b+c)+d) / (a+b+c+d)
    }
    else if (method == 37) {
        dis <- a / (a+b+c+d)
    }
    else if (method == 38) {
        dis <- (a*d) / sqrt((a+b)*(a + c)*(d + b)*(d + c))
    }
    else if (method == 39) {
        dis <- ((2*a)+(2*d)) / (a+d+(a+b+c+d))
    }
    else if (method == 40) {
        dis <- (a+d) / (b+c)
    }
    else if (method == 41) {
        dis <- log(((a+b+c+d) * (( abs((a*d)-(b*c)) - ( (a+b+c+d) / 2))^2) / ((a+b)*(a+c) *(b+d)*(c+d))))
    }
    else if (method == 42) {
        dis <- ((a*d)-(b*c)) / ((a*d)+(b*c))
    }
    else if (method == 43) {
        dis <- (4*((a*d) - (b*c))) / ((a+d)^2 + (b+c)^2)
    }
    else if (method == 44) {
        dis <- ((a+d) - (b+c)) / (a+b+c+d)
    }
    else if (method == 45) {
        dis <- (a*(a+b+c+d) - (2*max(a+b, a+c))) / (((a+b+c+d)*min(a+b, a+c)) - (2*max(a+b, a+c)))
    }
    else if (method == 46) {
        dis <- ((a+b+c+d)*((a*d) - (b*c))^2) / ((a+b)*(a+c)*(b+d)*(c+d))
    }
    else if (method == 47) {
        dis <- ((a*d) - (b*c)) / ((a+c)*(b+d))
    }
    else if (method == 48) {
        dis <- (a - ((a+b)*(a+c))) / ((a+b)*(a+c)*(b+d)*(c+d))
    }
    else if (method == 49) {
        dis <- a / (a+b)
    }
    else if (method == 50) {
        dis <- (3*a) / ((3*a) + b +c)
    }
    else if (method == 51) {
        dis <- (a / sqrt(min(a+b, a+c)*max(a+b, a+c))) - (1/(2*sqrt(min(a+b, a+c))))
    }
    else if (method == 52) {
        dis <- ((2*a) - (b+c)) / ((2*a) + b +c)
    }
    else if (method == 53) {
        dis <- a / ((2*a) + b + c)
    }
    else if (method == 54) {
        dis <- a / (2*b)
    }
    else if (method == 55) {
        dis <- a^2 / ((a+b)*(a+c))
    }
    else if (method == 56) {
        dis <- (a/(a+b)) + (a/(a+c))
    }
    else if (method == 57) {
        dis <- (sqrt(b + c) / (a+b+c+d))
    }
    else if (method == 58) {
        dis <- (sqrt(b + c) / sqrt(a+b+c+d))
    }
	dis <- as.dist(dis)
	attr(dis, "Size") <- anz
    attr(dis, "Labels") <- zeina
    attr(dis, "method") <- METHODS[method]
    attr(dis, "call") <- match.call()
    class(dis) <- "dist"
    if (listout) {
        dis <- liste(dis, entry=METHODS[method])
        dis$a <- a[row(a) > col(a)]
	    dis$b <- b[row(b) > col(b)]
	    dis$c <- c[row(c) > col(c)]
	    dis$d <- d[row(d) > col(d)]
    	}
    if (!is.null(coord)){
	   xydist <- liste(dist(coord), entry="distance")
	   dis <- cbind(xydist, as.vector(dis))
	   names(dis)[4] <- METHODS[method]
	   X <- (outer(coord[,1], coord[,1], FUN="+"))*0.5
	   Y <- (outer(coord[,2], coord[,2], FUN="+"))*0.5	   
	   dis$X <- X[row(X) > col(X)]
	   dis$Y <- Y[row(Y) > col(Y)]
	   dis$xdist <- dist(coord[,1])
	   dis$ydist <- dist(coord[,2])
	   dis$a <- a[row(a) > col(a)]
	   dis$b <- b[row(b) > col(b)]
	   dis$c <- c[row(c) > col(c)]
	   dis$d <- d[row(d) > col(d)]
	   if (!is.null(dn)) {
	       if(length(dn)==1){
	           dis <- dis[(dis$distance <= dn), ]
	       }
	       else{
	           dis <- dis[((dis$distance >= min(dn)) & (dis$distance <= max(dn))), ]
	       }
	   }
    }
    return(dis)
}





