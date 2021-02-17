### FILTER ###

# input: otu table (rows are taxa, columns are samples)
# output: filtered otu table
# extra: plots and realisations of pursued statistics

# Algorithm:
# 1: Remove empty rows/columns
# 2: Based on method run reduction algorithm
# 3: Compute statistics and plot results

# There are 4 different filters implemented based on method parameter

# a) simple
# Set in how many samples (proportion ) *treshold* taxa must be observed in at least *count* amount

# b) crude.var
# Algorithm calculates norm of crude covariance matrix of dataset as total variance measure. 
# In each step is removed taxa with least contribution to the total variance.
# Prameters: *keep* - proportion of taxa to keep; *keep.var* proportion of variance to keep.

# c) ait.var
# Same as crude.var but uses aitchison variance. It should work better,
# however a simple zero replacer must be used beforehand.
# Prameters: *keep* - proportion of taxa to keep; *var.keep* proportion of variance to keep.

# d) abundance
# As variance measure is used shannon index instead of variance compared to whole data value.


filter <- function(otu.table, method = "abundance", 
                  treshold = 0.2, count = 1,
                  keep = 0.8, var.keep = T){
  # check input
  if (!method %in% c("simple", "crude.var", "ait.var", "abundance")) {
    stop("Method must be one of: \"simple\", \"crude.var\", \"ait.var\", \"abudance\".")
  }
  
  if (treshold > 1 | treshold < 0) {
    stop("Parameter \"treshold\" must be number between 0 and 1.")
  }
  
  if (keep > 1 | keep < 0) {
    stop("Parameter \"keep\" must be number between 0 and 1.")
  }
  
  if (!is.logical(var.keep)) {
    stop("Paramater \"var.keep\" must be logical True / False.")
  }
  
  # Remove empty rows/cols
  r <- rowSums(otu.table) > 0
  s <- colSums(otu.table) > 0
  otu.table <- otu.table[r, s]
  
  # Calculate basic statistics of whole dataset
  zeros <- sum(otu.table == 0) / nrow(otu.table) / ncol(otu.table)
  crude.tot.var <- norm(crude.var(otu.table), type = "F")
  ait.tot.var <- norm(ait.var(otu.table), type = "F")
  s.i <- shannon.index(otu.table)
  D <- nrow(otu.table)
  all.stats<- c(crude.tot.var,  ait.tot.var, s.i, zeros, D)
  
  # Simple case
  if (method == "simple") {
    taxa.appearance <-
      apply(otu.table, MARGIN = 1, function(x)
        sum(x >= count))
    keep.t <- taxa.appearance >= treshold * ncol(otu.table)
    
    # filter data
    filtered<- otu.table[keep.t,]
    
    # compute new values of statistics
    zeros.f <- sum(filtered == 0) / nrow(otu.table) / ncol(otu.table)
    crude.tot.var.f <- norm(crude.var(filtered), type = "F")
    ait.tot.var.f <- norm(ait.var(filtered), type = "F")
    s.i.f <- shannon.index(filtered)
    D.f <- nrow(filtered)
    all.stats.f<- c(crude.tot.var.f,  ait.tot.var.f, s.i.f, zeros.f, D.f)
    
    # create table from them
    statistics <- rbind(all.stats, all.statsf)
    rownames(stats) <- c("Whole data", "Filtered data")
    colnames(stats) <- c("Removed row",
                         "Cruda variance lost", 
                         "Aitchison variance lost",
                         "Shannon index change", 
                         "Zeroes remaining in data",
                         "Reamining taxa")
    # print result
    print(stats)
    
    # create output
    output <- list("filtered" = filtered,
                   "statistics" = statistics)
    
    # it's done
    return(output)
    
      
  }else{
  
  # crude.var
  if (method == "crude.var") {
    info.lost <- rep(0, nrow(otu.table))
    for (i in 1:D) {
      info.lost[i] <-
        crude.tot.var - norm(crude.var(otu.table[-i, ]), type = "F")
    }
  }
  
  # ait.var
  if (method == "ait.var") {
    info.lost <- rep(0, nrow(otu.table))
    for (i in 1:D) {
      info.lost[i] <-
        ait.tot.var - norm(ait.var(otu.table[-i, ]), type = "F")
    }
  }
  
  # abundance
  if (method == "abundance") {
    info.lost <- rep(0, nrow(otu.table))
    for (i in 1:D) {
      info.lost[i] <- abs(s.i - shannon.index(otu.table[-i, ]))
    }
  }
  
  # removing one by one:
  if(method %in% c("crude.var", "ait.var", "abundance")){
    ord <- order(info.lost)
    cum.crude.var.lost <- rep(1, D)
    cum.ait.var.lost <- rep(1, D)
    cum.abundance.change <- rep(0, D)
    cum.zeros <- rep(0, D)
    for(i in 1:(D-1)){
      cum.crude.var.lost[i] <- 1- norm(crude.var(otu.table[-ord[1:i], ]), type = "F")  / crude.tot.var
      cum.ait.var.lost[i] <- 1 - norm(ait.var(otu.table[-ord[1:i], ]), type = "F") / ait.tot.var
      cum.abundance.change[i] <- shannon.index(otu.table[-ord[1:i], ]) / s.i
      cum.zeros[i] <- sum(otu.table[-ord[1:i], ] == 0) / ncol(otu.table) / (D-i)
    }
  }
  

  # stop criterions
  if(method == "crude.var" & var.keep){
    cut <- sum(cum.crude.var.lost < 1-keep)
    keep.t <- ord[cut:length(ord)]
  }
  
  if(method == "crude.var" & !var.keep){
    cut <- floor((1-keep)*D)
    keep.t <- ord[cut:length(ord)]
  }
  
  if(method == "ait.var" & var.keep){
    cut <- sum(cum.ait.var.lost < 1-keep)
    keep.t <- ord[cut:length(ord)]
  }
  
  if(method == "ait.var" & !var.keep){
    cut <- floor((1-keep)*D)
    keep.t <- ord[cut:length(ord)]
  }
  
  if(method == "abundance"){
    cut <- sum(cum.abundance.change > keep)-1
    keep.t <- ord[cut:length(ord)]
  }
  
  # filtered table
  filtered <- otu.table[keep.t, ]
  
  # plots
  dev.off()
  par(mfrow =c(2,2))
  plot(cum.crude.var.lost[1:(D-2)], type = "l",
       main = "Crude variance lost",
       xlab = "# of removed texa", ylab = "Variance lost [%]")
  abline(v=cut, col = "deeppink4", lty = 3)
  plot(cum.ait.var.lost, type = "l",
       main = "Aitchison variance lost",
       xlab = "# of removed taxa", ylab = "Variance lost [%]")
  abline(v=cut, col = "deeppink4", lty = 3)
  plot(cum.abundance.change, type = "l",
       main = "Shannon index change",
       xlab = "# of removed taxa", ylab = "Abundance change (1 is reference of full data)")
  abline(v=cut, col = "deeppink4", lty = 3)
  plot(cum.zeros, type = "l",
       main = "Zeroes remaining in data",
       xlab = "# of removed taxa", ylab = "Zeroes remaining [%]")
  abline(v=cut, col = "deeppink4", lty = 3)
  
  # table of statistics
  statistics<-cbind(ord, cum.crude.var.lost, cum.ait.var.lost, 
                    cum.abundance.change, cum.crude.var.lost, (D-1):0)
  colnames(statistics) <- c("Removed row",
                            "Crude variance lost", 
                            "Aitchison variance lost",
                            "Shannon index change", 
                            "Zeroes remaining in data",
                            "Remaining taxa")
  
  #print resuls
  print(statistics[cut, ])
  
  # create output list
  output <- list("filtered" = filtered,
                 "statistics.all" = all.stats,
                 "statistics" = statistics,
                 "cut" = cut)
  
  #it is done
  return(output)
  }
  }



## Some additional functions


# Computes aitinson covariance matrix of otu.table

ait.var <- function(otu.table){
  x <- as.matrix(otu.table)
  # init matrix
  tau <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  # simple zero replacer
  x[x == 0] <- 0.65
  # fill matrix
  for( i in 1:nrow(tau)){
    for( j in 1:ncol(tau)){
      tau[i,j] <- var(log((x[i,])/(x[j,])))
    }
  }
  # done
  return(tau)
}

crude.var <- function(otu.table){
  x <- as.matrix(otu.table)
  # simple zero replacer
  x[x == 0] <- 0.65
  x <- log(x)
  v <- cov(t(x), method = "spearman")
  # done
  return(v)
}

# Shannon index
shannon.index <- function(otu.table){
    N <- sum(otu.table)
    p <- rowSums(otu.table) / N
    return(-sum(p*log(p)))
}

# Simpson index
simpson.index <- function(otu.table){
  N <- sum(otu.table)
  p <- rowSums(otu.table) / N
  return(1/sum(p**2))
}
