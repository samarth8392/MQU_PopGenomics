### Creating an argument parser
library("optparse")

option_list = list(
  make_option(c("-f","--file"), type="character",default=NULL,help="chr file (tab)",metavar="character"),
  make_option(c("-o","--out"), type="character",default=NULL,help="Output",metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


#µ= (counts of mutated loci / sequence length) / 2t. The divergence time t
#t CI: (38 - 53 MYA)

mydata <- read.table(opt$file,header=T)

match <- sum(mydata$nmatch)
mismatch <- sum(mydata$nmismatch)
gap <- sum(mydata$cgap)

size <- unique(as.numeric(mydata$size2))

t_size <- sum(size)
tmin <- 38E6
tmax <- 53E6
test <- 46E6
tmed <- 39E6

true_mis <- mismatch-gap
mu_min <- ((true_mis/t_size)/2)/tmax
mu_max <- ((true_mis/t_size)/2)/tmin
mu_est <- ((true_mis/t_size)/2)/test
mu_med <- ((true_mis/t_size)/2)/tmed


cat(paste("Total matches =",match,"\n"), file=opt$out)
cat(paste("Total mismatches =",mismatch,"\n"), file=opt$out, append=T)
cat(paste("Total gaps =",gap,"\n"), file=opt$out,append=T)
cat(paste("Size =",t_size,"\n"), file=opt$out,append=T)
cat(paste("mu min =",mu_min,"\n"), file=opt$out,append=T)
cat(paste("mu max =",mu_max,"\n"), file=opt$out,append=T)
cat(paste("mu estimated =",mu_est,"\n"), file=opt$out,append=T)
cat(paste("mu median =",mu_med,"\n"), file=opt$out,append=T)
