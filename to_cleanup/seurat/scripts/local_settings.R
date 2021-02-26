#PATHS
data.dir <<- "~/Documents/primes_storage/data/"
data.dir2 <<- "/Volumes/easystore/primes_storage/data/"
# data.dir <<- "/Volumes/scqc/data/"
if (! exists("data.from.pg")) {
  output.dir <<- "~/Documents/primes_storage/output/"
  source.dir.prefix <<- "~/Documents/primes_storage/results/"
  # source.dir.prefix <<- "/Volumes/scqc/output/"
} else {
  output.dir <<- "~/Documents/primes_storage/output_pg/"
  #output.dir <<- "/Volumes/scqc/output_pg/"
  source.dir.prefix <<- output.dir
}


save.res.1 <<- FALSE
