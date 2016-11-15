source("packages.R")

file.vec <- Sys.glob(file.path("data", "*.txt"))
five.sparse.list <- list()
for(f in file.vec){
  one <- fread(f)
  setnames(one, c("mass.num", "count"))
  one[, mass.int := as.integer(mass.num*1e4)]
  int.tab <- table(one$mass.int)
  not.one <- int.tab[int.tab != 1]
  if(length(not.one)){
    print(not.one)
    stop("non.unique mass.int values")
  }
  sample.id <- sub(".txt$", "", basename(f))
  five.sparse.list[[f]] <- data.table(sample.id, one)
}
five.sparse <- do.call(rbind, five.sparse.list)

five.wide <- dcast(five.sparse, mass.num ~ sample.id, value.var="count")
w.vec <- diff(five.wide$mass.num)/2
five.wide[, mid.before := mass.num - c(w.vec[1], w.vec)]
five.wide[, mid.after := mass.num + c(w.vec, w.vec[length(w.vec)])]
five.na <- melt(
  five.wide,
  id.vars=c("mid.before", "mass.num", "mid.after"),
  variable.name="sample.id",
  value.name="count")
five.na[, count.approx := approx(mass.num, count, mass.num)$y, by=sample.id]
five.na[, count.zero := ifelse(is.na(count), 0L, count)]
not.same <- five.na[count.zero != count.approx,]
if(nrow(not.same)){
  print(not.same)
  stop("some approximations not zero")
}

mult <- 1e5
by.sample <- split(five.na, five.na$sample.id)
for(sample.id in names(by.sample)){
  too.many.rows <- by.sample[[sample.id]]
  not.same.as.prev <- too.many.rows[c(TRUE, diff(count.zero)!=0),]
  last.after <- too.many.rows[.N, mid.after]
  pos.vec <- as.integer(c(not.same.as.prev$mid.before, last.after) * mult)
  out.dt <- data.table(
    chrom="chrNA",
    start=sprintf("%d", pos.vec[-length(pos.vec)]),
    end=sprintf("%d", pos.vec[-1]),
    count=not.same.as.prev$count)
  out.dir <- file.path("problems", sample.id)
  out.bg <- sprintf(out.dir, "coverage.bedGraph")
  dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)
  fwrite(out.dt, out.bg, sep="\t", col.names=FALSE)
  problem <- out.dt[, data.table(
    chrom="chrNA",
    start=start[1],
    end=end[.N])]
  fwrite(problem, file.path(out.dir, "problem.bed"), sep="\t", col.names=FALSE)
  problem.target(out.dir)
}
