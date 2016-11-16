source("packages.R")
library(coseg)

file.vec <- Sys.glob(file.path("data", "*.txt"))
mult <- 1e5
five.sparse.list <- list()
for(f in file.vec){
  one <- fread(f)
  setnames(one, c("mass.num", "count"))
  one[, mass.int := as.integer(mass.num*mult)]
  int.tab <- table(one$mass.int)
  not.one <- int.tab[int.tab != 1]
  if(length(not.one)){
    print(not.one)
    stop("non.unique mass.int values")
  }
  sample.id <- sub(".txt$", "", basename(f))
  five.sparse.list[[sample.id]] <- data.table(sample.id, one)
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

by.sample <- split(five.na, five.na$sample.id)
mass.spec <- list()
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
  out.df <- data.frame(
    chromStart=pos.vec[-length(pos.vec)],
    chromEnd=pos.vec[-1],
    count=as.integer(not.same.as.prev$count))
  mass.spec[[sample.id]] <- out.df
  ##fit <- PeakSegFPOPchrom(out.df, 0)
  out.dir <- file.path("problems", sample.id)
  out.bg <- file.path(out.dir, "coverage.bedGraph")
  dir.create(out.dir, showWarnings=FALSE, recursive=TRUE)
  fwrite(out.dt, out.bg, sep="\t", col.names=FALSE)
  problem <- out.dt[, data.table(
    chrom="chrNA",
    start=start[1],
    end=end[.N])]
  fwrite(problem, file.path(out.dir, "problem.bed"), sep="\t", col.names=FALSE)
  t.info <- problem.target(out.dir)

  (result <- problem.PeakSegFPOP(out.dir, "31750"))

  pen.feasible <- t.info$models[status=="feasible", paste(min(penalty))]
  result <- problem.PeakSegFPOP(out.dir, pen.feasible)
  peak.dt <- result$segments[status=="peak",]
  lo.lim <- 7100000
  up.lim <- 7200000
  some.out <- out.dt[lo.lim < as.numeric(start) & as.numeric(end) < up.lim,]
  some.out[, chromStart := as.integer(start)]
  some.out[, chromEnd := as.integer(end)]
  some.out[, count := as.integer(count)]
  some.peaks <- peak.dt[lo.lim < chromStart & chromEnd < up.lim,]
  some.data <- five.sparse.list[[sample.id]][lo.lim < mass.int & mass.int < up.lim,]
  some.segs <- result$segments[lo.lim < chromStart & chromEnd < up.lim,]
  fit <- PeakSegPDPAchrom(some.out, 4L)
  ggplot()+
    geom_line(aes(mass.num, count), data=some.data)+
    geom_segment(aes(
      chromStart/mult, count,
      xend=chromEnd/mult, yend=count),
                 data=some.out)+
    geom_segment(aes(
      chromStart/mult, 0,
      xend=chromEnd/mult, yend=0),
                 data=some.peaks,
                 color="deepskyblue")+
    geom_point(aes(
      chromStart/mult, 0),
               data=some.peaks,
               color="deepskyblue")+
    geom_segment(aes(
      chromStart/mult, mean,
      xend=chromEnd/mult, yend=mean),
                 data=some.segs,
                 color="green")

  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(peaks ~ .)+
    geom_line(aes(mass.num, count), data=some.data, color="grey50")+
    geom_segment(aes(
      chromStart/mult, count,
      xend=chromEnd/mult, yend=count),
                 data=some.out)+
    geom_vline(aes(
      xintercept=chromStart/mult),
               data=data.table(fit$segments)[min(chromStart) < chromStart,],
               linetype="dashed",
               color="green")+
    geom_segment(aes(
      chromStart/mult, mean,
      xend=chromEnd/mult, yend=mean),
                 data=fit$segments,
                 color="green")

  ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(peaks ~ .)+
    geom_line(aes(mass.num, log10(count)),
              size=1.5,
              data=some.data, color="grey50")+
    geom_segment(aes(
      chromStart/mult, log10(count),
      xend=chromEnd/mult, yend=log10(count)),
                 data=some.out)+
    geom_vline(aes(
      xintercept=chromStart/mult),
               data=data.table(fit$segments)[min(chromStart) < chromStart,],
               linetype="dashed",
               color="green")+
    geom_segment(aes(
      chromStart/mult, log10(mean),
      xend=chromEnd/mult, yend=log10(mean)),
                 data=fit$segments,
                 color="green")

  pen.infeasible <- t.info$models[status=="infeasible", paste(max(penalty))]
  all.out <- data.table(out.dt)
  all.out[, chromStart := as.integer(start)]
  all.out[, chromEnd := as.integer(end)]
  all.out[, count := as.integer(count)]

}

lo.lim <- 25125650
up.lim <- 25915049
some.out <- all.out[lo.lim < as.numeric(start) & as.numeric(end) < up.lim,]
some.data <- five.sparse.list[[sample.id]][lo.lim < mass.int & mass.int < up.lim,]
some.segs.list <- list()
for(pen.str in c(pen.feasible, pen.infeasible)){
  result <- problem.PeakSegFPOP(out.dir, pen.str)
  some.segs.list[[paste(result$loss$peaks)]] <- data.table(
    segments=result$loss$segments,
    result$segments[lo.lim < chromEnd & chromStart < up.lim,])
}
some.segs <- do.call(rbind, some.segs.list)
gg <- ggplot()+
  ggtitle("209 segment model is more likely with new segments on one big peak\ninstead of on several smaller peaks")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ ., labeller=label_both)+
  geom_line(aes(mass.num, count), data=some.data, color="grey50")+
  geom_segment(aes(
    chromStart/mult, count,
    xend=chromEnd/mult, yend=count),
               data=some.out)+
  geom_vline(aes(
    xintercept=chromStart/mult),
             data=some.segs[min(chromStart) < chromStart,],
             linetype="dashed",
             color="green")+
  geom_segment(aes(
    chromStart/mult, mean,
    xend=chromEnd/mult, yend=mean),
               data=some.segs,
               color="green")
pdf("figure-205-209-segments.pdf", w=9)
print(gg)
dev.off()

lo.lim <- 25700000
up.lim <- 25750000
some.out <- all.out[lo.lim < as.numeric(start) & as.numeric(end) < up.lim,]
some.data <- five.sparse.list[[sample.id]][lo.lim < mass.int & mass.int < up.lim,]
some.segs.list <- list()
for(pen.str in c(pen.feasible, pen.infeasible)){
  result <- problem.PeakSegFPOP(out.dir, pen.str)
  pen.segs <- result$segments[lo.lim < chromEnd & chromStart < up.lim,]
  pen.segs[1, chromEnd := up.lim]
  pen.segs[.N, chromStart := lo.lim]
  some.segs.list[[paste(result$loss$peaks)]] <- data.table(
    segments=result$loss$segments,
    pen.segs)
}
some.segs <- do.call(rbind, some.segs.list)
some.segs[, equality.constraint := ifelse(c(diff(mean)==0, 0), "active", "inactive"), by=segments]
some.breaks <- some.segs[min(chromStart) < chromStart,]
gg <- ggplot()+
  ggtitle(paste(
    "Model with 205 segments has 3 segments in this region (0/2 equality constraints active)",
    "Model with 209 segments has 7 segments in this region (2/6 equality constraints active)",
    sep="\n"))+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(segments ~ ., labeller=label_both)+
  geom_line(aes(mass.num, count), data=some.data, color="grey50")+
  geom_segment(aes(
    chromStart/mult, count,
    xend=chromEnd/mult, yend=count),
               data=some.out)+
  scale_linetype_manual(values=c(active="dashed", inactive="dotted"))+
  geom_vline(aes(
    linetype=equality.constraint,
    xintercept=chromStart/mult),
             data=some.breaks,
             color="green")+
  geom_segment(aes(
    chromStart/mult, mean,
    xend=chromEnd/mult, yend=mean),
               data=some.segs,
               color="green")
pdf("figure-205-209-segments-zoom.pdf", w=9)
print(gg)
dev.off()



## save(mass.spec, file="~/R/cosegData/data/mass.spec.RData", compress="xz")
## prompt(mass.spec, file="~/R/cosegData/man/mass.spec.Rd")
