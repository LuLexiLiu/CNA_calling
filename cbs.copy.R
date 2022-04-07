Args <- commandArgs()
scriptdir=Args[6]
in_binstats=Args[7]
in_binGC=Args[8]
in_Badbins=Args[9]
samplename=Args[10]
outdir=Args[11]
callpeak=Args[12]
binflag=Args[15]

defiSD=as.numeric(Args[13])
nalpha=as.numeric(Args[14])
nnperm=10000
dir.create(outdir)

library("DNAcopy")
##library("Cairo")
##library("ggplot2")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

cbs.segment01 <- function(outdir, bad.bins, varbin.gc, varbin.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width) {

	gc <- read.table(varbin.gc, header=F)

	colnames(gc) = c("bin.chrom","bin.start.chrompos",
	    "bin.start.abspos","bin.end.chrompos",
	    "bin.length","mappable.positions",
	    "gc.content")

	bad <- read.table(bad.bins, header=F)

	genome_chrs = length(table(gc[,1]))

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- as.character(genome_chrs-1)
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- as.character(genome_chrs)
	chrom.numeric <- as.numeric(chrom.numeric)

    thisRatio <- read.table(varbin.data, header=F)

	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	thisRatio$ratio <- a / mean(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)

	a <- quantile(gc$bin.length, 0.985)
	thisRatioNobad <- thisRatio[which(bad[, 1] == 0), ]
	set.seed(25)
	CNA.object <- CNA(log(thisRatioNobad$lowratio, base=2), thisRatioNobad$chrom, thisRatioNobad$chrompos, data.type="logratio", sampleid=sample.name)
    smoothed.CNA.object <- smooth.CNA(CNA.object)
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width)
    thisShort <- segment.smoothed.CNA.object[[2]]
    m <- matrix(data=0, nrow=nrow(thisRatioNobad), ncol=1)
    prevEnd <- 0
    for (i in 1:nrow(thisShort)) {
        thisStart <- prevEnd + 1
        thisEnd <- prevEnd + thisShort$num.mark[i]
        m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
        prevEnd = thisEnd
    }
    thisRatioNobad$seg.mean.LOWESS <- m[, 1]
    #thisRatioNobad[,5:10] = round(thisRatioNobad[,5:10],2)
    #write.table(thisRatioNobad, sep="\t", file=paste(outdir, "/", sample.name, ".50k", sep=""), quote=F, row.names=F)
	#write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".50k.short", sep=""), quote=F, row.names=F)
	return(list(ratioNobad=thisRatioNobad, short=thisShort))
}

#Args <- commandArgs()
infos = cbs.segment01(
    outdir=outdir,
    bad.bins=in_Badbins,
    varbin.gc=in_binGC,
    varbin.data=in_binstats,
    sample.name=samplename,
    alt.sample.name="",
    alpha=nalpha,
    nperm=nnperm,
    undo.SD=defiSD,
    min.width=5)

peaks <- function(series, span=3, ties.method = "first") {
	if((span <- as.integer(span)) %% 2 != 1) stop("'span' must be odd")
	z <- embed(series, span)
	s <- span%/%2
	v <- max.col(z, ties.method=ties.method) == 1 + s
	pad <- rep(FALSE, s)
	result <- c(pad, v, pad)
	result
}

df <- infos$ratioNobad
dfs <- infos$short

starts <- c()
ends <- c()
prevEnd <- 0
len <- nrow(dfs)
for (j in 1:len) {
	thisStart = prevEnd + 1
	thisEnd = thisStart + dfs$num.mark[j] - 1
	starts <- c(starts, thisStart)
	ends <- c(ends, thisEnd)
	prevEnd = thisEnd
}

amat <- matrix(data=0, nrow=1500000, ncol=1)
counter <- 1
for (j in 1:(len-1)) {
	for (k in (j+1):len) {
		N <- round((starts[j] - ends[j] + 1) * (starts[k] - ends[k] + 1)/1000)
		D <- abs(2^dfs$seg.mean[j] - 2^dfs$seg.mean[k])
		if (N > 0) {
			amat[(counter:(counter+N-1)), 1] <- rep.int(D, N)
			counter <- counter+N
		}
	}
}

a3 <- amat[(1:counter),1]
a3.95 <- sort(a3)[round(.95*counter)]
a3d <- density(a3[which(a3 < a3.95)], n=1000)

if(callpeak=="callpeak"){
    cn1 <- a3d$x[which(peaks(as.vector(a3d$y), span=101))][2]
} else if(callpeak=="notcallpeak"){
    cn1 <- 0.5
} else {
    cn1 <- 1/as.double(callpeak)
}

##figurefile = paste(outdir,"/Peaks.png", sep="")
##CairoPNG(figurefile, width=500, height=300)
##p2 = ggplot(data.frame(cn_pairwise = a3[which(a3 < a3.95)]),
##	aes(cn_pairwise))+
##	geom_density()
##print(p2)
##dev.off()

df$cn.ratio <- df$lowratio / cn1
df$cn.seg <- df$seg.mean.LOWESS / cn1
df$copy.number <- round(df$cn.seg)

write.table(cbind(df[,1:4], format(df[,5:11], digits=3)), sep="\t", file=paste(outdir,"/","SD_",Args[13],"_alpha_",Args[14],"_Copy.",binflag,".txt", sep=""), quote=F, row.names=F)

