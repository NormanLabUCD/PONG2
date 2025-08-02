#######################################################################
#
# Package Name: PONG Version 0.9.1
#
# Description:
#   PONG -- KIR3DL1/S1 Genotype Imputation with Attribute Bagging
#   Modified from HIBAG
#
# Author: Laura Ann Leaton & Genelle F Harrison
# License: GPL-3
# Email: GenelleFH@gmail.com
#


#######################################################################
#
# the internal functions
#
#######################################################################

#######################################################################
# get the name of HLA gene
#

.hla_gene_name_string <- function(geno.name)
{
	stopifnot(is.character(geno.name))
	i <- grep(paste("\\b",c(
			# HLA Classic I genes
			"A", "B", "C", "E", "F", "G",
			# HLA Classic II genes
			"DMA", "DMB", "DOA", "DOB",
			"DRA", "DRB\\d", "DQA\\d", "DQB\\d", "DPA\\d", "DPB\\d",
			# HLA Classic III genes
			"LTA", "TNF", "LTB", "HSPA1L", "HSPA1A", "HSPA1B",
			"C2", "BF", "C4A", "C4B"),
		"\\b", sep="", collapse="|"),
		geno.name)
	geno.name[i] <- paste("HLA -", geno.name[i])
	geno.name
}


#######################################################################
# call function in parallel
#

.DynamicClusterCall <- function(cl, fun, combine.fun, msg.fn, n,
	stop.cluster, ...)
{
	# in order to use the internal functions accessed by ':::'
	# the functions are all defined in 'parallel/R/snow.R'

	.SendData <- parse(text="parallel:::sendData(con, list(type=type,data=value,tag=tag))")
	.RecvOneData <- parse(text="parallel:::recvOneData(cl)")

	postNode <- function(con, type, value = NULL, tag = NULL)
	{
		eval(.SendData)
	}
	sendCall <- function(con, fun, args, return = TRUE, tag = NULL)
	{
		postNode(con, "EXEC",
			list(fun = fun, args = args, return = return, tag = tag))
		NULL
	}
	recvOneResult <- function(cl)
	{
		v <- eval(.RecvOneData)
		list(value = v$value$value, node = v$node, tag = v$value$tag)
	}


	#################################################################
	# check
	stopifnot(is.null(cl) | inherits(cl, "cluster"))
	stopifnot(is.function(fun))
	stopifnot(is.function(combine.fun))
	stopifnot(is.function(msg.fn))
	stopifnot(is.numeric(n))
	stopifnot(is.logical(stop.cluster))

	val <- NULL
	if (!is.null(cl))
	{
		p <- length(cl)
		if (n > 0L && p)
		{
			## **** this closure is sending to all nodes
			argfun <- function(i) c(i, list(...))

			submit <- function(node, job)
				sendCall(cl[[node]], fun, argfun(job), tag = job)

			for (i in 1:min(n, p)) submit(i, i)
			for (i in 1:n)
			{
				d <- recvOneResult(cl)
				j <- i + min(n, p)

				stopflag <- FALSE
				if (j <= n)
				{
					submit(d$node, j)
				} else {
					if (stop.cluster)
					{
						parallel::stopCluster(cl[d$node])
						cl <- cl[-d$node]
						stopflag <- TRUE
					}
				}

				dv <- d$value
				if (inherits(dv, "try-error"))
				{
					if (stop.cluster)
						parallel::stopCluster(cl)
					stop("One node produced an error: ", as.character(dv))
				}

				msg.fn(d$node, dv)
				if (stopflag)
					message(sprintf("Stop \"job %d\".", d$node))

				val <- combine.fun(val, dv)
			}
		}
	} else {
		for (i in 1:n)
		{
			dv <- fun(i, ...)
			msg.fn(i, dv)
			val <- combine.fun(val, dv)
		}
	}

	val
}



#######################################################################
#
# the functions for SNP genotypes and haplotypes
#
#######################################################################

#
#
# hlaSNPGenoClass is a class of SNP genotypes
# list:
#     genotype -- a genotype matrix, ``# of SNPs'' X ``# of samples''
#     sample.id -- sample id
#     snp.id -- snp id
#     snp.position -- snp positions in basepair
#     snp.allele -- snp alleles, ``A allele/B allele''
#     assembly -- the human genome reference, such like "hg19"
#
#
# hlaSNPHaploClass is a class of SNP haplotypes
# list:
#     haplotype -- a haplotype matrix, ``# of SNPs'' X ``2 x # of samples''
#     sample.id -- sample id
#     snp.id -- snp id
#     snp.position -- snp positions in basepair
#     snp.allele -- snp alleles, ``A allele/B allele''
#     assembly -- the human genome reference, such like "hg19"
#
#


#######################################################################
# To create a "hlaSNPGenoClass" object (SNP genotype object)
#

hlaMakeSNPGeno <- function(genotype, sample.id, snp.id, snp.position,
	A.allele, B.allele, assembly=c("auto", "hg19", "hg38", "unknown"))
{
	# check
	stopifnot(is.matrix(genotype))
	stopifnot(length(snp.id) == nrow(genotype))
	stopifnot(length(sample.id) == ncol(genotype))
	stopifnot(length(snp.id) == length(snp.position))
	stopifnot(length(snp.id) == length(A.allele))
	stopifnot(length(snp.id) == length(B.allele))
	stopifnot(is.character(A.allele))
	stopifnot(is.character(B.allele))

	assembly <- match.arg(assembly)
	if (assembly == "auto")
	{
		message("using the default genome assembly (assembly=\"hg38\")")
		assembly <- "hg38"
	}

	rv <- list(genotype = genotype, sample.id = sample.id, snp.id = snp.id,
		snp.position = snp.position,
		snp.allele = paste(A.allele, B.allele, sep="/"),
		assembly = assembly)
	class(rv) <- "hlaSNPGenoClass"

	# valid snp.id
	flag <- is.na(rv$snp.id)
	if (any(flag))
	{
		warning(sprintf(
			"There are %d SNPs with missing SNP id, and they have been removed.",
			sum(flag)))
		rv <- hlaGenoSubset(rv, snp.sel=!flag)
	}
	# valid snp.position
	flag <- is.na(rv$snp.position)
	if (any(flag))
	{
		warning(sprintf(
			"There are %d SNPs with missing SNP positions, and they have been removed.",
			sum(flag)))
		rv <- hlaGenoSubset(rv, snp.sel=!flag)
	}

	return(rv)
}


#######################################################################
# To create a "hlaSNPGenoClass" object (SNP genotype object)
#

hlaMakeSNPHaplo <- function(haplotype, sample.id, snp.id, snp.position,
	A.allele, B.allele, assembly=c("auto", "hg19", "hg38", "unknown"))
{
	# check
	stopifnot(is.matrix(haplotype))
	stopifnot(length(snp.id) == nrow(haplotype))
	stopifnot(2*length(sample.id) == ncol(haplotype))
	stopifnot(length(snp.id) == length(snp.position))
	stopifnot(length(snp.id) == length(A.allele))
	stopifnot(length(snp.id) == length(B.allele))
	stopifnot(is.character(A.allele))
	stopifnot(is.character(B.allele))

	assembly <- match.arg(assembly)
	if (assembly == "auto")
	{
		message("using the default genome assembly (assembly=\"hg38\")")
		assembly <- "hg38"
	}

	rv <- list(haplotype = haplotype, sample.id = sample.id, snp.id = snp.id,
		snp.position = snp.position,
		snp.allele = paste(A.allele, B.allele, sep="/"),
		assembly = assembly)
	class(rv) <- "hlaSNPHaploClass"

	# valid snp.id
	flag <- is.na(rv$snp.id)
	if (any(flag))
	{
		warning(sprintf(
			"There are %d SNPs with missing SNP id, and they have been removed.",
			sum(flag)))
		rv <- hlaHaploSubset(rv, snp.sel=!flag)
	}
	# valid snp.position
	flag <- is.na(rv$snp.position)
	if (any(flag))
	{
		warning(sprintf(
			"There are %d SNPs with missing SNP positions, and they have been removed.",
			sum(flag)))
		rv <- hlaHaploSubset(rv, snp.sel=!flag)
	}

	return(rv)
}


#######################################################################
# To select a subset of SNP genotypes
#

hlaGenoSubset <- function(genoobj, samp.sel=NULL, snp.sel=NULL)
{
	# check
	stopifnot(inherits(genoobj, "hlaSNPGenoClass"))
	stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.integer(samp.sel))
	if (is.logical(samp.sel))
		stopifnot(length(samp.sel) == length(genoobj$sample.id))
	stopifnot(is.null(snp.sel) | is.logical(snp.sel) | is.integer(snp.sel))
	if (is.logical(snp.sel))
		stopifnot(length(snp.sel) == length(genoobj$snp.id))
	if (is.integer(samp.sel))
	{
		stopifnot(!any(is.na(samp.sel)))
		stopifnot(length(unique(samp.sel)) == length(samp.sel))
	}
	if (is.integer(snp.sel))
	{
		stopifnot(!any(is.na(snp.sel)))
		stopifnot(length(unique(snp.sel)) == length(snp.sel))
	}

	# subset
	if (is.null(samp.sel))
		samp.sel <- rep(TRUE, length(genoobj$sample.id))
	if (is.null(snp.sel))
		snp.sel <- rep(TRUE, length(genoobj$snp.id))
	rv <- list(genotype = genoobj$genotype[snp.sel, samp.sel],
		sample.id = genoobj$sample.id[samp.sel],
		snp.id = genoobj$snp.id[snp.sel],
		snp.position = genoobj$snp.position[snp.sel],
		snp.allele = genoobj$snp.allele[snp.sel],
		assembly = genoobj$assembly
	)
	class(rv) <- "hlaSNPGenoClass"
	return(rv)
}


#######################################################################
# To select a subset of SNP haplotypes
#

hlaHaploSubset <- function(haploobj, samp.sel=NULL, snp.sel=NULL)
{
	# check
	stopifnot(inherits(haploobj, "hlaSNPHaploClass"))
	stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.integer(samp.sel))
	if (is.logical(samp.sel))
		stopifnot(length(samp.sel) == length(haploobj$sample.id))
	stopifnot(is.null(snp.sel) | is.logical(snp.sel) | is.integer(snp.sel))
	if (is.logical(snp.sel))
		stopifnot(length(snp.sel) == length(haploobj$snp.id))
	if (is.integer(samp.sel))
		stopifnot(length(unique(samp.sel)) == length(samp.sel))
	if (is.integer(snp.sel))
		stopifnot(length(unique(snp.sel)) == length(snp.sel))

	# subset
	if (is.null(samp.sel))
		samp.sel <- rep(TRUE, length(haploobj$sample.id))
	if (is.numeric(samp.sel))
	{
		v <- samp.sel
		samp.sel <- rep(FALSE, length(haploobj$sample.id))
		samp.sel[v] <- TRUE
	}
	samp.sel <- rep(samp.sel, each=2)
	if (is.null(snp.sel))
		snp.sel <- rep(TRUE, length(haploobj$snp.id))

	rv <- list(haplotype = haploobj$haplotype[snp.sel, samp.sel],
		sample.id = haploobj$sample.id[samp.sel],
		snp.id = haploobj$snp.id[snp.sel],
		snp.position = haploobj$snp.position[snp.sel],
		snp.allele = haploobj$snp.allele[snp.sel],
		assembly = haploobj$assembly
	)
	class(rv) <- "hlaSNPHaploClass"
	return(rv)
}


#######################################################################
# To select a subset of SNP genotypes
#

hlaHaplo2Geno <- function(hapobj)
{
	stopifnot(inherits(hapobj, "hlaSNPHaploClass"))
	n <- dim(hapobj$haplotype)[2]
	rv <- list(
		genotype = hapobj$haplotype[, seq(1,n,2)] + hapobj$haplotype[, seq(2,n,2)],
		sample.id = hapobj$sample.id,
		snp.id = hapobj$snp.id,
		snp.position = hapobj$snp.position,
		snp.allele = hapobj$snp.allele,
		assembly = hapobj$assembly
	)
	class(rv) <- "hlaSNPGenoClass"
	return(rv)
}


#######################################################################
# To get the overlapping SNPs between target and template with
#   corrected strand.
#

hlaGenoSwitchStrand <- function(target, template,
	match.type=c("RefSNP+Position", "RefSNP", "Position"),
	same.strand=FALSE, verbose=TRUE)
{
	# check
	stopifnot(inherits(target, "hlaSNPGenoClass") |
		inherits(target, "hlaSNPHaploClass"))
	stopifnot(inherits(template, "hlaSNPGenoClass") |
		inherits(template, "hlaSNPHaploClass") |
		inherits(template, "hlaAttrBagClass") |
		inherits(template, "hlaAttrBagObj"))
	stopifnot(is.logical(same.strand))
	stopifnot(is.logical(verbose))
	match.type <- match.arg(match.type)

	# initialize
	s1 <- hlaSNPID(template, match.type)
	s2 <- hlaSNPID(target, match.type)

	flag <- TRUE
	if (length(s1) == length(s2))
	{
		if (all(s1 == s2, na.rm=TRUE))
		{
			s <- s1
			I1 <- I2 <- seq_len(length(s1))
			flag <- FALSE
		}
	}
	if (flag)
	{
		s <- intersect(s1, s2)
		if (length(s) <= 0) stop("There is no common SNP.")
		I1 <- match(s, s1); I2 <- match(s, s2)
	}

	# compute allele frequencies
	if (inherits(template, "hlaSNPGenoClass"))
	{
		template.afreq <- rowMeans(template$genotype, na.rm=TRUE) * 0.5
	} else if (inherits(template, "hlaSNPHaploClass"))
	{
		template.afreq <- rowMeans(template$haplotype, na.rm=TRUE)
	} else {
		template.afreq <- template$snp.allele.freq
	}
	if (inherits(target, "hlaSNPGenoClass"))
	{
		target.afreq <- rowMeans(target$genotype, na.rm=TRUE) * 0.5
	} else {
		target.afreq <- rowMeans(target$haplotype, na.rm=TRUE)
	}

	# call
	gz <- .C("PONG_AlleleStrand",
		template$snp.allele, template.afreq, I1,
		target$snp.allele, target.afreq, I2,
		same.strand, length(s), out=logical(length(s)),
		out.n.ambiguity=integer(1), out.n.mismatching=integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="PONG")
	if (gz$err != 0) stop(hlaErrMsg())

	if (verbose)
	{
		# switched allele pairs
		x <- sum(gz$out)
		if (x > 0)
		{
			if (x > 1)
			{
				a <- "are"; s <- "s"
			} else {
				a <- "is"; s <- ""
			}
			cat(sprintf(
				"There %s %d variant%s in total with switched allelic strand order%s.\n",
				a, x, s, s))
		} else {
			cat("No allelic strand orders are switched.\n")
		}

		# the number of ambiguity
		if (gz$out.n.ambiguity > 0)
		{
			if (gz$out.n.ambiguity > 1)
			{
				a <- "are"; s <- "s"
			} else {
				a <- "is"; s <- ""
			}
			cat(sprintf(
				"Due to stand ambiguity (such like C/G), the allelic strand order%s of %d variant%s %s determined by comparing allele frequencies.\n",
				s, gz$out.n.ambiguity, s, a))
		}

		# the number of mismatching
		if (gz$out.n.mismatching > 0)
		{
			if (gz$out.n.mismatching > 1)
			{
				a <- "are"; s <- "s"
			} else {
				a <- "is"; s <- ""
			}
			cat(sprintf(
				"Due to mismatching alleles, the allelic strand order%s of %d variant%s %s determined by comparing allele frequencies.\n",
				s, gz$out.n.mismatching, s, a))
		}
	}

	# result
	if (inherits(target, "hlaSNPGenoClass"))
	{
		geno <- target$genotype[I2, ]
		if (is.vector(geno))
			geno <- matrix(geno, ncol=1)
		for (i in which(gz$out)) geno[i, ] <- 2 - geno[i, ]
		rv <- list(genotype = geno)
		rv$sample.id <- target$sample.id
		rv$snp.id <- target$snp.id[I2]
		rv$snp.position <- target$snp.position[I2]
		rv$snp.allele <- template$snp.allele[I1]
		rv$assembly <- template$assembly
		class(rv) <- "hlaSNPGenoClass"

	} else {
		haplo <- target$haplotype[I2, ]
		if (is.vector(haplo))
			haplo <- matrix(haplo, ncol=1)
		for (i in which(gz$out)) haplo[i, ] <- 1 - haplo[i, ]
		rv <- list(haplotype = haplo)
		rv$sample.id <- target$sample.id
		rv$snp.id <- target$snp.id[I2]
		rv$snp.position <- target$snp.position[I2]
		rv$snp.allele <- template$snp.allele[I1]
		rv$assembly <- template$assembly
		class(rv) <- "hlaSNPHaploClass"
	}

	return(rv)
}


#######################################################################
# To get the information of SNP ID and position
#

hlaSNPID <- function (obj, type = c("Position", "Pos+Allele", "RefSNP+Position", 
    "RefSNP")) 
{
    stopifnot(inherits(obj, "hlaSNPGenoClass") | inherits(obj, 
        "hlaAttrBagClass") | inherits(obj, "hlaAttrBagObj"))
    type <- match.arg(type)
    switch(type, Position = obj$snp.position, `Pos+Allele` = paste(obj$snp.position, 
        obj$snp.allele, sep = "-"), `RefSNP+Position` = paste(obj$snp.id, 
        obj$snp.position, sep = "-"), RefSNP = obj$snp.id, invisible())
}



#######################################################################
# To combine two SNP genotype dataset
#

hlaGenoCombine <- function(geno1, geno2,
	match.type=c("RefSNP+Position", "RefSNP", "Position"),
	allele.check=TRUE, same.strand=FALSE, verbose=TRUE)
{
	# check
	stopifnot(inherits(geno1, "hlaSNPGenoClass"))
	stopifnot(inherits(geno2, "hlaSNPGenoClass"))
	stopifnot(is.logical(allele.check))
	stopifnot(is.logical(same.strand))
	stopifnot(is.logical(verbose))
	match.type <- match.arg(match.type)

	if (allele.check)
	{
		tmp2 <- hlaGenoSwitchStrand(geno2, geno1, match.type,
			same.strand, verbose)
		tmp1 <- hlaGenoSubset(geno1, snp.sel=
			match(hlaSNPID(tmp2, match.type), hlaSNPID(geno1, match.type)))
	} else {
		s1 <- hlaSNPID(geno1, match.type)
		s2 <- hlaSNPID(geno2, match.type)
		set <- unique(intersect(s1, s2))
		tmp1 <- hlaGenoSubset(geno1, snp.sel=match(set, s1))
		tmp2 <- hlaGenoSubset(geno2, snp.sel=match(set, s2))
	}

	rv <- list(genotype = cbind(tmp1$genotype, tmp2$genotype),
		sample.id = c(tmp1$sample.id, tmp2$sample.id),
		snp.id = tmp1$snp.id, snp.position = tmp1$snp.position,
		snp.allele = tmp1$snp.allele,
		assembly = tmp1$assembly)
	colnames(rv$genotype) <- NULL
	rownames(rv$genotype) <- NULL

	class(rv) <- "hlaSNPGenoClass"
	return(rv)
}


#######################################################################
# Convert to PLINK PED format
#

hlaGeno2PED <- function(geno, out.fn)
{
	# check
	stopifnot(inherits(geno, "hlaSNPGenoClass"))
	stopifnot(is.character(out.fn))

	# MAP file
	rv <- data.frame(chr=rep(6, length(geno$snp.id)), rs=geno$snp.id,
		morgan=rep(0, length(geno$snp.id)), bp=geno$snp.position,
		stringsAsFactors=FALSE)
	write.table(rv, file=paste(out.fn, ".map", sep=""),
		row.names=FALSE, col.names=FALSE, quote=FALSE)

	# PED file
	n <- length(geno$sample.id)
	m <- matrix("", nrow=n, ncol=2*length(geno$snp.id))
	for (i in 1:length(geno$snp.id))
	{
		allele <- unlist(strsplit(geno$snp.allele[i], "/"))
		g <- geno$genotype[i, ] + 1
		m[, 2*(i-1)+1] <- allele[c(2, 1, 1)[g]]
		m[, 2*(i-1)+2] <- allele[c(2, 2, 1)[g]]
	}
	rv <- cbind(Family=geno$sample.id, Ind=geno$sample.id,
		Paternal=rep(0, n), Maternal=rep(0, n),
		Sex=rep(0, n), Pheno=rep(-9, n), m)
	write.table(rv, file=paste(out.fn, ".ped", sep=""),
		row.names=FALSE, col.names=FALSE, quote=FALSE)

	# return
	return(invisible(NULL))
}


#######################################################################
# Convert from PLINK BED format
#

hlaBED2Geno <- function(bed.fn, fam.fn, bim.fn, rm.invalid.allele=FALSE,
	import.chr="KIR", assembly=c("auto", "hg19", "hg38", "unknown"), ##KIRedit
	verbose=TRUE)
{
	# check
	stopifnot(is.character(bed.fn) & (length(bed.fn)==1))
	stopifnot(is.character(fam.fn) & (length(fam.fn)==1))
	stopifnot(is.character(bim.fn) & (length(bim.fn)==1))
	stopifnot(is.character(import.chr))
	stopifnot(is.logical(rm.invalid.allele) & (length(rm.invalid.allele)==1))
	stopifnot(is.logical(verbose) & (length(verbose)==1))

	assembly <- match.arg(assembly)
	if (assembly == "auto")
	{
		message("using the default genome assembly (assembly=\"hg38\")")
		warning("Please explicitly specify the argument 'assembly=\"hg19\"' or 'assembly=\"hg38\"'.")
		assembly <- "hg38"
	} else if (assembly == "unknown")
	{
		warning("Please explicitly specify the argument 'assembly=\"hg19\"' or 'assembly=\"hg38\"'.")
	}

	# detect bed.fn
	bed <- .C("PONG_BEDFlag", bed.fn, snporder=integer(1), err=integer(1),
		NAOK=TRUE, PACKAGE="PONG")
	if (bed$err != 0) stop(hlaErrMsg())
	if (verbose)
	{
		cat("Open \"", bed.fn, sep="")
		if (bed$snporder == 0)
			cat("\" in the individual-major mode.\n")
		else
			cat("\" in the SNP-major mode.\n")
	}

	# read fam.fn
	famD <- read.table(fam.fn, header=FALSE, stringsAsFactors=FALSE)
	names(famD) <- c("FamilyID", "InvID", "PatID", "MatID", "Sex", "Pheno")
	if (length(unique(famD$InvID)) == dim(famD)[1])
	{
		sample.id <- famD$InvID
	} else {
		sample.id <- paste(famD$FamilyID, famD$InvID, sep="-")
		if (length(unique(sample.id)) != dim(famD)[1])
			stop("IDs in PLINK bed are not unique!")
	}
	if (verbose)
		cat("Open \"", fam.fn, "\".\n", sep="")

	# read bim.fn
	bimD <- read.table(bim.fn, header=FALSE, stringsAsFactors=FALSE)
	names(bimD) <- c("chr", "snp.id", "map", "pos", "allele1", "allele2")

	# chromosome
	chr <- bimD$chr; chr[is.na(chr)] <- ""
	# position
	snp.pos <- bimD$pos
	snp.pos[!is.finite(snp.pos)] <- 0

	# snp.id
	snp.id <- bimD$snp.id
	if (length(snp.id) != length(unique(snp.id)))
		stop("The SNP IDs in the PLINK binary file should be unique!")

	# snp allele
	snp.allele <- paste(bimD$allele1, bimD$allele2, sep="/")
	if (verbose)
		cat("Open \"", bim.fn, "\".\n", sep="")

	# SNP selection
	if (length(import.chr) == 1)
	{
		if (import.chr == "KIR")
		{
			if (assembly %in% c("hg19", "NCBI36"))
			{
				snp.flag <- (chr==19) & (50247563<=snp.pos) & (snp.pos<=59128983) ##KIRedit
			} else if (assembly %in% c("hg38", "NCBI37"))
			{
				###snp.flag <- (chr==19) & (55247563<=snp.pos) & (snp.pos<=55361930) ##KIR region alone
				snp.flag <- (chr==19) & (46457117<=snp.pos) & (snp.pos<=58617616) ##KIR +500000
			} else {
				stop("Invalid genome assembly.")
			}
			n.snp <- as.integer(sum(snp.flag))
			if (verbose)
			{
				cat(sprintf("Import %d SNPs within the KIR gene cluster on chromosome 19.\n",
					n.snp))
			}
			import.chr <- NULL
		} else if (import.chr == "")
		{
			n.snp <- length(snp.id)
			snp.flag <- rep(TRUE, n.snp)
			if (verbose)
				cat(sprintf("Import %d SNPs.\n", n.snp))
			import.chr <- NULL
		}
	}
	if (!is.null(import.chr))
	{
		snp.flag <- (chr %in% import.chr) & (snp.pos>0)
		n.snp <- as.integer(sum(snp.flag))
		if (verbose)
		{
			cat(sprintf("Import %d SNPs from chromosome %s.\n", n.snp,
				paste(import.chr, collapse=",")))
		}
	}
	if (n.snp <= 0) stop("There is no SNP imported.")

	# call the C function
	rv <- .C("PONG_ConvBED", bed.fn, length(sample.id), length(snp.id), n.snp,
		(bed$snporder==0), snp.flag, verbose,
		geno = matrix(as.integer(0), nrow=n.snp, ncol=length(sample.id)),
		err=integer(1), NAOK=TRUE, PACKAGE="PONG")
	if (rv$err != 0) stop(hlaErrMsg())

	# result
	v <- list(genotype = rv$geno, sample.id = sample.id, snp.id = snp.id[snp.flag],
		snp.position = snp.pos[snp.flag], snp.allele = snp.allele[snp.flag],
		assembly = assembly)
	class(v) <- "hlaSNPGenoClass"

	# remove invalid snps
	if (rm.invalid.allele)
	{
		# check duplicated SNP ID
		flag <- duplicated(v$snp.id)
		if (any(flag))
		{
			if (verbose)
			{
				cat(sprintf("%d SNPs with duplicated ID have been removed.\n",
					sum(flag)))
			}
			v <- hlaGenoSubset(v, snp.sel=!flag)
		}

		# check invalid alleles
		snp.allele <- v$snp.allele
		snp.allele[is.na(snp.allele)] <- "?/?"
		flag <- sapply(strsplit(snp.allele, "/"),
			function(x)
			{
				if (length(x) == 2)
				{
					all(x %in% c("A", "G", "C", "T"))
				} else {
					FALSE
				}
			}
		)
		if (any(!flag) & verbose)
			cat(sprintf("%d SNPs with invalid alleles have been removed.\n", sum(!flag)))

		# get a subset
		v <- hlaGenoSubset(v, snp.sel=flag)
	}

	return(v)
}


#######################################################################
# Convert from SNP GDS format (SNPRelate)
#

hlaGDS2Geno <- function(gds.fn, rm.invalid.allele=FALSE,
	import.chr="xMHC", assembly=c("auto", "hg19", "hg38", "unknown"),
	verbose=TRUE)
{
	# library
	eval(parse(text='
		if (!require(SNPRelate))
			stop("The SNPRelate package should be installed.")
	'))

	# check
	stopifnot(is.character(gds.fn) & is.vector(gds.fn))
	stopifnot(length(gds.fn) == 1)

	stopifnot(is.logical(rm.invalid.allele) & is.vector(rm.invalid.allele))
	stopifnot(length(rm.invalid.allele) == 1)

	stopifnot(is.character(import.chr))

	stopifnot(is.logical(verbose) & is.vector(verbose))
	stopifnot(length(verbose) == 1)

	assembly <- match.arg(assembly)
	if (assembly == "auto")
	{
		message("using the default genome assembly (assembly=\"hg38\")")
		warning("Please explicitly specify the argument 'assembly=\"hg19\"' or 'assembly=\"hg38\"'.")
		assembly <- "hg38"
	} else if (assembly == "unknown")
	{
		warning("Please explicitly specify the argument 'assembly=\"hg19\"' or 'assembly=\"hg38\"'.")
	}


	####  open the GDS SNP file  ####

	chr <- NULL
	snp.pos <- NULL
	snp.id <- NULL

	eval(parse(text='gfile <- snpgdsOpen(gds.fn)'))
	on.exit(eval(parse(text='snpgdsClose(gfile)')))

	# snp.id
	eval(parse(text='
		snp.id <- read.gdsn(index.gdsn(gfile, "snp.id"))
		if (!is.null(index.gdsn(gfile, "snp.rs.id", silent=TRUE)))
		{
			snp.rsid <- read.gdsn(index.gdsn(gfile, "snp.rs.id"))
		} else
			snp.rsid <- snp.id
	'))

	# chromosome
	eval(parse(text='
		chr <- read.gdsn(index.gdsn(gfile, "snp.chromosome"))
	'))

	# position
	eval(parse(text='
		snp.pos <- read.gdsn(index.gdsn(gfile, "snp.position"))
		snp.pos[!is.finite(snp.pos)] <- 0
	'))

	# SNP selection
	if (length(import.chr) == 1)
	{
		if (import.chr == "xMHC")
		{
			if (assembly %in% c("hg19", "NCBI36"))
			{
				snp.flag <- (chr==6) & (25759242<=snp.pos) & (snp.pos<=33534827)
			} else if (assembly %in% c("hg38", "NCBI37"))
			{
				snp.flag <- (chr==6) & (25651242<=snp.pos) & (snp.pos<=33544122)
			} else {
				stop("Invalid genome assembly.")
			}
			n.snp <- as.integer(sum(snp.flag))
			if (verbose)
			{
				cat(sprintf("Import %d SNPs within the xMHC region on chromosome 6.\n",
					n.snp))
			}
			import.chr <- NULL
		} else if (import.chr == "")
		{
			n.snp <- length(snp.id)
			snp.flag <- rep(TRUE, n.snp)
			if (verbose)
				cat(sprintf("Import %d SNPs.\n", n.snp))
			import.chr <- NULL
		}
	}
	if (!is.null(import.chr))
	{
		snp.flag <- (chr %in% import.chr) & (snp.pos>0)
		n.snp <- as.integer(sum(snp.flag))
		if (verbose)
		{
			cat(sprintf("Import %d SNPs from chromosome %s.\n", n.snp,
				paste(import.chr, collapse=",")))
		}
	}
	if (n.snp <= 0) stop("There is no SNP imported.")

	# result
	eval(parse(text='
		v <- list(genotype = snpgdsGetGeno(gfile, snp.id=snp.id[snp.flag],
				snpfirstdim=TRUE, verbose=FALSE),
			sample.id = read.gdsn(index.gdsn(gfile, "sample.id")),
			snp.id = snp.rsid[snp.flag],
			snp.position = snp.pos[snp.flag],
			snp.allele = read.gdsn(index.gdsn(gfile, "snp.allele"))[snp.flag],
			assembly = assembly)
	'))
	class(v) <- "hlaSNPGenoClass"

	# remove invalid snps
	if (rm.invalid.allele)
	{
		# check duplicated SNP ID
		flag <- duplicated(v$snp.id)
		if (any(flag))
		{
			if (verbose)
			{
				cat(sprintf("%d SNPs with duplicated ID have been removed.\n",
					sum(flag)))
			}
			v <- hlaGenoSubset(v, snp.sel=!flag)
		}

		# check invalid alleles
		snp.allele <- v$snp.allele
		snp.allele[is.na(snp.allele)] <- "?/?"
		flag <- sapply(strsplit(snp.allele, "/"),
			function(x)
			{
				if (length(x) == 2)
				{
					all(x %in% c("A", "G", "C", "T"))
				} else {
					FALSE
				}
			}
		)
		if (any(!flag) & verbose)
		{
			cat(sprintf("%d SNPs with invalid alleles have been removed.\n",
				sum(!flag)))
		}
		v <- hlaGenoSubset(v, snp.sel=flag)
	}

	return(v)
}



#######################################################################
# Summarize a "hlaSNPGenoClass" object
#
#' @export
summary.hlaSNPGenoClass <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(inherits(object, "hlaSNPGenoClass"))
	geno <- object

	fn <- function(x)
	{
		sprintf("min: %g, max: %g, mean: %g, median: %g, sd: %g",
			min(x, na.rm=TRUE), max(x, na.rm=TRUE),
			mean(x, na.rm=TRUE), median(x, na.rm=TRUE), sd(x, na.rm=TRUE))
	}

	rv <- list(mr.snp = hlaGenoMRate(geno), mr.samp = hlaGenoMRate_Samp(geno),
		maf = hlaGenoMFreq(geno),
		allele = table(geno$snp.allele))

	if (show)
	{
		cat("SNP genotypes: \n")
		cat(sprintf("\t%d samples X %d SNPs\n",
			length(geno$sample.id), length(geno$snp.id)))
		cat(sprintf("\tSNPs range from %dbp to %dbp",
			min(geno$snp.position, na.rm=TRUE),
			max(geno$snp.position, na.rm=TRUE)))
		if (!is.null(geno$assembly))
			cat(" on ", geno$assembly, "\n", sep="")
		else
			cat("\n")

		# missing rate for SNP
		cat(sprintf("Missing rate per SNP:\n\t%s\n", fn(rv$mr.snp)))
		# missing rate for sample
		cat(sprintf("Missing rate per sample:\n\t%s\n", fn(rv$mr.samp)))

		# minor allele frequency
		cat(sprintf("Minor allele frequency:\n\t%s\n", fn(rv$maf)))

		# allele information
		cat("Allelic information:")
		print(rv$allele)
	}

	# return
	return(invisible(rv))
}

#######################################################################
# Summarize a "hlaSNPHaploClass" object
#
#' @export
summary.hlaSNPHaploClass <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(inherits(object, "hlaSNPHaploClass"))
	haplo <- object

	fn <- function(x)
	{
		sprintf("min: %g, max: %g, mean: %g, median: %g, sd: %g",
			min(x, na.rm=TRUE), max(x, na.rm=TRUE),
			mean(x, na.rm=TRUE), median(x, na.rm=TRUE), sd(x, na.rm=TRUE))
	}

	rv <- list(mr.snp = hlaGenoMRate(haplo), mr.samp = hlaGenoMRate_Samp(haplo),
		maf = hlaGenoMFreq(haplo),
		allele = table(haplo$snp.allele))

	if (show)
	{
		cat("SNP Haplotypes: \n")
		cat(sprintf("\t%d samples X %d SNPs\n",
			length(haplo$sample.id), length(haplo$snp.id)))
		cat(sprintf("\tSNPs range from %dbp to %dbp",
			min(haplo$snp.position, na.rm=TRUE), max(haplo$snp.position, na.rm=TRUE)))
		if (!is.null(haplo$assembly))
			cat(" on ", haplo$assembly, "\n", sep="")
		else
			cat("\n")

		# missing rate for SNP
		cat(sprintf("Missing rate per SNP:\n\t%s\n", fn(rv$mr.snp)))
		# missing rate for sample
		cat(sprintf("Missing rate per sample:\n\t%s\n", fn(rv$mr.samp)))

		# minor allele frequency
		cat(sprintf("Minor allele frequency:\n\t%s\n", fn(rv$maf)))

		# allele information
		cat("Allelic information:")
		print(rv$allele)
	}

	# return
	return(invisible(rv))
}




#######################################################################
#
# the function list for genotypes and haplotypes
#
#######################################################################

#######################################################################
# To the allele frequencies from genotypes or haplotypes
#

hlaGenoAFreq <- function(obj)
{
	# check
	stopifnot(inherits(obj, "hlaSNPGenoClass") |
		inherits(obj, "hlaSNPHaploClass"))
    if (inherits(obj, "hlaSNPGenoClass"))
    {
        rowMeans(obj$genotype, na.rm=TRUE) * 0.5
    } else {
        rowMeans(obj$haplotype, na.rm=TRUE)
    }
}


#######################################################################
# To the minor allele frequencies from genotypes or haplotypes
#

hlaGenoMFreq <- function(obj)
{
	# check
	stopifnot(inherits(obj, "hlaSNPGenoClass") |
		inherits(obj, "hlaSNPHaploClass"))
    if (inherits(obj, "hlaSNPGenoClass"))
    {
        F <- rowMeans(obj$genotype, na.rm=TRUE) * 0.5
    } else {
        F <- rowMeans(obj$haplotype, na.rm=TRUE)
    }
    pmin(F, 1-F)
}


#######################################################################
# To the missing rates from genotypes or haplotypes per SNP
#

hlaGenoMRate <- function(obj)
{
	# check
	stopifnot(inherits(obj, "hlaSNPGenoClass") |
		inherits(obj, "hlaSNPHaploClass"))
    if (inherits(obj, "hlaSNPGenoClass"))
    {
		rowMeans(is.na(obj$genotype))
	} else {
		rowMeans(is.na(obj$haplotype))
	}
}


#######################################################################
# To the missing rates from genotypes or haplotypes per sample
#

hlaGenoMRate_Samp <- function(obj)
{
	# check
	stopifnot(inherits(obj, "hlaSNPGenoClass") |
		inherits(obj, "hlaSNPHaploClass"))
    if (inherits(obj, "hlaSNPGenoClass"))
    {
        colMeans(is.na(obj$genotype))
    } else {
        F <- colMeans(is.na(obj$haplotype))
        F[seq(1, length(F), 2)] + F[seq(2, length(F), 2)]
    }
}


#######################################################################
#
# the function list for KIR3DL1/S1 types
#
#######################################################################

#######################################################################
# To get the starting and ending positions in basepair for HLA loci
#

hlaLociInfo <- function(assembly =
	c("auto", "auto-silent", "hg19", "hg38", "unknown"))
	{
		# check
		assembly <- match.arg(assembly)
		if (assembly == "auto")
		{
			message("using the default genome assembly (assembly=\"hg38\")")
			assembly <- "hg38"
		} else if (assembly == "auto-silent")
		{
			assembly <- "hg38"
		}


		if (assembly == "hg19")
		{
			#print(sprintf("GeneInfo_%s.txt", assembly))
			fn <- system.file("extdata", sprintf("GeneInfo_%s.txt", assembly), package = "PONG")
			if (file.exists(fn)){
				v <- read.table(fn, header = TRUE, stringsAsFactors = FALSE)
				rownames(v) <- v$name
				v[, -1L]
				# return
			
				return(v)
			}
			
		}else if (assembly == "hg38")
		{
			fn <- system.file("doc", sprintf("GeneInfo_%s.txt", assembly), package = "PONG")
			if (file.exists(fn)){
				v <- read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
				rownames(v) <- v$name
				v[, -1L]
				# return
				return(v)
			}
		}else {
		stop("Unknown human genome reference in 'assembly'!")
		}

	} 


#######################################################################
# Limit the resolution of KIR3DL1/S1 alleles
#

hlaAlleleDigit <- function(obj, max.resolution="6-digit", rm.suffix=FALSE)
{
	# check
	stopifnot(inherits(obj, "hlaAlleleClass") | is.character(obj))
	stopifnot(is.logical(rm.suffix))
	if (is.character(obj))
		stopifnot(is.vector(obj))
	stopifnot(max.resolution %in% c("2-digit", "3-digit","4-digit", "5-digit","6-digit", "8-digit",
		"allele", "protein", "2", "4", "6", "8", "full", ""))

	if (!(max.resolution %in% c("full", "")))
	{
		if (is.character(obj))
		{
			len <- c(1, 2, 3, 4, 5, 6, 1, 2, 1, 2, 3, 4)
			names(len) <- c("2-digit", "3-digit","4-digit", "5-digit","6-digit", "8-digit",
				"allele", "protein", "2", "4", "6", "8")
			maxlen <- len[[as.character(max.resolution)]]

			obj <- sapply(strsplit(obj, ":"), FUN =
					function(s, idx) {
						if (any(is.na(s)))
						{
							NA
						} else {
							if (length(idx) < length(s)) s <- s[idx]
							if (length(idx) == length(s))
							{
								if (rm.suffix & (nchar(s[length(s)])>0))
								{
									z <- unlist(strsplit(s[length(s)], ""))
									for (i in 1:length(z))
									{
										if (!(z[i] %in% as.character(0:9)))
										{
											if (i > 1)
												z <- z[1:(i-1)]
											break
										}
									}
									s[length(s)] <- paste(z, collapse="")
								}
							}
							paste(s, collapse=":")
						}
					},
				idx = 1:maxlen)
		} else {
			rv <- list(locus = obj$locus,
				pos.start = obj$pos.start, pos.end = obj$pos.end,
				value = data.frame(sample.id = obj$value$sample.id,
					allele1 = hlaAlleleDigit(obj$value$allele1, max.resolution),
					allele2 = hlaAlleleDigit(obj$value$allele2, max.resolution),
					stringsAsFactors=FALSE),
				assembly = obj$assembly
			)
			if ("prob" %in% names(obj$value))
				rv$value$prob <- obj$value$prob
			class(rv) <- "hlaAlleleClass"
			obj <- rv
		}
	}

	return(obj)
}


#######################################################################
# Get unique KIR3DL1/S1 alleles
#

hlaUniqueAllele <- function(hla)
{
	# check
	stopifnot(is.character(hla) | inherits(hla, "hlaAlleleClass"))

	if (is.character(hla))
	{
		hla <- hla[!is.na(hla)]
		hla <- unique(hla)
		rv <- .C("PONG_SortAlleleStr", length(hla), hla,
			out = character(length(hla)),
			err = integer(1), NAOK = TRUE, PACKAGE = "PONG")
		if (rv$err != 0) stop(hlaErrMsg())
		rv$out
	} else {
		hlaUniqueAllele(as.character(c(hla$value$allele1, hla$value$allele2)))
	}
}


#######################################################################
# To make a class of KIR3DL1/S1 alleles
#

hlaAllele <- function(sample.id, H1, H2, max.resolution="", locus="any",
	assembly=c("auto", "auto-silent", "hg19", "hg38", "unknown"),
	locus.pos.start=NA, locus.pos.end=NA, prob=NULL, na.rm=TRUE)
{
	# check
	stopifnot(is.vector(sample.id))
	stopifnot(is.vector(H1) & is.character(H1))
	stopifnot(is.vector(H2) & is.character(H2))
	stopifnot(length(sample.id) == length(H1))
	stopifnot(length(sample.id) == length(H2))
	stopifnot(max.resolution %in% c("2-digit", "3-digit","4-digit", "5-digit","6-digit",
		"8-digit", "allele", "protein", "2", "4", "6", "8", "full", ""))

	HLAinfo <- hlaLociInfo(assembly)
	if (!is.null(prob))
		stopifnot(length(sample.id) == length(prob))

	# build
	H1[H1 == ""] <- NA
	H1 <- hlaAlleleDigit(H1, max.resolution)
	H2 <- hlaAlleleDigit(H2, max.resolution)

	if (locus %in% names(HLAinfo$pos.HLA.start))
	{
		if (!is.finite(locus.pos.start))
			locus.pos.start <- HLAinfo$pos.HLA.start[[locus]]
		if (!is.finite(locus.pos.end))
			locus.pos.end <- HLAinfo$pos.HLA.end[[locus]]
	}

	# remove missing values
	if (na.rm)
		flag <- (!is.na(H1)) & (!is.na(H2))
	else
		flag <- rep(TRUE, length(sample.id))

	# result
	rv <- list(locus = locus,
		pos.start = locus.pos.start, pos.end = locus.pos.end,
		value = data.frame(sample.id = sample.id[flag],
			allele1 = H1[flag], allele2 = H2[flag],
			stringsAsFactors=FALSE),
		assembly = HLAinfo$assembly
	)
	if (!is.null(prob))
		rv$value$prob <- prob[flag]
	class(rv) <- "hlaAlleleClass"
	return(rv)
}


#######################################################################
# To make a class of KIR3DL1/S1 alleles
#

hlaAlleleSubset <- function(hla, samp.sel=NULL)
{
	# check
	stopifnot(inherits(hla, "hlaAlleleClass"))
	stopifnot(is.null(samp.sel) | is.logical(samp.sel) | is.integer(samp.sel))
	if (is.logical(samp.sel))
		stopifnot(length(samp.sel) == dim(hla$value)[1])
	if (is.integer(samp.sel))
		stopifnot(length(unique(samp.sel)) == length(samp.sel))

	# result
	if (is.null(samp.sel))
		samp.sel <- rep(TRUE, dim(hla$value)[1])
	rv <- list(locus = hla$locus,
		pos.start = hla$pos.start, pos.end = hla$pos.end,
		value = hla$value[samp.sel, ],
		assembly = hla$assembly
	)

	if (!is.null(hla$postprob))
	{
		rv$postprob <- hla$postprob[, samp.sel]
		if (is.vector(rv$postprob))
		{
			rv$postprob <- matrix(rv$postprob, ncol=1)
		}
	}

	class(rv) <- "hlaAlleleClass"
	return(rv)
}


#######################################################################
# To combine two classes of KIR3DL1/S1 alleles
#

hlaCombineAllele <- function(H1, H2)
{
	# check
	stopifnot(inherits(H1, "hlaAlleleClass"))
	stopifnot(inherits(H2, "hlaAlleleClass"))
	stopifnot(length(intersect(H1$sample.id, H2$sample.id)) == 0)
	stopifnot(H1$locus == H2$locus)
	stopifnot(H1$pos.start == H2$pos.start)
	stopifnot(H1$pos.end == H2$pos.end)

	id <- c("sample.id", "allele1", "allele2")

	# result
	rv <- list(locus = H1$locus,
		pos.start = H1$pos.start, pos.end = H1$pos.end,
		value = rbind(H1$value[, id], H2$value[, id]),
		assembly = H1$assembly
	)
	rownames(rv$value) <- NULL
	if (!is.null(H1$value$prob) & !is.null(H2$value$prob))
	{
		rv$value$prob <- c(H1$value$prob, H2$value$prob)
	}

	if (!is.null(H1$postprob) & !is.null(H2$postprob))
	{
		rv$postprob <- cbind(H1$postprob, H2$postprob)
	}

	class(rv) <- "hlaAlleleClass"
	return(rv)
}


#######################################################################
# To compare KIR3DL1/S1 alleles
#

hlaCompareAllele <- function(TrueHLA, PredHLA, allele.limit=NULL,
	call.threshold=NaN, max.resolution="", output.individual=FALSE,
	verbose=TRUE)
{
	# check
	stopifnot(inherits(TrueHLA, "hlaAlleleClass"))
	stopifnot(inherits(PredHLA, "hlaAlleleClass"))
	stopifnot(is.null(allele.limit) |
		(is.vector(allele.limit) & is.character(allele.limit)) |
		inherits(allele.limit, "hlaAttrBagClass") |
		inherits(allele.limit, "hlaAttrBagObj"))
	stopifnot(max.resolution %in% c("2-digit", "3-digit","4-digit", "5-digit","6-digit",
		"8-digit", "allele", "protein", "2", "4", "6", "8", "full", ""))
	stopifnot(is.logical(output.individual))
	stopifnot(is.logical(verbose))

	# get the common samples
	samp <- intersect(TrueHLA$value$sample.id, PredHLA$value$sample.id)
	if ((length(samp) != length(TrueHLA$value$sample.id)) |
		(length(samp) != length(PredHLA$value$sample.id)))
	{
		if (verbose)
		{
			message("Calling 'hlaCompareAllele': there are ", length(samp),
				" individuals in common.\n")
		}
	}
	# True KIR3DL1/S1
	flag <- match(samp, TrueHLA$value$sample.id)
	if (length(samp) != length(TrueHLA$value$sample.id))
	{
		TrueHLA <- hlaAlleleSubset(TrueHLA, flag)
	} else {
		if (!all(flag == 1:length(TrueHLA$value$sample.id)))
			TrueHLA <- hlaAlleleSubset(TrueHLA, flag)
	}
	# Predicted KIR3DL1/S1
	flag <- match(samp, PredHLA$value$sample.id)
	if (length(samp) != length(PredHLA$value$sample.id))
	{
		PredHLA <- hlaAlleleSubset(PredHLA, flag)
	} else {
		if (!all(flag == 1:length(PredHLA$value$sample.id)))
			PredHLA <- hlaAlleleSubset(PredHLA, flag)
	}

	# init
	flag <- !is.na(TrueHLA$value$allele1) & !is.na(TrueHLA$value$allele2) &
		!is.na(PredHLA$value$allele1) & !is.na(PredHLA$value$allele2)
	ts1 <- TrueHLA$value$allele1[flag]; ts2 <- TrueHLA$value$allele2[flag]
	ps1 <- PredHLA$value$allele1[flag]; ps2 <- PredHLA$value$allele2[flag]
	samp.id <- TrueHLA$value$sample.id[flag]

	# call threshold
	if (is.finite(call.threshold))
	{
		prob <- PredHLA$value$prob
		if (!is.null(prob)) prob <- prob[flag]
	} else {
		prob <- NULL
	}

	# allele limitation
	if (!is.null(allele.limit))
	{
		if (inherits(allele.limit, "hlaAttrBagClass") |
			inherits(allele.limit, "hlaAttrBagObj"))
		{
			allele <- hlaUniqueAllele(allele.limit$hla.allele)
			TrainFreq <- allele.limit$hla.freq
			TrainNum <- allele.limit$n.samp
		} else {
			allele <- hlaUniqueAllele(as.character(allele.limit))
			TrainFreq <- NULL
			TrainNum <- NaN
		}
	} else {
		allele <- hlaUniqueAllele(c(ts1, ts2))
		TrainFreq <- NULL
		TrainNum <- NaN
	}

	# max resolution
	if (!(max.resolution %in% c("full", "")))
	{
		ts1 <- hlaAlleleDigit(ts1, max.resolution)
		ts2 <- hlaAlleleDigit(ts2, max.resolution)
		ps1 <- hlaAlleleDigit(ps1, max.resolution)
		ps2 <- hlaAlleleDigit(ps2, max.resolution)
		tmp <- hlaAlleleDigit(allele, max.resolution)
		allele <- hlaUniqueAllele(tmp)
		if ((length(tmp) != length(allele)) & !is.null(TrainFreq))
		{
			x <- rep(0, length(allele))
			for (i in 1:length(allele))
				x[i] <- sum(TrainFreq[tmp == allele[i]])
			TrainFreq <- x
		}
	}

	# allele filter
	flag <- (ts1 %in% allele) & (ts2 %in% allele)
	ts1 <- ts1[flag]; ts2 <- ts2[flag]
	ps1 <- ps1[flag]; ps2 <- ps2[flag]
	samp.id <- samp.id[flag]
	if (!is.null(prob)) prob <- prob[flag]

	# init ...
	cnt.ind <- 0; cnt.haplo <- 0; cnt.call <- 0
	n <- length(ts1)
	m <- length(allele)

	TrueNum <- rep(0, m); names(TrueNum) <- allele
	TrueNumAll <- rep(0, m); names(TrueNumAll) <- allele
	PredNum <- rep(0, m+1); names(PredNum) <- c(allele, "...")
	confusion <- matrix(0.0, nrow = m+1, ncol = m,
		dimnames = list(Predict=names(PredNum), True=names(TrueNum)))
	WrongTab <- NULL

	# for PredNum
	fn <- function(x, LT)
		{ if (x %in% LT) { return(x) } else { return("...") } }

	acc.array <- rep(NaN, n)
	ind.truehla <- character(n)
	ind.predhla <- character(n)

	if (n > 0)
	{
		for (i in 1:n)
		{
			# increase
			TrueNumAll[[ ts1[i] ]] <- TrueNumAll[[ ts1[i] ]] + 1
			TrueNumAll[[ ts2[i] ]] <- TrueNumAll[[ ts2[i] ]] + 1

			# probability cut-off
			if (is.null(prob))
				flag <- TRUE
			else
				flag <- (prob[i] >= call.threshold)
			if (flag)
			{
				# update TrueNum and PredNum
				TrueNum[[ ts1[i] ]] <- TrueNum[[ ts1[i] ]] + 1
				TrueNum[[ ts2[i] ]] <- TrueNum[[ ts2[i] ]] + 1
				PredNum[[ fn(ps1[i], allele) ]] <- PredNum[[ fn(ps1[i], allele) ]] + 1
				PredNum[[ fn(ps2[i], allele) ]] <- PredNum[[ fn(ps2[i], allele) ]] + 1

				# correct count of individuals
				if ( ((ts1[i]==ps1[i]) & (ts2[i]==ps2[i])) |
					((ts2[i]==ps1[i]) & (ts1[i]==ps2[i])) )
				{
					cnt.ind <- cnt.ind + 1
				}

				# correct count of haplotypes
				s <- c(ts1[i], ts2[i]); p <- c(ps1[i], ps2[i])
				ind.truehla[i] <- paste(s[order(s)], collapse="/")
				ind.predhla[i] <- paste(p[order(p)], collapse="/")
				
				hnum <- 0
				if ((s[1]==p[1]) | (s[1]==p[2]))
				{
					if (s[1]==p[1]) { p[1] <- "" } else { p[2] <- "" }
					confusion[s[1], s[1]] <- confusion[s[1], s[1]] + 1
					cnt.haplo <- cnt.haplo + 1
					hnum <- hnum + 1
				}
				if ((s[2]==p[1]) | (s[2]==p[2]))
				{
					confusion[s[2], s[2]] <- confusion[s[2], s[2]] + 1
					cnt.haplo <- cnt.haplo + 1
					hnum <- hnum + 1
				}
				acc.array[i] <- 0.5*hnum

				# for confusion matrix
				s <- c(ts1[i], ts2[i]); p <- c(ps1[i], ps2[i])
				if (hnum == 1)
				{
					if ((s[1]==p[1]) | (s[1]==p[2]))
					{
						if (s[1]==p[1])
						{
							confusion[fn(p[2], allele), s[2]] <-
								confusion[fn(p[2], allele), s[2]] + 1
						} else {
							confusion[fn(p[1], allele), s[2]] <-
								confusion[fn(p[1], allele), s[2]] + 1
						}
					} else {
						if (s[2]==p[1])
						{
							confusion[fn(p[2], allele), s[1]] <-
								confusion[fn(p[2], allele), s[1]] + 1
						} else {
							confusion[fn(p[1], allele), s[1]] <-
								confusion[fn(p[1], allele), s[1]] + 1
						}
					}
				} else if (hnum == 0)
				{
					WrongTab <- cbind(WrongTab,
						c(s, fn(p[1], allele), fn(p[2], allele)))
				}

				# the number of calling
				cnt.call <- cnt.call + 1
			}
		}
	}

	# overall
	overall <- data.frame(total.num.ind = n,
		crt.num.ind = cnt.ind, crt.num.haplo = cnt.haplo,
		acc.ind = cnt.ind/cnt.call, acc.haplo = 0.5*cnt.haplo/cnt.call,
		call.threshold = call.threshold)
	if (is.finite(call.threshold))
	{
		overall$n.call <- cnt.call
		overall$call.rate <- cnt.call / n
	} else {
		overall$n.call <- n
		overall$call.rate <- 1.0
		overall$call.threshold <- 0
	}

	# confusion matrix
	if (is.null(WrongTab))
	{
		nw <- as.integer(0)
	} else {
		nw <- ncol(WrongTab)
	}
	rv <- .C("PONG_Confusion", as.integer(m), confusion,
		nw, match(WrongTab, names(PredNum)) - as.integer(1),
		out = matrix(0.0, nrow=m+1, ncol=m, dimnames=
			list(Predict=names(PredNum), True=names(TrueNum))),
		tmp = double((m+1)*m),
		err = integer(1), NAOK = TRUE, PACKAGE = "PONG")
	if (rv$err != 0) stop(hlaErrMsg())
	confusion <- round(rv$out, 2)

	# detail -- sensitivity and specificity
	detail <- data.frame(allele = allele, stringsAsFactors=FALSE)
	if (!is.null(TrainFreq))
	{
		detail$train.num <- 2 * TrainFreq * TrainNum
		detail$train.freq <- TrainFreq
	}
	detail$valid.num <- TrueNumAll
	detail$valid.freq <- TrueNumAll / sum(TrueNumAll)
	detail$call.rate <- TrueNum / TrueNumAll

	sens <- diag(confusion) / TrueNum
	spec <- 1 - (PredNum[1:m] - diag(confusion)) / (2*cnt.call - TrueNum)
	detail$accuracy <- (sens*TrueNum + spec*(2*cnt.call - TrueNum)) / (2*cnt.call)
	detail$sensitivity <- sens
	detail$specificity <- spec
	detail$ppv <- diag(confusion) / rowSums(confusion)[1:m]
	detail$npv <- 1 - (TrueNum - diag(confusion)) / (2*n - rowSums(confusion)[1:m])

	detail$call.rate[!is.finite(detail$call.rate)] <- 0
	detail[detail$call.rate<=0,
		c("sensitivity", "specificity", "ppv", "npv", "accuracy")] <- NaN

	# get miscall
	rv <- confusion; diag(rv) <- 0
	m.max <- apply(rv, 2, max); m.idx <- apply(rv, 2, which.max)
	s <- names(PredNum)[m.idx]; s[m.max<=0] <- NA
	p <- m.max / apply(rv, 2, sum)
	detail <- cbind(detail, miscall=s, miscall.prop=p, stringsAsFactors=FALSE)
	rownames(detail) <- NULL

	# output
	rv <- list(overall=overall, confusion=confusion, detail=detail)
	if (output.individual)
	{
		rv$individual <- data.frame(sample.id=samp.id,
			true.hla=ind.truehla, pred.hla=ind.predhla,
			accuracy=acc.array, stringsAsFactors=FALSE)
	}
	return(rv)
}


#######################################################################
# Return a sample list satisfying the filter conditions
#

hlaSampleAllele <- function(TrueHLA, allele.limit = NULL, max.resolution="")
{
	# check
	stopifnot(inherits(TrueHLA, "hlaAlleleClass"))
	stopifnot(is.null(allele.limit) | is.vector(allele.limit) |
		inherits(allele.limit, "hlaAttrBagClass") |
		inherits(allele.limit, "hlaAttrBagObj"))
	stopifnot(max.resolution %in% c("2-digit", "3-digit","4-digit", "5-digit","6-digit",
		"8-digit", "allele", "protein", "2", "4", "6", "8", "full", ""))

	# init
	flag <- !is.na(TrueHLA$value$allele1) & !is.na(TrueHLA$value$allele2)
	ts1 <- TrueHLA$value$allele1[flag]
	ts2 <- TrueHLA$value$allele2[flag]

	# max resolution
	if (!(max.resolution %in% c("full", "")))
	{
		ts1 <- hlaAlleleDigit(ts1, max.resolution)
		ts2 <- hlaAlleleDigit(ts2, max.resolution)
	}

	# allele limitation
	if (!is.null(allele.limit))
	{
		if (inherits(allele.limit, "hlaAttrBagClass") |
			inherits(allele.limit, "hlaAttrBagObj"))
		{
			allele <- levels(factor(allele.limit$hla.allele))
		} else {
			allele <- levels(factor(allele.limit))
		}
		if (!(max.resolution %in% c("full", "")))
		{
			allele <- hlaAlleleDigit(allele, max.resolution)
		}
		flag[flag] <- (ts1 %in% allele) & (ts2 %in% allele)
	}

	# return
	return(TrueHLA$value$sample.id[flag])
}


#######################################################################
# Divide the list of KIR3DL1/S1 types to the training and validation sets
#

hlaSplitAllele <- function(HLA, train.prop=0.5)
{
	# check
	stopifnot(inherits(HLA, "hlaAlleleClass"))

	train.set <- NULL
	H <- HLA
	while (dim(H$value)[1] > 0)
	{
		v <- summary(H, show=FALSE)
		if (dim(v)[1] > 1)
		{
			v <- v[order(v[, "count"]), ]
		}

		allele <- rownames(v)[1]
		samp.id <- H$value$sample.id[ H$value$allele1==allele |
			H$value$allele2==allele ]
		samp.id <- as.character(samp.id)

		n.train <- ceiling(length(samp.id) * train.prop)
		train.sampid <- sample(samp.id, n.train)
		train.set <- c(train.set, train.sampid)

		H <- hlaAlleleSubset(H, samp.sel=
			match(setdiff(H$value$sample.id, samp.id), H$value$sample.id))
	}

	train.set <- train.set[order(train.set)]
	list(
		training = hlaAlleleSubset(HLA,
			samp.sel = match(train.set, HLA$value$sample.id)),
		validation = hlaAlleleSubset(HLA,
			samp.sel = match(setdiff(HLA$value$sample.id, train.set),
				HLA$value$sample.id)
			)
	)
}


#######################################################################
# To select SNPs in the flanking region of a specified LRC locus
#

hlaGenoSubsetFlank <- function (genoobj, locus = "any", flank.bp = 500000L, assembly=c("auto", "hg19", "hg38", "unknown"), 
    pos.mid = NA_integer_) 
{
    stopifnot(inherits(genoobj, "hlaSNPGenoClass"))
    stopifnot(is.character(locus), length(locus) == 1L, !is.na(locus))
    stopifnot(is.numeric(flank.bp), length(flank.bp) == 1L, is.finite(flank.bp))
    stopifnot(is.character(assembly), length(assembly) == 1L, 
        !is.na(assembly))
    stopifnot(is.numeric(pos.mid), length(pos.mid) == 1L)
    if (locus != "any") {
        assembly <- match.arg(assembly)
        HLAInfo <- hlaLociInfo(assembly)
        ID <- rownames(HLAInfo)
        if (!(locus %in% ID)) 
            stop(paste("'locus' should be one of", paste(ID, 
                collapse = ", ")))
        pos.start <- HLAInfo[locus, "start"] - flank.bp
        pos.end <- HLAInfo[locus, "end"] + flank.bp
        if (!is.na(pos.mid)) 
            warning("'pos.mid' is ignored when 'locus' is specified.", 
                immediate. = TRUE)
        if (!is.null(HLAInfo$suggest.pos)) {
            i <- HLAInfo[locus, "suggest.pos"]
            if (!is.na(i) && i < 0L) {
                warning("The position of '", locus, "' may not be appropriate for building the model.", 
                  immediate. = TRUE)
            }
        }
    }
    else {
        if (!is.finite(pos.mid)) 
            stop("'pos.mid' should be specified.")
        pos.start <- pos.mid - flank.bp
        pos.end <- pos.mid + flank.bp
    }
    if (is.finite(pos.start) & is.finite(pos.end)) {
        flag <- (pos.start <= genoobj$snp.position) & (genoobj$snp.position <= 
            pos.end)
        flag[is.na(flag)] <- FALSE
        rv <- hlaGenoSubset(genoobj, snp.sel = flag)
    }
    else stop("No position information for ", locus, " on ", 
        assembly)
    rv
}

#######################################################################
# Summary a "hlaAlleleClass" object
#
#' @export
summary.hlaAlleleClass <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(inherits(object, "hlaAlleleClass"))
	hla <- object

	HUA <- hlaUniqueAllele(c(hla$value$allele1, hla$value$allele2))
	HLA <- factor(match(c(hla$value$allele1, hla$value$allele2), HUA))
	levels(HLA) <- HUA	
	count <- table(HLA)
	freq <- prop.table(count)
	rv <- cbind(count=count, freq=freq)

	# get the number of unique genotypes
	m <- data.frame(a1 = hla$value$allele1,
		a2 = hla$value$allele2, stringsAsFactors=FALSE)
	lst <- apply(m, 1, FUN=function(x) {
		if (!is.na(x[1]) & !is.na(x[2]))
		{
			if (x[1] <= x[2])
				paste(x[1], x[2], sep="/")
			else
				paste(x[2], x[1], sep="/")
		} else
			NA
	})
	unique.n.geno <- nlevels(factor(lst))

	if (show)
	{
		cat("Gene: ", .hla_gene_name_string(hla$locus), "\n", sep="")
		cat(sprintf("Range: [%dbp, %dbp]", hla$pos.start, hla$pos.end))
		if (!is.null(hla$assembly))
			cat(" on ", hla$assembly, "\n", sep="")
		else
			cat("\n")
		cat(sprintf("# of samples: %d\n", dim(hla$value)[1]))
		cat(sprintf("# of unique KIR3DL1/S1 alleles: %d\n", length(count)))
		cat(sprintf("# of unique KIR3DL1/S1 genotypes: %d\n", unique.n.geno))

		p <- hla$value$prob
		if (!is.null(p))
		{
			z <- table(cut(p, breaks=c(0, 0.25, 0.5, 0.75, 1),
				right=FALSE, include.lowest=TRUE))
			z[] <- sprintf("%d (%0.1f%%)", z, prop.table(z)*100)
			names(attr(z, "dimnames")) <- "Posterior probability:"
			print(z)
		}
	}

	# return
	return(invisible(rv))
}




##########################################################################
##########################################################################
#
# Attribute Bagging method -- HIBAG algorithm used for PONG
#

##########################################################################
# To fit an attribute bagging model for predicting
#

hlaAttrBagging <- function(hla, snp, nclassifier=100,
	mtry=c("sqrt", "all", "one"), prune=TRUE, rm.na=TRUE,
	verbose=TRUE, verbose.detail=FALSE)
{
	# check
	stopifnot(inherits(hla, "hlaAlleleClass"))
	stopifnot(inherits(snp, "hlaSNPGenoClass"))
	stopifnot(is.character(mtry) | is.numeric(mtry))
	stopifnot(is.logical(verbose))
	stopifnot(is.logical(verbose.detail))
	if (verbose.detail) verbose <- TRUE

	# get the common samples
	samp.id <- intersect(hla$value$sample.id, snp$sample.id)

	# hla types
	samp.flag <- match(samp.id, hla$value$sample.id)
	hla.allele1 <- hla$value$allele1[samp.flag]
	hla.allele2 <- hla$value$allele2[samp.flag]
	if (rm.na)
	{
		if (any(is.na(c(hla.allele1, hla.allele2))))
		{
			warning("There are missing KIR3DL1/S1 alleles, and the corresponding samples have been removed.")
			flag <- is.na(hla.allele1) | is.na(hla.allele2)
			samp.id <- setdiff(samp.id, hla$value$sample.id[samp.flag[flag]])
			samp.flag <- match(samp.id, hla$value$sample.id)
			hla.allele1 <- hla$value$allele1[samp.flag]
			hla.allele2 <- hla$value$allele2[samp.flag]
		}
	} else {
		if (any(is.na(c(hla.allele1, hla.allele2))))
		{
			stop("There are missing KIR3DL1/S1 alleles!")
		}
	}

	# SNP genotypes
	samp.flag <- match(samp.id, snp$sample.id)
	snp.geno <- snp$genotype[, samp.flag]
	storage.mode(snp.geno) <- "integer"

	tmp.snp.id <- snp$snp.id
	tmp.snp.position <- snp$snp.position
	tmp.snp.allele <- snp$snp.allele

	# remove mono-SNPs
	snpsel <- rowMeans(snp.geno, na.rm=TRUE)
	snpsel[!is.finite(snpsel)] <- 0
	snpsel <- (0 < snpsel) & (snpsel < 2)
	if (sum(!snpsel) > 0)
	{
		snp.geno <- snp.geno[snpsel, ]
		if (verbose)
		{
			a <- sum(!snpsel)
			cat(sprintf("%s monomorphic SNP%s ha%s been removed.\n",
				a, if (a>1) "s" else "", if (a>1) "ve" else "s"))
		}
		tmp.snp.id <- tmp.snp.id[snpsel]
		tmp.snp.position <- tmp.snp.position[snpsel]
		tmp.snp.allele <- tmp.snp.allele[snpsel]
	}

	if (length(samp.id) <= 0)
		stop("There is no common sample between 'hla' and 'snp'.")
	if (length(dim(snp.geno)[1]) <= 0)
		stop("There is no valid SNP markers.")


	###################################################################
	# initialize ...

	n.snp <- dim(snp.geno)[1]      # Num. of SNPs
	n.samp <- dim(snp.geno)[2]     # Num. of samples
	HUA <- hlaUniqueAllele(c(hla.allele1, hla.allele2))
	H <- factor(match(c(hla.allele1, hla.allele2), HUA))
	levels(H) <- HUA
	n.hla <- nlevels(H)
	H1 <- as.integer(H[1:n.samp]) - as.integer(1)
	H2 <- as.integer(H[(n.samp+1):(2*n.samp)]) - as.integer(1)

	# create an attribute bagging object
	rv <- .C("PONG_Training", n.snp, n.samp, snp.geno, n.hla,
		H1, H2, AB=integer(1), err=integer(1),
		NAOK = TRUE, PACKAGE = "PONG")
	if (rv$err != 0) stop(hlaErrMsg())
	ABmodel <- rv$AB

	# number of variables randomly sampled as candidates at each split
	mtry <- mtry[1]
	if (is.character(mtry))
	{
		if (mtry == "sqrt")
		{
			mtry <- ceiling(sqrt(n.snp))
		} else if (mtry == "all")
		{
			mtry <- n.snp
		} else if (mtry == "one")
		{
			mtry <- as.integer(1)
		} else {
			stop("Invalid mtry!")
		}
	} else if (is.numeric(mtry))
	{
		if (is.finite(mtry))
		{
			if ((0 < mtry) & (mtry < 1)) mtry <- n.snp*mtry
			mtry <- ceiling(mtry)
			if (mtry > n.snp) mtry <- n.snp
		} else {
			mtry <- ceiling(sqrt(n.snp))
		}
	} else {
		stop("Invalid mtry value!")
	}
	if (mtry <= 0) mtry <- as.integer(1)

	if (verbose)
	{
		cat("Build a PONG model with", nclassifier, "individual classifiers:\n")
		cat("# of SNPs randomly sampled as candidates for each selection: ",
			mtry, "\n", sep="")
		cat("# of SNPs: ", n.snp, ", # of samples: ", n.samp, "\n", sep="")
		cat("# of unique KIR3DL1/S1 alleles: ", n.hla, "\n", sep="")
	}


	###################################################################
	# training ...
	# add new individual classifers
	rv <- .C("PONG_NewClassifiers", ABmodel, as.integer(nclassifier),
		as.integer(mtry), as.logical(prune), verbose, verbose.detail,
		err=integer(1), NAOK = TRUE, PACKAGE = "PONG")
	if (rv$err != 0) stop(hlaErrMsg())

	# output
	rv <- list(n.samp = n.samp, n.snp = n.snp, sample.id = samp.id,
		snp.id = tmp.snp.id, snp.position = tmp.snp.position,
		snp.allele = tmp.snp.allele,
		snp.allele.freq = 0.5*rowMeans(snp.geno, na.rm=TRUE),
		hla.locus = hla$locus, hla.allele = levels(H),
		hla.freq = prop.table(table(H)),
		assembly = as.character(snp$assembly)[1],
		model = ABmodel,
		appendix = list())
	if (is.na(rv$assembly)) rv$assembly <- "unknown"

	class(rv) <- "hlaAttrBagClass"
	return(rv)
}


##########################################################################
# To fit an attribute bagging model for predicting
#

hlaParallelAttrBagging <- function(cl, hla, snp, auto.save="",
	nclassifier=100, mtry=c("sqrt", "all", "one"), prune=TRUE, rm.na=TRUE,
	stop.cluster=FALSE, verbose=TRUE)
{
	# check
	stopifnot(is.null(cl) | inherits(cl, "cluster"))
	stopifnot(inherits(hla, "hlaAlleleClass"))
	stopifnot(inherits(snp, "hlaSNPGenoClass"))

	stopifnot(is.character(auto.save) & (length(auto.save)==1))
	stopifnot(is.numeric(nclassifier))
	stopifnot(is.character(mtry) | is.numeric(mtry))
	stopifnot(is.logical(prune))
	stopifnot(is.logical(rm.na))
	stopifnot(is.logical(stop.cluster))
	stopifnot(is.logical(verbose))

	if (!is.null(cl))
	{
	    if (!require(parallel, warn.conflicts=FALSE))
			stop("The `parallel' package should be installed.")
	}

	if (verbose)
	{
		if (!is.null(cl))
		{
			cat(sprintf(
				"Build a PONG model of %d individual classifier%s in parallel with %d node%s:\n",
				nclassifier, if (nclassifier>1) "s" else "",
				length(cl), if (length(cl)>1) "s" else ""
			))
		} else {
			cat(sprintf(
				"Build a PONG model of %d individual classifier%s:\n",
				nclassifier, if (nclassifier>1) "s" else ""
			))
		}
		if (auto.save != "")
			cat("The model is autosaved in '", auto.save, "'.\n", sep="")
	}

	# set random number
	if (!is.null(cl))
	{
		RNGkind("L'Ecuyer-CMRG")
		rand <- .Random.seed
		parallel::clusterSetRNGStream(cl)
	}

	ans <- local({
		total <- 0

		.DynamicClusterCall(cl,
			fun = function(job, hla, snp, mtry, prune, rm.na)
			{
				library(PONG)
				model <- hlaAttrBagging(hla=hla, snp=snp, nclassifier=1,
					mtry=mtry, prune=prune, rm.na=rm.na,
					verbose=FALSE, verbose.detail=FALSE)
				mobj <- hlaModelToObj(model)
				hlaClose(model)
				return(mobj)
			},
			combine.fun = function(obj1, obj2)
			{
				if (is.null(obj1))
					mobj <- obj2
				else if (is.null(obj2))
					mobj <- obj1
				else
					mobj <- hlaCombineModelObj(obj1, obj2)
				if (auto.save != "")
					save(mobj, file=auto.save)
				if (verbose & !is.null(mobj))
				{
					z <- summary(mobj, show=FALSE)
					cat(sprintf("  --  average out-of-bag accuracy: %0.2f%%, sd: %0.2f%%, min: %0.2f%%, max: %0.2f%%\n",
						z$info["accuracy", "Mean"], z$info["accuracy", "SD"],
						z$info["accuracy", "Min"], z$info["accuracy", "Max"]))
				}
				mobj
			},
			msg.fn = function(job, obj)
			{
				if (verbose)
				{
					z <- summary(obj, show=FALSE)
					total <<- total + 1
					cat(date(), sprintf(
						", %4d, job %3d, # of SNPs: %g, # of haplotypes: %g, accuracy: %0.1f%%\n",
						total, as.integer(job), z$info["num.snp", "Mean"],
						z$info["num.haplo", "Mean"],
						z$info["accuracy", "Mean"]), sep="")
				}
			},
			n = nclassifier, stop.cluster = stop.cluster,
			hla=hla, snp=snp, mtry=mtry, prune=prune, rm.na=rm.na
		)
	})

	if (!is.null(cl) & !stop.cluster)
	{
		parallel::nextRNGStream(rand)
		parallel::nextRNGSubStream(rand)
	}

	# return
	if (auto.save == "")
		return(hlaModelFromObj(ans))
	else
		return(invisible())
}


##########################################################################
# To fit an attribute bagging model for predicting
#

hlaClose <- function(model)
{
	# check
	stopifnot(inherits(model, "hlaAttrBagClass"))

	# class handler
	rv <- .C("PONG_Close", model$model, err=integer(1),
		NAOK = TRUE, PACKAGE = "PONG")
	if (rv$err != 0) stop(hlaErrMsg())

	# output
	return(invisible(NULL))
}


#######################################################################
# Check missing SNP predictors
#

hlaCheckSNPs <- function(model, object,
	match.type=c("RefSNP+Position", "RefSNP", "Position"), verbose=TRUE)
{
	# check
	stopifnot(inherits(model, "hlaAttrBagClass") |
		inherits(model, "hlaAttrBagObj"))
	stopifnot((is.vector(object) & is.character(object)) |
		inherits(object, "hlaSNPGenoClass"))
	match.type <- match.arg(match.type)
	stopifnot(is.logical(verbose))

	# initialize
	if (inherits(model, "hlaAttrBagClass"))
		model <- hlaModelToObj(model)

	# show information
	if (verbose)
	{
		cat("The PONG model:\n")
		cat(sprintf("\tThere are %d SNP predictors in total.\n",
			length(model$snp.id)))
		cat(sprintf("\tThere are %d individual classifiers.\n",
			length(model$classifiers)))
	}

	if (is.vector(object))
	{
		target.snp <- as.character(object)
		src.snp <- hlaSNPID(model, match.type)
	} else {
		target.snp <- hlaSNPID(object, match.type)
		src.snp <- hlaSNPID(model, match.type)
	}

	NumOfSNP <- integer(length(model$classifiers))
	NumOfValidSNP <- integer(length(model$classifiers))

	# enumerate each classifier
	for (i in 1:length(model$classifiers))
	{
		v <- model$classifiers[[i]]
		flag <- src.snp[v$snpidx] %in% target.snp
		NumOfSNP[i] <- length(v$snpidx)
		NumOfValidSNP[i] <- sum(flag)
	}

	rv <- data.frame(NumOfValidSNP = NumOfValidSNP, NumOfSNP = NumOfSNP,
		fraction = NumOfValidSNP/NumOfSNP)

	if (verbose)
	{
		cat("Summarize the missing fractions of SNP predictors per classifier:\n")
		print(summary(1 - rv$fraction))
	}

	# output
	invisible(rv)
}


#######################################################################
# Predict KIR3DL1/S1 alleles from unphased SNP data
#

predict.hlaAttrBagClass <- function (object, snp, cl = FALSE, type = c("response+dosage", 
    "response", "prob", "response+prob"), vote = c("prob", "majority"), 
    allele.check = TRUE, match.type = c("Position", "Pos+Allele", 
        "RefSNP+Position", "RefSNP"), same.strand = FALSE, verbose = TRUE, 
    verbose.match = TRUE, ...) 
{
    stopifnot(inherits(object, "hlaAttrBagClass"))
    hlaPredict(object, snp, cl, type, vote, allele.check, match.type, 
        same.strand, verbose, verbose.match)
}

kirPredict <- function (object, snp, cl = FALSE, type = c("response+dosage", 
    "response", "prob", "response+prob"), vote = c("prob", "majority"), 
    allele.check = TRUE, match.type = c("Position", "Pos+Allele", 
        "RefSNP+Position", "RefSNP"), same.strand = FALSE, verbose = TRUE, 
    verbose.match = TRUE) 
{
    stopifnot(inherits(object, "hlaAttrBagClass"))
    stopifnot(is.logical(cl) | is.numeric(cl) | inherits(cl, 
        "cluster"))
    stopifnot(is.logical(allele.check), length(allele.check) == 
        1L)
    stopifnot(is.logical(same.strand), length(same.strand) == 
        1L)
    stopifnot(is.logical(verbose), length(verbose) == 1L)
    stopifnot(is.logical(verbose.match), length(verbose.match) == 
        1L)
    type <- match.arg(type)
    vote <- match.arg(vote)
    match.type <- match.arg(match.type)
    vote_method <- match(vote, c("prob", "majority"))
    if (inherits(cl, "cluster")) {
        if (!requireNamespace("parallel", quietly = TRUE)) 
            stop("The `parallel' package should be installed.")
        if (length(cl) <= 1L) 
            cl <- FALSE
    }
    if (!is.null(object$appendix$warning)) 
        warning(object$appendix$warning, immediate. = TRUE)
    if (verbose) {
        #CNum <- .Call(HIBAG_GetNumClassifiers, object$model)
		CNum <- .C("PONG_GetNumClassifiers", object$model, CNum = integer(1),
			err=integer(1), NAOK=TRUE, PACKAGE="PONG")

        s <- object$hla.allele
        if (length(s) > 3L) 
            s <- c(s[1:3], "...")
        cat("PONG model for ", .hla_gene_name_string(object$hla.locus), 
            ":\n", "    ", CNum$CNum, " individual classifier(s)", 
            "\n", "    ", length(object$snp.id), " SNPs\n", "    ", 
            length(object$hla.allele), " unique HLA alleles: ", 
            paste(s, collapse = ", "), "\n", sep = "")
        
		cat("Prediction:\n")
        if (vote_method == 1L) 
            cat("    based on the averaged posterior probabilities\n")
        else cat("    by voting from all individual classifiers\n")
        if (inherits(cl, "cluster")) {
            cat("    run in parallel with ", length(cl), " compute node", 
                length(cl), "\n", sep = "")
        }
    }
    if (!inherits(snp, "hlaSNPGenoClass")) {
        stopifnot(is.numeric(snp), is.vector(snp) | is.matrix(snp))
        if (is.vector(snp)) {
            stopifnot(length(snp) == object$n.snp)
            snp <- matrix(snp, ncol = 1L)
        }
        else {
            stopifnot(nrow(snp) == object$n.snp)
        }
        geno.sampid <- 1L:ncol(snp)
        assembly <- "auto-silent"
    }
    else {
        model.assembly <- as.character(object$assembly)[1L]
        if (is.na(model.assembly)) 
            model.assembly <- "unknown"
        geno.assembly <- as.character(snp$assembly)[1L]
        if (is.na(geno.assembly)) 
            geno.assembly <- "unknown"
        refstr <- sprintf("Model assembly: %s, SNP assembly: %s", 
            model.assembly, geno.assembly)
        if (verbose) 
            cat(refstr, "\n", sep = "")
        if (model.assembly != geno.assembly) {
            if (any(c(model.assembly, geno.assembly) %in% "unknown")) {
                if (verbose) 
                  message("The human genome references might not match!")
                if (geno.assembly == "unknown") 
                  assembly <- model.assembly
                else assembly <- geno.assembly
            }
            else {
                warning("The human genome references do not match! ", 
                  refstr, ".", immediate. = TRUE)
                assembly <- model.assembly
            }
        }
        else {
            if (model.assembly != "unknown") 
                assembly <- model.assembly
            else assembly <- "auto"
        }
        if (verbose && verbose.match) {
            cat("Matching the SNPs between the model and the test data:\n")
            tab <- NULL
            for (tp in c("Position", "Pos+Allele", "RefSNP+Position", 
                "RefSNP")) {
                s1 <- hlaSNPID(object, tp)
                s2 <- hlaSNPID(snp, tp)
                s <- unique(intersect(s1, s2))
                mcnt <- length(s1) - length(s)
                d <- data.frame(c1 = paste0("   ", tp), c2 = sprintf("%d (%.1f%%)", 
                  mcnt, mcnt/length(s1) * 100), c3 = ifelse(tp == 
                  match.type, "*being used", ""), stringsAsFactors = FALSE)
                tab <- rbind(tab, d)
            }
            names(tab) <- c("match.type=\"--\"", "  missing SNPs #", 
                "")
            tab[1L, 3L] <- paste(tab[1L, 3L], "[1]")
            tab[2L, 3L] <- paste(tab[2L, 3L], "[2]")
            print(tab, row.names = FALSE)
            cat("      [1]: useful if ambiguous strands on array-based platforms\n")
            cat("      [2]: suggested if the model and test data have been matched to the same reference genome\n")
            s <- object$appendix$platform
            if (is.null(s)) 
                s <- "not applicable"
            else s <- paste(s, collapse = ",")
            cat("    Model platform: ", s, "\n", sep = "")
        }
        else if (verbose) {
            cat("Using match.type='", match.type, "' for SNP matching\n", 
                sep = "")
        }
        geno.sampid <- snp$sample.id
        obj.id <- hlaSNPID(object, match.type)
        geno.id <- hlaSNPID(snp, match.type)
        flag <- FALSE
        if (length(obj.id) == length(geno.id)) {
            if (all(obj.id == geno.id)) 
                flag <- TRUE
        }
        if (flag) {
            if (allele.check) {
                snp <- hlaGenoSwitchStrand(snp, object, match.type, 
                  same.strand, verbose)$genotype
            }
            else {
                snp <- snp$genotype
            }
        }
        else {

			#print(match(obj.id, geno.id))
            snp.sel <- match(obj.id, geno.id)
            snp.sel[duplicated(snp.sel)] <- NA_integer_
            snp.allele <- snp$snp.allele
            snp.allele[is.na(snp.allele)] <- ""
            tmp <- list(genotype = snp$genotype[snp.sel, , drop = FALSE], 
                sample.id = snp$sample.id, snp.id = object$snp.id, 
                snp.position = object$snp.position, snp.allele = snp.allele[snp.sel], 
                assembly = snp$assembly)
            flag <- is.na(tmp$snp.allele)
            tmp$snp.allele[flag] <- object$snp.allele[match(tmp$snp.id[flag], 
                object$snp.id)]
            class(tmp) <- "hlaSNPGenoClass"
            missing.cnt <- sum(flag)
            if (missing.cnt == length(obj.id)) {
                stop("There is no overlapping of SNPs!")
            }
            else if (missing.cnt > 0.5 * length(obj.id)) {
                warning("More than 50% of SNPs are missing!", 
                  immediate. = TRUE)
            }
            if (allele.check) {
                snp <- hlaGenoSwitchStrand(tmp, object, match.type, 
                  same.strand, verbose)$genotype
            }
            else {
                snp <- tmp$genotype
            }
        }
    }
    if (dim(snp)[1L] != object$n.snp) {
        stop("The number of SNPs is not valid, and it maybe due to duplicated 'snp.id' ", 
            "or incorrect dimension of genotype matrix.")
    }
    n.samp <- dim(snp)[2L]
    n.hla <- length(object$hla.allele)
    if (verbose) {
        cat(sprintf("# of samples: %d\n", n.samp))
		#Workaround for the following error: Error in .Call(HIBAG_Predict_Resp, object$model, as.integer(snp), n.samp,  :
		#  "PONG_Kernel_Version" not available for .Call() for package "PONG"
        #cat("CPU flags: ", .C("PONG_Kernel_Version", PACKAGE="PONG")[[2L]][1L], "\n", sep = "")
    }
	
    if (!inherits(cl, "cluster")) {
        nthread <- 1L
        if (isTRUE(cl)) 
            nthread <- as.integer(defaultNumThreads())
        if (is.numeric(cl)) 
            nthread <- cl[1L]
        if (is.na(nthread) || nthread < 1L) 
            nthread <- 1L
        pm <- attr(cl, "proc_ptr")
        if (type %in% c("response", "response+dosage", "response+prob")) {
            if (type == "response") 
			{
                rv <- .C("PONG_Predict_Resp", object$model, as.integer(snp),
					n.samp, as.integer(vote_method), as.logical(verbose),
					H1=integer(n.samp), H2=integer(n.samp), prob=double(n.samp),
					err=integer(1), NAOK=TRUE, PACKAGE="PONG")
            }
            else if (type == "response+dosage") {
                rv <- .C("PONG_Predict_Resp_dosage", object$model, as.integer(snp),
					n.samp, as.integer(vote_method), as.logical(verbose),
					H1=integer(n.samp), H2=integer(n.samp), prob=double(n.samp), 
					postprob=matrix(NaN, nrow=n.hla*(n.hla+1)/2, ncol=n.samp),
					err=integer(1), NAOK=TRUE, PACKAGE="PONG")
            }
            else {
                rv <- .C("PONG_Predict_Resp_Prob", object$model, as.integer(snp),
					n.samp, as.integer(vote_method), as.logical(verbose),
					H1=integer(n.samp), H2=integer(n.samp), prob=double(n.samp), 
					postprob=matrix(NaN, nrow=n.hla*(n.hla+1)/2, ncol=n.samp),
					err=integer(1), NAOK=TRUE, PACKAGE="PONG")
            }
            res <- hlaAllele(geno.sampid, H1 = object$hla.allele[rv$H1 + 
                1L], H2 = object$hla.allele[rv$H2 + 1L], locus = object$hla.locus, 
                prob = rv$prob, na.rm = FALSE, assembly = assembly)
            res$value$matching <- rv$matching
            if (!is.null(rv$dosage)) {
                res$dosage <- rv$dosage
                rownames(res$dosage) <- object$hla.allele
                colnames(res$dosage) <- geno.sampid
            }
            if (!is.null(rv$postprob)) {
                res$postprob <- rv$postprob
                colnames(res$postprob) <- geno.sampid
                m <- outer(object$hla.allele, object$hla.allele, 
                  paste, sep = "/")
                rownames(res$postprob) <- m[lower.tri(m, diag = TRUE)]
            }
            NA.cnt <- sum(is.na(res$value$allele1) | is.na(res$value$allele2))
        }
        else {
            rv <- .Call(HIBAG_Predict_Resp_Prob, object$model, 
                as.integer(snp), n.samp, vote_method, nthread, 
                verbose, pm)
            names(rv) <- c("H1", "H2", "prob", "matching", "dosage", 
                "postprob")
            res <- rv$postprob
            colnames(res) <- geno.sampid
            m <- outer(object$hla.allele, object$hla.allele, 
                paste, sep = "/")
            rownames(res) <- m[lower.tri(m, diag = TRUE)]
            NA.cnt <- sum(colSums(res) <= 0L, na.rm = TRUE)
        }
    }
    else {
        rv <- parallel::clusterApply(cl = cl, parallel::splitIndices(n.samp, 
            length(cl)), fun = function(i, mobj, snp, type, vote) {
            if (length(i) > 0L) {
                library(HIBAG)
                m <- hlaModelFromObj(mobj)
                on.exit(hlaClose(m))
                hlaPredict(m, snp[, i], type = type, vote = vote, 
                  verbose = FALSE)
            }
            else NULL
        }, mobj = hlaModelToObj(object), snp = snp, type = type, 
            vote = vote)
        if (type %in% c("response", "response+dosage", "response+prob")) {
            res <- rv[[1L]]
            for (i in 2L:length(rv)) {
                if (!is.null(rv[[i]])) 
                  res <- hlaCombineAllele(res, rv[[i]])
            }
            res$value$sample.id <- geno.sampid
            if (!is.null(res$dosage)) 
                colnames(res$dosage) <- geno.sampid
            if (!is.null(res$postprob)) 
                colnames(res$postprob) <- geno.sampid
            NA.cnt <- sum(is.na(res$value$allele1) | is.na(res$value$allele2))
        }
        else {
            res <- rv[[1L]]
            for (i in 2L:length(rv)) {
                if (!is.null(rv[[i]])) 
                  res <- cbind(res, rv[[i]])
            }
            colnames(res) <- geno.sampid
            NA.cnt <- sum(colSums(res) <= 0L, na.rm = TRUE)
        }
    }
    if (NA.cnt > 0L) {
        warning("No prediction output for ", NA.cnt, " individual", 
            .plural(NA.cnt), " (possibly due to missing SNPs).", 
            immediate. = TRUE)
    }
    res
}


#######################################################################
# Merge predictions by voting
#

hlaPredMerge <- function(..., weight=NULL, equivalence=NULL)
{
	# check "..."
	pdlist <- list(...)
	if (length(pdlist) <= 0)
		stop("No object is passed to 'hlaPredMerge'.")
	for (i in 1:length(pdlist))
	{
		if (!inherits(pdlist[[i]], "hlaAlleleClass"))
			stop("The object(s) passed to 'hlaPredMerge' should be 'hlaAlleleClass'.")
		if (is.null(pdlist[[i]]$postprob))
		{
			stop("The object(s) passed to 'hlaPredMerge' should have a field of 'postprob',",
				" returned by 'predict(..., type=\"response+prob\", vote=\"majority\")'.")
		}
	}

	# check equivalence
	stopifnot(is.null(equivalence) | is.data.frame(equivalence))
	if (is.data.frame(equivalence))
	{
		if (ncol(equivalence) != 2)
		{
			stop("'equivalence' should have two columns: ",
				"the first for new equivalent alleles, and the second for the ",
				"alleles possibly existed in the object(s) passed to 'hlaPredMerge'.")
		}
	}

	# check locus and sample.id
	samp.id <- pdlist[[1]]$value$sample.id
	locus <- pdlist[[1]]$locus
	for (i in 1:length(pdlist))
	{
		if (!identical(samp.id, pdlist[[i]]$value$sample.id))
			stop("The sample IDs should be the same.")
		if (!identical(locus, pdlist[[i]]$locus))
			stop("The locus should be the same.")
	}

	# check weight
	if (!is.null(weight))
	{
		stopifnot(is.numeric(weight) & is.vector(weight))
		if (length(pdlist) != length(weight))
			stop("Invalid 'weight'.")
		weight <- abs(weight)
		weight <- weight / sum(weight)
	} else {
		weight <- rep(1/length(pdlist), length(pdlist))
	}

	#############################################################
	# replace function
	#

	replace <- function(allele)
	{
		if (!is.null(equivalence))
		{
			i <- match(allele, equivalence[, 2])
			flag <- !is.na(i)
			i <- i[flag]
			allele[flag] <- equivalence[i, 1]
		}
		allele
	}


	# all different alleles
	hla.allele <- NULL
	for (i in 1:length(pdlist))
	{
		h <- unique(unlist(strsplit(rownames(pdlist[[i]]$postprob), "/")))
		hla.allele <- unique(c(hla.allele, replace(h)))
	}
	hla.allele <- hlaUniqueAllele(hla.allele)
	n.hla <- length(hla.allele)
	n.samp <- length(samp.id)

	prob <- matrix(0.0, nrow=n.hla*(n.hla+1)/2, ncol=n.samp)
	m <- outer(hla.allele, hla.allele, function(x, y) paste(x, y, sep="/"))
	m <- m[lower.tri(m, diag=TRUE)]

	# for-loop
	for (i in 1:length(pdlist))
	{
		p <- pdlist[[i]]$postprob
		h <- replace(unlist(strsplit(rownames(p), "/")))

		h1 <- h[seq(1, length(h), 2)]; h2 <- h[seq(2, length(h), 2)]
		j1 <- match(paste(h1, h2, sep="/"), m)
		j2 <- match(paste(h2, h1, sep="/"), m)
		j1[is.na(j1)] <- j2[is.na(j1)]

		# check
		stopifnot(!any(is.na(j1)))

		p <- p * weight[i]
		for (j in seq_len(length(j1)))
			prob[j1[j], ] <- prob[j1[j], ] + p[j, ]
	}
	colnames(prob) <- samp.id
	rownames(prob) <- m

	pb <- apply(prob, 2, max)
	pt <- unlist(strsplit(m[apply(prob, 2, which.max)], "/"))
	assembly <- pdlist[[1]]$assembly
	if (is.null(assembly)) assembly <- "auto"

	rv <- hlaAllele(samp.id,
		H1 = pt[seq(2, length(pt), 2)],
		H2 = pt[seq(1, length(pt), 2)],
		locus = locus,
		locus.pos.start = pdlist[[1]]$pos.start,
		locus.pos.end = pdlist[[1]]$pos.end,
		prob = pb, na.rm = FALSE,
		assembly = assembly)
	rv$postprob <- prob
	rv
}




#######################################################################
# summarize the "hlaAttrBagClass" object
#
#' @export
summary.hlaAttrBagClass <- function(object, show=TRUE, ...)
{
	obj <- hlaModelToObj(object)
	summary(obj, show=show)
}


#######################################################################
# Save the parameters in a model of attribute bagging
#

hlaModelToObj <- function(model)
{
	# check
	stopifnot(inherits(model, "hlaAttrBagClass"))

	# call, get the number of classifiers
	rv <- .C("PONG_GetNumClassifiers", model$model, CNum = integer(1),
		err=integer(1), NAOK=TRUE, PACKAGE="PONG")
	if (rv$err != 0) stop(hlaErrMsg())

	# for each tree
	res <- vector("list", rv$CNum)
	for (i in 1:length(res))
	{
		# call, get the number of haplotypes
		rv <- .C("PONG_Idv_GetNumHaplo", model$model, as.integer(i),
			NumHaplo = integer(1), NumSNP = integer(1),
			err=integer(1), NAOK=TRUE, PACKAGE="PONG")
		if (rv$err != 0) stop(hlaErrMsg())

		# call, get freq. and haplotypes
		rv <- .C("PONG_Classifier_GetHaplos", model$model, as.integer(i),
			freq=double(rv$NumHaplo), hla=integer(rv$NumHaplo),
			haplo=character(rv$NumHaplo), snpidx = integer(rv$NumSNP),
			samp.num = integer(model$n.samp), acc = double(1),
			err=integer(1), NAOK=TRUE, PACKAGE="PONG")
		if (rv$err != 0) stop(hlaErrMsg())

		res[[i]] <- list(
			samp.num = rv$samp.num,
			haplos = data.frame(freq = rv$freq, hla = model$hla.allele[rv$hla],
				haplo = rv$haplo, stringsAsFactors=FALSE),
			snpidx = rv$snpidx,
			outofbag.acc = rv$acc)
	}

	rv <- list(n.samp = model$n.samp, n.snp = model$n.snp,
		sample.id = model$sample.id, snp.id = model$snp.id,
		snp.position = model$snp.position, snp.allele = model$snp.allele,
		snp.allele.freq = model$snp.allele.freq,
		hla.locus = model$hla.locus,
		hla.allele = model$hla.allele, hla.freq = model$hla.freq,
		assembly = model$assembly,
		classifiers = res,
		appendix <- model$appendix)
	class(rv) <- "hlaAttrBagObj"
	return(rv)
}


#######################################################################
# To combine two model objects of attribute bagging
#

hlaCombineModelObj <- function(obj1, obj2)
{
	# check
	stopifnot(inherits(obj1, "hlaAttrBagObj"))
	stopifnot(inherits(obj2, "hlaAttrBagObj"))
	stopifnot(identical(obj1$hla.locus, obj2$hla.locus))
	stopifnot(identical(obj1$snp.id, obj2$snp.id))
	stopifnot(identical(obj1$hla.allele, obj2$hla.allele))
	stopifnot(identical(obj1$assembly, obj2$assembly))

	samp.id <- unique(c(obj1$sample.id, obj2$sample.id))
	if (!is.null(obj1$appendix) | !is.null(obj2$appendix))
	{
		appendix <- list(
			platform =
				unique(c(obj1$appendix$platform, obj2$appendix$platform)),
			information =
				unique(c(obj1$appendix$information, obj2$appendix$information)),
			warning =
				unique(c(obj1$appendix$warning, obj2$appendix$warning))
		)
	} else {
		appendix <- NULL
	}

	rv <- list(n.samp = length(samp.id), n.snp = obj1$n.snp,
		sample.id = samp.id, snp.id = obj1$snp.id,
		snp.position = obj1$snp.position, snp.allele = obj1$snp.allele,
		snp.allele.freq = (obj1$snp.allele.freq + obj2$snp.allele.freq)*0.5,
		hla.locus = obj1$hla.locus, hla.allele = obj1$hla.allele,
		hla.freq = (obj1$hla.freq + obj2$hla.freq)*0.5,
		assembly = obj1$assembly,
		classifiers = c(obj1$classifiers, obj2$classifiers),
		appendix = appendix)
	class(rv) <- "hlaAttrBagObj"
	return(rv)
}


#######################################################################
# To get the top n individual classifiers
#

hlaSubModelObj <- function(obj, n)
{
	# check
	stopifnot(inherits(obj, "hlaAttrBagObj"))
	obj$classifiers <- obj$classifiers[1:n]
	return(obj)
}


#######################################################################
# To get a "hlaAttrBagClass" class
#

hlaModelFromObj <- function(obj)
{
	# check
	stopifnot(inherits(obj, "hlaAttrBagObj"))

	# create an attribute bagging object
	rv <- .C("PONG_New",
		as.integer(obj$n.samp), as.integer(obj$n.snp), length(obj$hla.allele),
		model = integer(1), err=integer(1), NAOK=TRUE, PACKAGE="PONG")
	if (rv$err != 0) stop(hlaErrMsg())
	ABmodel <- rv$model

	# add individual classifiers
	for (tree in obj$classifiers)
	{
		hla <- match(tree$haplos$hla, obj$hla.allele) - 1
		if (any(is.na(hla)))
			stop("Invalid KIR3DL1/S1 alleles in the individual classifier.")
		if (is.null(tree$samp.num))
			snum <- rep(as.integer(1), obj$n.samp)
		else
			snum <- tree$samp.num
		rv <- .C("PONG_NewClassifierHaplo", ABmodel, length(tree$snpidx),
			as.integer(tree$snpidx-1), as.integer(snum), dim(tree$haplos)[1],
			as.double(tree$haplos$freq), as.integer(hla),
			as.character(tree$haplos$haplo), as.double(tree$outofbag.acc),
			err=integer(1), NAOK=TRUE, PACKAGE="PONG")
		if (rv$err != 0) stop(hlaErrMsg())
	}

	# output
	rv <- list(n.samp = obj$n.samp, n.snp = obj$n.snp,
		sample.id = obj$sample.id, snp.id = obj$snp.id,
		snp.position = obj$snp.position, snp.allele = obj$snp.allele,
		snp.allele.freq = obj$snp.allele.freq,
		hla.locus = obj$hla.locus, hla.allele = obj$hla.allele,
		hla.freq = obj$hla.freq,
		assembly = as.character(obj$assembly)[1],
		model = ABmodel,
		appendix = obj$appendix)
	if (is.na(rv$assembly)) rv$assembly <- "unknown"

	class(rv) <- "hlaAttrBagClass"
	return(rv)
}


#######################################################################
# summarize the "hlaAttrBagObj" object
#
#' @export
summary.hlaAttrBagObj <- function(object, show=TRUE, ...)
{
	# check
	stopifnot(inherits(object, "hlaAttrBagObj"))
	obj <- object

	if (show)
	{
		cat("Gene: ", .hla_gene_name_string(obj$hla.locus), "\n", sep="")
		cat("Training dataset:", obj$n.samp, "samples X",
			length(obj$snp.id), "SNPs\n")
		cat("\t# of KIR3DL1/S1 alleles: ", length(obj$hla.allele), "\n", sep="")
	}

	# summarize ...
	snpset <- NULL
	outofbag.acc <- rep(NaN, length(obj$classifiers))
	numsnp <- rep(NA, length(obj$classifiers))
	numhaplo <- rep(NA, length(obj$classifiers))
	snp.hist <- rep(0, length(obj$snp.id))
	for (i in 1:length(obj$classifiers))
	{
		outofbag.acc[i] <- obj$classifiers[[i]]$outofbag.acc
		numsnp[i] <- length(obj$classifiers[[i]]$snpidx)
		numhaplo[i] <- length(obj$classifiers[[i]]$haplos$hla)
		snp.hist[obj$classifiers[[i]]$snpidx] <- snp.hist[obj$classifiers[[i]]$snpidx] + 1
		snpset <- unique(c(snpset, obj$classifiers[[i]]$snpidx))
	}
	snpset <- snpset[order(snpset)]
	outofbag.acc <- outofbag.acc * 100

	info <- data.frame(
		Mean = c(mean(numsnp), mean(numhaplo), mean(outofbag.acc)),
		SD = c(sd(numsnp), sd(numhaplo), sd(outofbag.acc)),
		Min = c(min(numsnp), min(numhaplo), min(outofbag.acc)),
		Max = c(max(numsnp), max(numhaplo), max(outofbag.acc))
	)
	rownames(info) <- c("num.snp", "num.haplo", "accuracy")

	if (show)
	{
		cat("\t# of individual classifiers: ", length(obj$classifiers), "\n", sep="")
		cat("\ttotal # of SNPs used: ", length(snpset), "\n", sep="")
		cat(sprintf(
			"\taverage # of SNPs in an individual classifier: %0.2f, sd: %0.2f, min: %d, max: %d\n",
			mean(numsnp), sd(numsnp), min(numsnp), max(numsnp)))
		cat(sprintf(
			"\taverage # of haplotypes in an individual classifier: %0.2f, sd: %0.2f, min: %d, max: %d\n",
			mean(numhaplo), sd(numhaplo), min(numhaplo), max(numhaplo)))
		cat(sprintf(
			"\taverage out-of-bag accuracy: %0.2f%%, sd: %0.2f%%, min: %0.2f%%, max: %0.2f%%\n",
			mean(outofbag.acc), sd(outofbag.acc), min(outofbag.acc),
			max(outofbag.acc)))

		if (is.null(obj$assembly))
			cat("Genome assembly: unknown\n")
		else
			cat("Genome assembly: ", obj$assembly, "\n", sep="")

		if (!is.null(obj$appendix$platform))
			cat("Platform:", obj$appendix$platform, "\n")
		if (!is.null(obj$appendix$information))
			cat("Information:", obj$appendix$information, "\n")
		if (!is.null(obj$appendix$warning))
			message(obj$appendix$warning)
	}

	rv <- list(num.classifier = length(obj$classifiers),
		num.snp = length(snpset),
		snp.id = obj$snp.id, snp.position = obj$snp.position,
		snp.hist = snp.hist, info = info)
	return(invisible(rv))
}




##########################################################################
# to finalize the PONG model
#

hlaPublish <- function(mobj, platform=NULL, information=NULL, warning=NULL,
	rm.unused.snp=TRUE, anonymize=TRUE, verbose=TRUE)
{
	# check
	stopifnot(inherits(mobj, "hlaAttrBagObj") |
		inherits(mobj, "hlaAttrBagClass"))
	stopifnot(is.null(platform) | is.character(platform))
	stopifnot(is.null(information) | is.character(information))
	stopifnot(is.null(warning) | is.character(warning))
	stopifnot(is.logical(rm.unused.snp) & (length(rm.unused.snp)==1))
	stopifnot(is.logical(verbose))
	if (inherits(mobj, "hlaAttrBagClass"))
		mobj <- hlaModelToObj(mobj)

	# additional information
	if (is.null(platform))
		platform <- mobj$appendix$platform
	if (is.null(information))
		information <- mobj$appendix$information
	if (is.null(warning))
		warning <- mobj$appendix$warning
	mobj$appendix <- list(
		platform=platform, information=information, warning=warning)

	# remove unused SNPs
	if (rm.unused.snp)
	{
		# get frequency of use for SNPs
		snp.hist <- rep(0, length(mobj$snp.id))
		for (i in 1:length(mobj$classifiers))
		{
			idx <- mobj$classifiers[[i]]$snpidx
			snp.hist[idx] <- snp.hist[idx] + 1
		}

		flag <- (snp.hist > 0)
		if (sum(flag) < mobj$n.snp)
		{
			if (verbose)
			{
				cnt <- mobj$n.snp - sum(flag)
				cat(sprintf("Remove %d unused SNP%s.\n", cnt, if (cnt>1) "s" else ""))
			}
		}
		mobj$n.snp <- sum(flag)
		mobj$snp.id <- mobj$snp.id[flag]
		mobj$snp.position <- mobj$snp.position[flag]
		mobj$snp.allele <- mobj$snp.allele[flag]
		mobj$snp.allele.freq <- mobj$snp.allele.freq[flag]

		idx.list <- rep(0, length(flag))
		idx.list[flag] <- 1:mobj$n.snp

		for (i in 1:length(mobj$classifiers))
		{
			mobj$classifiers[[i]]$snpidx <-
				idx.list[ mobj$classifiers[[i]]$snpidx ]
		}
	}

	# anonymize
	if (anonymize)
	{
		mobj$sample.id <- NULL
		for (i in 1:length(mobj$classifiers))
		{
			mobj$classifiers[[i]]$samp.num <- NULL
		}
	}

	# output
	return(mobj)
}


##########################################################################
# to get a model object of attribute bagging from a list of files
#

hlaModelFiles <- function(fn.list, action.missingfile=c("ignore", "stop"),
	verbose=TRUE)
{
	# check
	stopifnot(is.character(fn.list))
	stopifnot(is.logical(verbose))
	action.missingfile <- match.arg(action.missingfile)

	# for-loop
	rv <- NULL
	for (fn in fn.list)
	{
		if (file.exists(fn))
		{
			tmp <- get(load(fn))
			if (is.null(rv))
			{
				rv <- tmp
			} else {
				rv <- hlaCombineModelObj(rv, tmp)
			}
		} else {
			s <- sprintf("There is no '%s'.", fn)
			if (action.missingfile == "stop")
			{
				stop(s)
			} else {
				if (verbose) message(s)
			}
		}
	}
	rv
}


##########################################################################
# Out-of-bag estimation of overall accuracy, per-allele sensitivity, etc
#

hlaOutOfBag <- function(model, hla, snp, call.threshold=NaN, verbose=TRUE)
{
	# check
	stopifnot(inherits(model, "hlaAttrBagObj") |
		inherits(model, "hlaAttrBagClass"))
	stopifnot(inherits(hla, "hlaAlleleClass"))
	stopifnot(inherits(snp, "hlaSNPGenoClass"))

	stopifnot(is.numeric(call.threshold) & is.vector(call.threshold))
	stopifnot(length(call.threshold) == 1)

	stopifnot(is.logical(verbose) & is.vector(verbose))
	stopifnot(length(verbose) == 1)


	######################################################
	# initialize ...
	if (inherits(model, "hlaAttrBagClass"))
	{
		model <- hlaModelToObj(model)
		if (verbose) print(model)
	}

	# map samples
	if (is.null(model$sample.id))
		stop("There is no sample ID in the model.")
	samp.idx <- match(model$sample.id, snp$sample.id)
	if (any(is.na(samp.idx)))
		stop("Some of sample.id in the model do not exist in SNP genotypes.")
	hla.samp.idx <- match(model$sample.id, hla$value$sample.id)
	if (any(is.na(hla.samp.idx)))
		stop("Some of sample.id in the model do not exist in HLA types.")

	# map SNPs
	snp.idx <- match(model$snp.id, snp$snp.id)
	if (any(is.na(snp.idx)))
		stop("Some of snp.id in the model do not exist in SNP genotypes.")

	# genotypes and the number of classifiers
	geno <- snp$genotype[snp.idx, samp.idx]
	nclass <- length(model$classifiers)

	# the returned value
	ans <- NULL

	# the column names of details
	nm1 <- c("allele", "train.num", "train.freq")
	nm2 <- c("call.rate", "accuracy", "sensitivity", "specificity", "ppv", "npv")

	# for-loop
	for (i in 1:nclass)
	{
		mx <- model
		mx$classifiers <- mx$classifiers[i]
		s <- mx$classifiers[[1]]$samp.num
		if (is.null(s))
			stop("There is no bootstrap sample index.")

		tmp.model <- hlaModelFromObj(mx)
		tmp.geno <- geno[, s == 0]
		v <- predict(tmp.model, tmp.geno, verbose=FALSE)
		hlaClose(tmp.model)
		v$value$sample.id <- mx$sample.id[s == 0]
		pam <- hlaCompareAllele(hla, v, allele.limit=mx,
			call.threshold=call.threshold, verbose=FALSE)

		if (!is.null(ans))
		{
			# overall
			ans$overall <- ans$overall + pam$overall

			# confusion matrix
			ans$confusion <- ans$confusion + pam$confusion

			# details
			pam$detail <- pam$detail[, nm2]
			ans$n.detail <- ans$n.detail + !is.na(pam$detail)
			pam$detail[is.na(pam$detail)] <- 0
			ans$detail <- ans$detail + pam$detail
		} else {
			ans <- pam
			ans$detailhead <- ans$detail[, nm1]
			colnames(ans$detailhead) <- c("allele", "valid.num", "valid.freq")
			ans$detail <- ans$detail[, nm2]
			ans$n.detail <- !is.na(ans$detail)
			ans$detail[is.na(ans$detail)] <- 0
		}

		if (verbose)
		{
			cat(date(), sprintf(", passing the %d/%d classifiers.\n",
				i, nclass), sep="")
		}
	}

	# average
	ans$overall <- ans$overall / nclass
	ans$confusion <- ans$confusion / nclass
	ans$detail <- ans$detail / ans$n.detail

	# get miscall
	rv <- ans$confusion; diag(rv) <- 0
	m.max <- apply(rv, 2, max); m.idx <- apply(rv, 2, which.max)
	s <- rownames(ans$confusion)[m.idx]; s[m.max<=0] <- NA
	p <- m.max / apply(rv, 2, sum)

	# output
	ans$detail <- cbind(ans$detailhead, ans$detail,
		miscall=s, miscall.prop=p, stringsAsFactors=FALSE)
	ans$detailhead <- NULL
	ans$n.detail <- NULL
	ans
}


##########################################################################
# to create a report for evaluating accuracies
#

hlaReport <- function(object, export.fn="", type=c("txt", "tex", "html"),
	header=TRUE)
{
	# check
	stopifnot(is.list(object))
	stopifnot(is.data.frame(object$detail))
	stopifnot(is.character(export.fn))
	type <- match.arg(type)
	stopifnot(is.logical(header))

	# create an output file
	if (export.fn != "")
	{
		f <- file(export.fn, "wt")
		on.exit(close(f))
	} else {
		f <- ""
	}

	###############################################################
	# text format

	d <- data.frame(Allele=object$detail$allele, stringsAsFactors=FALSE)

	if (!is.null(object$detail$train.num))
		d$NumTrain <- object$detail$train.num
	if (!is.null(object$detail$train.freq))
	{
		d$FreqTrain <- sprintf("%0.4f", object$detail$train.freq)
		d$FreqTrain[d$NumTrain <= 0] <- "0"

		L1 <- c("Allele", "Num.", "Freq.", "Num.", "Freq.",
			"CR", "ACC", "SEN", "SPE", "PPV", "NPV", "Miscall")
		L2 <- c("", "Train", "Train", "Valid.", "Valid.",
			"(%)", "(%)", "(%)", "(%)", "(%)", "(%)", "(%)")
		L2a <- c("", "Train", "Train", "Valid.", "Valid.",
			"(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)")
	} else {
		L1 <- c("Allele", "Num.", "Freq.",
			"CR", "ACC", "SEN", "SPE", "PPV", "NPV", "Miscall")
		L2 <- c("", "Valid.", "Valid.",
			"(%)", "(%)", "(%)", "(%)", "(%)", "(%)", "(%)")
		L2a <- c("", "Valid.", "Valid.",
			"(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)", "(\\%)")
	}

	d$NumValid <- object$detail$valid.num
	d$FreqValid <- sprintf("%0.4f", object$detail$valid.freq)
	d$FreqValid[d$NumValid <= 0] <- "0"

	d$CR <- sprintf("%0.1f", object$detail$call.rate*100)
	d$CR[object$detail$call.rate < 0.0005] <- "--"

	d$ACC <- sprintf("%0.1f", object$detail$accuracy*100)
	d$ACC[!is.finite(object$detail$accuracy)] <- "--"

	d$SEN <- sprintf("%0.1f", object$detail$sensitivity*100)
	d$SEN[!is.finite(object$detail$sensitivity)] <- "--"

	d$SPE <- sprintf("%0.1f", object$detail$specificity*100)
	d$SPE[!is.finite(object$detail$specificity)] <- "--"

	d$PPV <- sprintf("%0.1f", object$detail$ppv*100)
	d$PPV[!is.finite(object$detail$ppv)] <- "--"

	d$NPV <- sprintf("%0.1f", object$detail$npv*100)
	d$NPV[!is.finite(object$detail$npv)] <- "--"

	d$Miscall <- sprintf("%s (%0.0f)",
		object$detail$miscall, object$detail$miscall.prop*100)
	d$Miscall[is.na(object$detail$miscall) | !is.finite(object$detail$miscall.prop)] <- "--"


	if (type == "txt")
	{
		cat(L1, file=f, sep="\t", append=TRUE)
		cat("\n", file=f, append=TRUE)
		cat(L2, file=f, sep="\t", append=TRUE)
		cat("\n", file=f, append=TRUE)
		cat("----\n", file=f, append=TRUE)
		cat(sprintf("Overall accuracy: %0.1f%%, Call rate: %0.1f%%\n",
			object$overall$acc.haplo*100, object$overall$call.rate*100),
			file=f, append=TRUE)
		write.table(d, file=f, append=TRUE, quote=FALSE, sep="\t",
			row.names=FALSE, col.names=FALSE)

	} else if (type == "tex")
	{
		if (header)
		{
			cat("\\title{Imputation Evaluation}", "",
				"\\documentclass[12pt]{article}", "",
				"\\usepackage{fullpage}",
				"\\usepackage{longtable}", "",
				"\\begin{document}", "",
				"\\maketitle", "",
				"\\setlength{\\LTcapwidth}{6.4in}", "",
				file=f, append=TRUE, sep="\n")
		}

		cat("% -------- BEGIN TABLE --------\n", file=f, append=TRUE)
		if (!is.null(object$detail$train.freq))
		{
			cat("\\begin{longtable}{rrrrr | rrrrrrl}\n", file=f, append=TRUE)
		} else {
			cat("\\begin{longtable}{rrr | rrrrrrl}\n", file=f, append=TRUE)
		}
		cat(c("\\caption{The sensitivity (SEN), specificity (SPE),",
			"positive predictive value (PPV), negative predictive value",
			"(NPV) and call rate (CR)}\n"), file=f, append=TRUE, sep=" ")
		cat("\\label{tab:accuracy} \\\\\n", file=f, append=TRUE)

		cat(L1, file=f, sep=" & ", append=TRUE)
		cat(" \\\\\n", file=f, append=TRUE)
		cat(L2a, file=f, sep=" & ", append=TRUE)
		cat(" \\\\\n", file=f, append=TRUE)
		cat("\\hline\\hline\n\\endfirsthead\n", file=f, append=TRUE)

		if (!is.null(object$detail$train.freq))
		{
			cat("\\multicolumn{12}{c}{{\\normalsize \\tablename\\ \\thetable{} -- Continued from previous page}} \\\\\n",
				file=f, append=TRUE)
		} else {
			cat("\\multicolumn{10}{c}{{\\normalsize \\tablename\\ \\thetable{} -- Continued from previous page}} \\\\\n",
				file=f, append=TRUE)
		}

		cat(L1, file=f, sep=" & ", append=TRUE)
		cat(" \\\\\n", file=f, append=TRUE)
		cat(L2a, file=f, sep=" & ", append=TRUE)
		cat(" \\\\\n", file=f, append=TRUE)
		cat("\\hline\\hline\n\\endhead\n\\hline\n", file=f, append=TRUE)

		if (!is.null(object$detail$train.freq))
		{
			cat("\\multicolumn{12}{r}{Continued on next page ...} \\\\\n",
				file=f, append=TRUE)
		} else {
			cat("\\multicolumn{10}{r}{Continued on next page ...} \\\\\n",
				file=f, append=TRUE)
		}

		cat("\\hline\n\\endfoot\n\\hline\\hline\n\\endlastfoot\n",
			file=f, append=TRUE)

		if (!is.null(object$detail$train.freq))
		{
			cat(sprintf(
				"\\multicolumn{12}{l}{\\it Overall accuracy: %0.1f\\%%, Call rate: %0.1f\\%%} \\\\\n",
				object$overall$acc.haplo*100, object$overall$call.rate*100),
				file=f, append=TRUE)
		} else {
			cat(sprintf(
				"\\multicolumn{10}{l}{\\it Overall accuracy: %0.1f\\%%, Call rate: %0.1f\\%%} \\\\\n",
				object$overall$acc.haplo*100, object$overall$call.rate*100),
				file=f, append=TRUE)
		}

		write.table(d, file=f, append=TRUE, quote=FALSE, sep=" & ",
			row.names=FALSE, col.names=FALSE, eol=" \\\\\n")

		cat("\\end{longtable}\n% -------- END TABLE --------\n",
			file=f, append=TRUE)

		if (header)
		{
			cat("\n\\end{document}\n", file=f, append=TRUE)
		}
	} else if (type == "html")
	{
		if (header)
		{
			cat("<!DOCTYPE html>",
				"<html>",
				"<head>",
				"  <title>Imputation Evaluation</title>",
				"</head>",
				"<body>",
				file=f, append=TRUE, sep="\n")
		}

		cat("<h1>Imputation Evaluation</h1>",
			"<p></p>",
			"<h3><b>Table 1:</b> The sensitivity (SEN), specificity (SPE), positive predictive value (PPV), negative predictive value (NPV) and call rate (CR).</h3>",
			file=f, append=TRUE, sep="\n")

		cat("<table id=\"TB-Acc\" class=\"tabular\" border=\"1\"  CELLSPACING=\"1\">",
			"<tr>", 
				paste(paste("<th>", L1, " ", L2, "</th>", sep=""), collapse=" "),
			"</tr>",
			"<tr>",
			paste("<td colspan=\"", length(L1), "\">", sep=""),
			sprintf("<i> Overall accuracy: %0.1f%%, Call rate: %0.1f%% </i>",
				object$overall$acc.haplo*100, object$overall$call.rate*100),
			"</td>",
			"</tr>",
			file=f, append=TRUE, sep="\n")

		for (i in 1:nrow(d))
		{
			cat("<tr>", 
				paste(paste("<td>", d[i, ], "</td>", sep=""), collapse=" "),
				"</tr>",
				file=f, append=TRUE, sep="\n")
		}

		cat("</table>\n", file=f, append=TRUE)

		if (header)
		{
			cat("\n</body>\n</html>\n", file=f, append=TRUE)
		}
	}

	invisible()
}


##########################################################################
##########################################################################
#
# Linkage Disequilibrium
#

##########################################################################
# To calculate linkage disequilibrium between KIR locus and SNP markers
#

hlaGenoLD <- function(hla, geno)
{
	# check
	stopifnot(inherits(hla, "hlaAlleleClass"))
	if (inherits(geno, "hlaSNPGenoClass"))
	{
		stopifnot(dim(hla$value)[1] == length(geno$sample.id))
		if (any(hla$value$sample.id != geno$sample.id))
		{
			hla <- hlaAlleleSubset(hla, samp.sel =
				match(geno$sample.id, hla$value$sample.id))
		}
		geno <- geno$genotype
	} else if (is.matrix(geno))
	{
		stopifnot(is.numeric(geno))
		stopifnot(dim(hla$value)[1] == dim(geno)[2])
	} else if (is.vector(geno))
	{
		stopifnot(is.numeric(geno))
		stopifnot(dim(hla$value)[1] == length(geno))
		geno <- matrix(geno, ncol=1)
	} else {
		stop("geno should be `hlaSNPGenoClass', a vector or a matrix.")
	}

	# KIR3DL1/S1 alleles indicators
	alleles <- unique(c(hla$value$allele1, hla$value$allele2))
	alleles <- alleles[order(alleles)]
	allele.mat <- matrix(as.integer(0),
		nrow=length(hla$value$allele1), ncol=length(alleles))
	for (i in 1:length(alleles))
	{
		allele.mat[, i] <- (hla$value$allele1==alleles[i]) +
			(hla$value$allele2==alleles[i])
	}

	apply(geno, 1,
		function(x, allele.mat) {
			suppressWarnings(mean(cor(x, allele.mat, use="pairwise.complete.obs")^2, na.rm=TRUE))
		},
		allele.mat=allele.mat)
}




##########################################################################
##########################################################################
#
# Visualization
#

##########################################################################
# To visualize an attribute bagging model
#
#' @export
plot.hlaAttrBagClass <- function(x, ...)
{
	obj <- hlaModelToObj(x)
	plot(obj, ...)
}
#' @export
print.hlaAttrBagClass <- function(x, ...)
{
	obj <- hlaModelToObj(x)
	print(obj)
	invisible()
}

##########################################################################
# To visualize an attribute bagging model
# Here
#' @export
plot.hlaAttrBagObj <- function(x, xlab=NULL, ylab=NULL,
	locus.color="red", locus.lty=2, locus.cex=1.25,
	assembly=c("auto", "hg19", "hg38", "unknown"), ...)
{
	# check
	stopifnot(inherits(x, "hlaAttrBagObj"))
	assembly <- match.arg(assembly)

	# the starting and ending positions of HLA locus
	if (assembly == "auto")
	{
		if (!is.null(x$assembly))
			assembly <- x$assembly
	}
	info <- hlaLociInfo(assembly)
	pos.start <- info$pos.HLA.start[[x$hla.locus]]/1000
	pos.end <- info$pos.HLA.end[[x$hla.locus]]/1000

	# summary of the attribute bagging model
	desp <- summary(x, show=FALSE)

	# x - label, y - label
	if (is.null(xlab)) xlab <- "SNP Position (KB)"
	if (is.null(ylab)) ylab <- "Frequency of Use"

	# draw
	plot(x$snp.position/1000, desp$snp.hist, xlab=xlab, ylab=ylab, ...)
	abline(v=pos.start, col=locus.color, lty=locus.lty)
	abline(v=pos.end, col=locus.color, lty=locus.lty)
	text((pos.start + pos.end)/2, max(desp$snp.hist), paste("HLA", x$hla.locus, sep="-"),
		col=locus.color, cex=locus.cex)
}
#' @export
print.hlaAttrBagObj <- function(x, ...)
{
	summary(x)
	invisible()
}

#######################################################################
# To get the error message
#

hlaErrMsg <- function()
{
	rv <- .C("PONG_ErrMsg", msg=character(1), NAOK=TRUE, PACKAGE="PONG")
	rv$msg
}



#######################################################################
# Internal R library functions
#######################################################################
#' @export
.onAttach <- function(lib, pkg)
{
	# initialize PONG
	rv <- .C("PONG_Init", SSE.Flag=integer(1), PACKAGE="PONG")
	registerS3method("predict", "hlaAttrBagClass", predict.hlaAttrBagClass)
	# information
	packageStartupMessage(
		"PONG (KIR3DL1/S1 Genotype Imputation with Attribute Bagging): v1.0.0")
	if (rv$SSE.Flag != 0)
		packageStartupMessage("Supported by Streaming SIMD Extensions 2 (SSE2)")

	TRUE
}
#' @export
.Last.lib <- function(libpath)
{
	# finalize PONG
	rv <- .C("PONG_Done", PACKAGE="PONG")
	TRUE
}
