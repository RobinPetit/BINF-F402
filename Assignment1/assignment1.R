#!/usr/bin/Rscript

# Lib necessary to use LaTeX in plots
library('latex2exp');
library('gplots');  # for heatmap.2

INFINIUM_PATH  <- 'Infinium450k_raw_data.txt';
ANNO_MINI_PATH <- 'HumanMethylation450_anno_mini.csv';
NB_CASES <- 3
NB_CONTROLS <- NB_CASES
CONTROLS <- 1:NB_CONTROLS
CASES <- (1:NB_CASES) + NB_CONTROLS;
ALPHA_PSEUDO_COUNT <- 100  # see: Du et al., Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis, 2010
ALL_CHROMOSOMES <- c(1:22, 'X', 'Y');
CHROMOSOMES_TO_TEST <- 19:21;

LOG_P_VALUE_THRESHOLD <- -3;
DELTA_BETA_THRESHOLD <- .2;

# For Infinium array, `signal A` is the unmethylated probe and
# `signal B` is the methylated probe

compute.beta <- function(signalA, signalB) {
	# Note the pseudo count $\alpha = 100$ to avoid dividing by 0
	return (signalB / (signalA + signalB + ALPHA_PSEUDO_COUNT));
}

# Get the column corresponding to the required sample
get.sample.signal <- function(table, sample.id, signalA=T) {
	if(signalA) {
		return (as.vector(table[[2*sample.id-1]]));
	} else {
		return (as.vector(table[[2*sample.id]]));
	}
}

# Returns the $\beta$ values of a given sample
get.beta <- function(infinium, sample.id) {
	return (compute.beta(
		get.sample.signal(infinium, sample.id, T),
		get.sample.signal(infinium, sample.id, F)
	));
}

# Plot the beta distribution of each sample from 1 to 6
plot.beta.distribution.infinium <- function(infinium) {
	for(sample.id in 1:(NB_CASES+NB_CONTROLS)) {
		plot(
			density(get.beta(infinium, sample.id)),
			xlim=c(-.2, 1.2),
			ylim=c(0, 5),
			xlab=TeX('$\\beta$'),
			ylab=TeX('Density of $\\beta$'),
			main=TeX(paste('density of $\\beta$ distribution of sample', sample.id)),
		);
	}
}

# Returns the subtable corresponding to the required chromosome
extract.chromosome <- function(annotations, chr) {
	return (annotations[annotations$CHR == chr,]);
}

# Make the barplots of the average methylation level of each chromosome
plot.mean.methylation.level <- function(infinium, annotations) {
	means <- matrix(rep(0, 2*length(ALL_CHROMOSOMES)), nrow=2, ncol=length(ALL_CHROMOSOMES));
	idx <- 1;
	for(chromosome in ALL_CHROMOSOMES) {
		intersection <- intersect(rownames(infinium), rownames(extract.chromosome(annotations, chromosome)));
		chromosome.infinium.subtable <- infinium[intersection,];
		for(case.sample.id in CASES)
			means[1,idx] <- means[1,idx] + mean(get.beta(chromosome.infinium.subtable, case.sample.id))/NB_CASES;
		for(control.sample.id in CONTROLS)
			means[2,idx] <- means[2,idx] + mean(get.beta(chromosome.infinium.subtable, control.sample.id))/NB_CONTROLS;
		idx <- idx+1;
	}
	barplot(means, names.arg=ALL_CHROMOSOMES, horiz=T, beside=T, legend=c('Cases', 'Controls'), las=1);
	barplot(means[,CHROMOSOMES_TO_TEST], names.arg=CHROMOSOMES_TO_TEST, beside=T, legend=c('Cases', 'Controls'), las=1);
}

# Plot the histograms of the distribution and log-distribution of the number of probes per gene
plot.number.of.genes.per.probe <- function(annotations) {
	non.null.genes.associations <- as.vector(annotations[annotations$UCSC_RefGene_Name != '',]$UCSC_RefGene_Name);
	genes.counter <- table(unlist(sapply(non.null.genes.associations, function(s) {unique(strsplit(s, ';'))})));
	hist(genes.counter, main='Distribution of the number of probes per gene',
		sub=sprintf('Average number of probes per gene is %g', mean(genes.counter)),
		xlab='Number of probes per gene',
		ylab='Frequency');
	plot(log(table(genes.counter)), main='log plot of the distribution of the number of probes per gene',
		xlab='Number of probes per gene',
		ylab='log(Frequency)');
}

# Get the $\Delta\beta$ vector of a (control, sample) couple
get.delta.beta <- function(infinium, case.sample.id, control.sample.id) {
	return (get.beta(infinium, control.sample.id) - get.beta(infinium, case.sample.id));
}

# Plots the density of $\Delta\beta$ for each (control, sample) couple
plot.delta.beta.distribution <- function(infinium) {
	# plot \Delta \beta distribution for each couple control/case
	for(control.sample.id in CONTROLS) {
		case.sample.id <- control.sample.id + NB_CONTROLS;
		infinium$de
		plot(
			density(get.delta.beta(infinium, case.sample.id, control.sample.id)),
			main=TeX(sprintf('$\\Delta\\beta$ (control - case) for couple #%d', control.sample.id)),
			xlab=TeX('$\\Delta\\beta$'),
			ylab='Density',
			xlim=c(-1, 1),
			ylim=c(0, 5.5)
		);
	}
	# Plot \Delta \beta for the mean of all samples
	delta.betas <- matrix(0, nrow=nrow(infinium), ncol=NB_CONTROLS+NB_CASES);
	for(control.sample.id in CONTROLS) {
		delta.betas[,control.sample.id] <- get.beta(infinium, control.sample.id);
		case.sample.id <- control.sample.id + NB_CONTROLS;
		delta.betas[,case.sample.id] <- get.beta(infinium, case.sample.id);
	}
	infinium$delta.beta <- apply(delta.betas[,CONTROLS], 1, mean) - apply(delta.betas[,CASES], 1, mean);
	plot(
		density(infinium$delta.beta),
		main=TeX('mean $\\Delta\\beta$ (controls - cases)'),
		xlab=TeX('$\\Delta\\beta$'),
		ylab='Density',
		xlim=c(-1, 1),
		ylim=c(0, 5.5)
	);
	infinium$p.values <- apply(delta.betas, 1,
		function(row) {
			return (t.test(row[CONTROLS], row[CASES])$p.value);
		}
	);
	infinium$log.p.value <- log(infinium$p.value, 10);
	plot(
		density(infinium$p.values),
		main=TeX('$p$-value density'),
		xlab=TeX('$p$-value'),
		ylab='density'
	);
	plot(
		density(-infinium$log.p.value),
		main=TeX('$-\\log_{10}(p$-value$)$ density'),
		xlab=TeX('$-\\log_{10}(p$-value$)$'),
		ylab='density'
	);
	return (infinium);
}

plot.volcano <- function(infinium) {
	significant.idx <- as.logical((abs(infinium$delta.beta) >= .2) * (infinium$log.p.value <= -3));
	significant <- infinium[significant.idx,];
	non.significant <- infinium[!significant.idx,];
	plot(
		c(non.significant$delta.beta, significant$delta.beta),
		-c(non.significant$log.p.value, significant$log.p.value),
		col=c(rep('black', nrow(non.significant)), rep('red', nrow(significant))),
		xlab=TeX('$\\Delta\\beta$'),
		ylab=TeX('$-\\log_{10}(p$-value$)$')
	)
	lines(c(min(infinium$delta.beta), max(infinium$delta.beta)), rep(-LOG_P_VALUE_THRESHOLD, 2), col='grey', lty=2)
	lines(rep(DELTA_BETA_THRESHOLD, 2), c(min(-infinium$log.p.value), max(-infinium$log.p.value)), col='grey', lty=2)
	lines(rep(-DELTA_BETA_THRESHOLD, 2), c(min(-infinium$log.p.value), max(-infinium$log.p.value)), col='grey', lty=2)
}

get.top.delta.beta <- function(infinium, n=30) {
	return (infinium[order(abs(infinium$delta.beta), decreasing=T)[1:n],]);
}

plot.heatmap <- function(infinium) {
	top <- get.top.delta.beta(infinium);
	par(mar=c(0, 1, 3, 1))
	heatmap.2(t(as.matrix(top[,1:(2*(NB_CONTROLS+NB_CASES))])), Rowv=NA, Colv=NA, scale='none', cexRow=.8, cexCol=.75,
		margins=c(6, 8), dendogram='none', trace='none');
	for(i in 1:2) {
		heatmap.2(t(as.matrix(top[,seq(i, 2*(NB_CONTROLS+NB_CASES),2)])), Rowv=NA, Colv=NA, scale='none', cexRow=.8, cexCol=.75,
			margins=c(6, 8), dendogram='none', trace='none');
	}
}

# Use first column as row names
infinium <- read.table(INFINIUM_PATH, header=T, dec=',', row.names=1);
annotations <- read.table(ANNO_MINI_PATH, header=T, sep=',', row.names=1);

# Plots
plot.beta.distribution.infinium(infinium);
plot.mean.methylation.level(infinium, annotations);
plot.number.of.genes.per.probe(annotations);
infinium <- plot.delta.beta.distribution(infinium);
plot.volcano(infinium);
plot.heatmap(infinium);
