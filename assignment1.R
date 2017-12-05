#!/usr/bin/Rscript

# Lib necessary to use LaTeX in plots
library('latex2exp');

INFINIUM_PATH  <- 'Infinium450k_raw_data.txt';
ANNO_MINI_PATH <- 'HumanMethylation450_anno_mini.csv';
NB_CASES <- 3
NB_CONTROLS <- 3
CONTROLS <- 1:NB_CONTROLS
CASES <- (1:NB_CASES) + NB_CONTROLS;
ALL_CHROMOSOMES <- c(1:22, 'X', 'Y');
CHROMOSOMES_TO_TEST <- 19:21;

# For Infinium array, `signal A` is the unmethylated probe and
# `signal B` is the methylated probe

compute.beta <- function(signalA, signalB) {
	# Note the pseudo count $\\alpha = 1$ to avoid dividing by 0
	return (signalB / (signalA + signalB + 1));
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
		for(case.sampel.id in CASES)
			means[1,idx] <- means[1,idx] + mean(get.beta(chromosome.infinium.subtable, case.sampel.id))/NB_CASES;
		for(control.sample.id in CONTROLS)
			means[2,idx] <- means[2,idx] + mean(get.beta(chromosome.infinium.subtable, control.sample.id))/NB_CONTROLS;
		idx <- idx+1;
	}
	barplot(means, names.arg=ALL_CHROMOSOMES, horiz=T, beside=T, legend=c('Cases', 'Controls'), las=1);
	barplot(means[,CHROMOSOMES_TO_TEST], names.arg=CHROMOSOMES_TO_TEST, beside=T, legend=c('Cases', 'Controls'), las=1);
}

# ...
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

# Use first column as row names
infinium <- read.table(INFINIUM_PATH, header=T, dec=',', row.names=1);
annotations <- read.table(ANNO_MINI_PATH, header=T, sep=',', row.names=1);
plot.beta.distribution.infinium(infinium);
plot.mean.methylation.level(infinium, annotations);
plot.number.of.genes.per.probe(annotations);
