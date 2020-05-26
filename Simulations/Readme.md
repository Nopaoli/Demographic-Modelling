# Analyses of simulated data

Here you can find all necessary data to replicate the analyses of simulated data from the paper Momigliano et al. 2020.

## MS Simulations

Within this subfolder you have several folders where the codes to recreate the 64 combinations of ancestral size expansion and bottlegrowth as well as IM and SC scanerios with symmetric and asymmetric migration

Within the "100k_Loci folder" you have the codes to to recreate the simulations for 100 k 36 bp loci , assuming a mutation rate of 1x10<sup>-8</sup> and no recombination within locus. This is the equivalent of sampling roughly 100 k 2b-RAD loci. This is achived by running 100k simulations with a &theta; of 0.0288, where \theta is based on the reference population size N<sub>0</sub> which we set to be the same as N<sub>1</sub> at 20 000 individuals.  So &theta; is defined as 4N<sub>0</sub>&mu;. The parameter &mu; is the locus mutation rate, in this case since the locus lenght is 36 bp &mu; is 36x10<sup>e-8</sup> and  &theta; is 4x20000x6x10<sup>e-8</sup> which is 0.0288.

Within the "100k_Loci folder" you have the codes to to recreate the simulations for 1 million 36 bp loci. In this case &theta; is ten times the value above, 0.288

When analyses the data from MS simulations in &delta;&alpha;&delta;&iota; or *moments* simply read use the commands
dadi.Spectrum.from_ms_file (infile, average=False)
moments.Spectrum.from_ms_file (infile, average=False)

## Dadi Models

Within the Dadi_models folder you have the the basic IM and SC models used to analyse the simulations with  &delta;&alpha;&delta;&iota;, as well as functions for optmizing and parsing the results. These are modified versions of Daniel Portik's &delta;&alpha;&delta;&iota; pipeline. 

## Moments models

Within the Moments_models folder you have the 8 models (IM, SC; IM_B, SC_B, IM_AE, SC_AE, IM_AE_B, SC_AE_B) used for analysing teh simulated data with *moments* as well as functions to optimize and summarize the results. The optimization functions are a modfied verions from Daniel Portik's Daniel Portik's &delta;&alpha;&delta;&iota; pipeline that works in *moments*.
