# Analyses of simulated data

Here you can find all necessary commands and models to replicate the analyses of simulated data.

## MS simulations

Within this subfolder you have several folders where the codes to recreate the 64 combinations of ancestral size expansion and bottlegrowth as well as IM and SC scanerios with symmetric and asymmetric migration.

Within the "100k_Loci folder" you have the codes to to recreate the simulations for 100 k 36 bp loci , assuming a mutation rate of 1x10<sup>-8</sup> and no recombination within locus. This is the equivalent of sampling roughly 100 k 2b-RAD loci. This is achived by running 100k simulations with a &theta; of 0.0288, where \theta is based on the reference population size N<sub>0</sub> which we set to be the same as N<sub>1</sub> at 20 000 individuals.  So &theta; is defined as 4N<sub>0</sub>&mu;. The parameter &mu; is the locus mutation rate, in this case since the locus lenght is 36 bp &mu; is 36x10<sup>e-8</sup> and  &theta; is 4x20000x6x10<sup>e-8</sup> which is 0.0288.

Within the "1_Million_Loci folder" you have the codes to to recreate the simulations for 1 million 36 bp loci. In this case &theta; is ten times the value above, 0.288

When analysing the data from MS simulations in *dadi* or *moments* simply read use the commands
dadi.Spectrum.from_ms_file (infile, average=False)
moments.Spectrum.from_ms_file (infile, average=False)

## *dadi*  models

Within the Dadi_models folder you have the the basic IM and SC models used to analyse the simulations with  *dadi*, as well as functions for optmizing and parsing the results. These are modified versions of Daniel Portik's*dadi* pipeline (https://github.com/dportik/dadi_pipeline)

## *moments* models

Within the Moments_models folder you have the 8 models (IM, SC; IM_B, SC_B, IM_AE, SC_AE, IM_AE_B, SC_AE_B) used for analysing the simulated data with *moments* as well as functions to optimize and summarize the results. The optimization functions are a modfied verions from Daniel Portik's dadi pipeline that works in *moments*.
