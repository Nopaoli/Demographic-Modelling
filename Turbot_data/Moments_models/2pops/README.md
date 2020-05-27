# Two-population models for *moments*

Here are the 32 two-population models tested in paper. 

They represent a fully orthogonal combination of the following scenarios: 

**Gene Flow**= Isolation with Migration (IM) vs Secondary Contact (SC) mdodels, both with contemporary asymmetric migration

**Genomic islands models** = models with (2M models)  or without heterogeneous migration rates across the genome 
**Backgroundselection/sweeps models** = models with (2N models)  or without heterogeneous Ne  across the genome 
**Ancestral expansion models**= models with (AE models) or without Ne change in the ancestral population before the split
**"Founder" models** = models with (B models) or without a bottleneck followed by growth at the time of split in the Baltic Sea population

In total, this results in 16 IM models and 16 SC models. 

The models assume the data is provided in the form of a &delta;&alpha;&delta;&iota; data dictionaryformat, from which the folded jAFS is computed. The data dictionary containing 1 SNP per 2b-RAD locus is given in the "Turbot_data" folder.  

For running multiple optimizations, i used a bunch of bash scripts that replicate what is done by the "Optimize_Functions.py" from Daniel Portik's pipeline and parse the data. The simplest (but not the fastes) way to run the optimization i pefromed in the manuscript is to copy-paste the model definition and parameter bounds into the optimization scripts i used for the simulated data. 

