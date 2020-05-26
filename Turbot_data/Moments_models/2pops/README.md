# Two-population models for *Moments*

Here are the 32 two-population models tested in paper. 

They represent a fully orthoginal combination of teh following scenarios: 

**Gene Flow**= Isolation with Migration (IM) vs Secondary Contact (SC) mdodels, both with contemporary asymmetric migration

**Genomic islands models** = models with (2M models)  or without heterogeneous migration rates across the genome 
**Backgroundselection/sweeps models** = models with (2N models)  or without heterogeneous Ne  across the genome 
**Ancestral expansion models**= models with (AE models) or without Ne change in the ancestral population before the split
**"Founder" models** = models with (B models) or without a bottleneck followed by growth at the time of split in the Baltic Sea population

In total, this results in 16 IM models and 16 SC models. 

The models use a data dictiory in &delta;&alpha;&delta;&iota; format, from which the folded jAFS is computed

For running multiple optimizations, i used a nunch of bash scripts that replicate what is done by the "Optimize_Functions.py" from Daniel Portik's pipeline and parse the data.

