# 'Multi-individual Microsatelite identification' - (MiMi)
### A tool to improve the design of novel microsatellite panels from genome next-generation sequencing data.

#### What does it do?
MiMi is a script that attempts to build on the microsatellite markers design process [pal_finder](https://sourceforge.net/projects/palfinder/)
by increasing the rate at which markers amplify by PCR and also allows the user to select polymorphic loci from the data. It does this by
using the genome data from several individuals of the same species, rather than from a single indivdual which is more common in the microsatellite
design process.


#### What does it allow me to do?
MiMi allows you to visualise three important pieces of information which are not available when designing microsatellite markers from the genome of a single individual.

Firstly, you can select primer pairs which show strong sequence conservation across several individuals. This gives a much higher rate of PCR success and should allow for a 
reduction in the frequency of null alleles (in theory...).

![Figure1 - strong sequence conservation](/images/fig1.png)


Secondly, by looking at the microsatellite locus itself in several individuals you can select loci which are polymiorphic in the amount of microsatellite repeats
and avoid loci where all individuals appear to have the same number of repeats.

![Figure2 - variable number of repeats](/images/fig3.png)

Finally, you can detect whether a potential microsatellite marker contains other fragment length altering polymorphisms, outside of those caused by the change in number
of motif repeats. Insertion/deletions in the flanking regions are clearly visible and can be avoided when designing your microsatellite panel. Mutations of thus sort
are an important source of error and would otherwise be very difficult to detect in a panel designed in a single individual.

![Figure3 - insertion/deletion mutation](/images/fig2.png)

#### Installation Instructions
This was made and tested on [Ubuntu Linux](https://www.ubuntu.com/) (currently 18.04) but *should* also work on OSX (YMMV).

### Dependencies
You **must** have [Biopython](https://biopython.org/) and [PANDAseq](https://biopython.org/) installed.

#### Who made this?
[Graeme Fox](graemefox.github.io)




