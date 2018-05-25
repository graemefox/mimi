**.----------------.  .----------------.  .----------------.  .----------------. **
**| .--------------. || .--------------. || .--------------. || .--------------. |**
**| | ____    ____ | || |     _____    | || | ____    ____ | || |     _____    | |**
**| ||_   \  /   _|| || |    |_   _|   | || ||_   \  /   _|| || |    |_   _|   | |**
**| |  |   \/   |  | || |      | |     | || |  |   \/   |  | || |      | |     | |**
**| |  | |\  /| |  | || |      | |     | || |  | |\  /| |  | || |      | |     | |**
**| | _| |_\/_| |_ | || |     _| |_    | || | _| |_\/_| |_ | || |     _| |_    | |**
**| ||_____||_____|| || |    |_____|   | || ||_____||_____|| || |    |_____|   | |**
**| |              | || |              | || |              | || |              | |**
**| '--------------' || '--------------' || '--------------' || '--------------' |**
**'----------------'  '----------------'  '----------------'  '----------------' **

# 'multi individual microsatelite identification' - (mimi)
### A tool to improve the design of novel microsatellite panels from genome next-generation sequencing data.


#### What does it do?
mimi is a Python script that attempts to build on the microsatellite markers design process [pal_finder](https://sourceforge.net/projects/palfinder/)
by increasing the rate at which markers amplify by PCR and also allows the user to select polymorphic loci from the data. It does this by
using the genome data from several individuals of the same species, rather than from a single indivdual which is more common in the microsatellite
design process.


#### What does it allow me to do?
mimi allows you to visualise three important pieces of information which are not available when designing microsatellite markers from the genome of a single individual.

1) You can select primer pairs which show strong sequence conservation across several individuals. This gives a much higher rate of PCR success and should allow for a
reduction in the frequency of null alleles (in theory...). (fig. 1).

**Fig 1**

A potential primer highlighted in red, showing 100% sequence conservaion in four individuals above.

![Figure1 - strong sequence conservation](/images/fig1.png)


2) By looking at the microsatellite locus itself in several individuals you can select loci which are polymiorphic in the amount of microsatellite repeats
and avoid loci where all individuals appear to have the same number of repeats. (fig. 2).

**Fig. 2**

A "CATA" microsatellite motif in five individuals highlighted in red. The number of repeats varies between individuals with strong sequence conservation either side of the microsatellite.

![Figure2 - variable number of repeats](/images/fig3.png)

3) You can detect whether a potential microsatellite marker contains other fragment length altering polymorphisms, outside of those caused by the change in number
of motif repeats. Insertion/deletions in the flanking regions are clearly visible and can be avoided when designing your microsatellite panel. Mutations of thus sort
are an important source of error and would otherwise be very difficult to detect in a panel designed in a single individual. (fig. 3).

**Fig. 3**

A section of flanking sequence (in between a priming region and the microsatellite itsel). One individual has a significant deletion mutation (highlighted in red) in otherwise strong
sequence conservation.

![Figure3 - insertion/deletion mutation](/images/fig2.png)

### Versions
v0.1 - This is the initial version which should establish core functionality. Bugs are likely to be found and improvements/modifications will come in the future.
This was made by a humble biology PhD student, not a "programmer", so please be patient!

### Installation Instructions
This was made and tested on [Ubuntu Linux](https://www.ubuntu.com/) (currently 18.04) but *should* also work on OSX (Your mileage may vary).

#### Dependencies
You **must** have [Biopython](https://biopython.org/) and [PANDAseq](https://biopython.org/) installed.
PANDAseq **must** be in your $PATH.

ie. Run the following from a terminal on your system.
```
pandaseq
```

If PANDAseq is correctly installed, you will see:

```
You must supply both forward and reverse reads.
Too confused to continue.
Try -h for help.
```

#### Set up the configuration file
There is a mimi_config.txt file which contains the parameters used to control mimi. Open it in a **plain text** editor (not a word processor) and change the following fields to accomodate your data:

```
Amount of individuals sequenced:
number_of_samples = 3
```

MiMi will run on anything >1 sample but you will not get particularly meaningful results with a small number. I recommend eight samples (eight individuals).

```
Proportion of individuals in which a microsatellite loci should occur (default 0.5)
proportion_of_individuals = 0.5
```

"proportion_of_individuals" controls in how many of your sequence files a microsatellite locus must be discovered to pass through the filters. The default value is 0.5
but if you find you are not getting many results, this can be reduced.

```
Each Individual should have two input files: paths to "R1" and "R2" file in the
following format:


Input your real data below:
input1_R1 = /path/to/individual_1_R1.fastq
input1_R2 = /path/to/individual_1_R2.fastq
input2_R1 = /path/to/individual_2_R1.fastq
input2_R2 = /path/to/individual_2_R2.fastq
```

These are the paths to your raw FASTQ files. Each individual should have two entries, one for the forward of the pair of reads (R1) and one for the reverse of the pair of reads (R2).
Extend the list if neccessary by adding new rows. Single-end sequencing reads are not supported.

```
Prior to running mimi, each individual should have been processed with pal_finder and the
additional filtering and PANDAseq QC steps by Griffiths et al (unpublished)
with the path to the output file given below in the following format:

input1_pal_finder = /path/to/individual_1_pal_finder_output.txt
input2_pal_finder = /path/to/individual_2_pal_finder_output.txt
```

These are the paths to the output files from the Griffiths et al (2016) workflow. This must have been performed on your data prior to running MiMi as the process depends on
these output files for the loci and primer sequences.

```
Path to pal_finder_v0.02.04 script and config file:

pal_finder_path = /path/to/pal_finder_v0.02.04.pl
pal_finder_config = /path/to/pal_finder/config.txt
```

Finally, your pal_finder installation shipped two files: "pal_finder_v0.02.04.pl" and "config.txt". Provide the paths for these two files. Do not alter the config file as this is accessed and modified
by the mimi script. Please note there are two config files; one which is accessed by MiMi (mimi_config.txt by default) and one which is accessed by pal_finder (config.txt by default).
 These are different and both are required.

#### Run the script
You may need to give mimi permission to run:
```
sudo chmod +x ./mimi.py
```



#### Testing with the demo data
When you clone the mimi repo, it supplies some demo data to test your installation.

Run the script and pass the demo_config file:
```
./mimi.py -c /demo_data/demo_config.txt
```

Run the script and pass the MiMi config file:
```
./mimi.py -c /path/to/config.txt
```
On my modest desktop machine, using data from a single Miseq run (approx 12-16Gb) the MiMi process runs in approximately four hours.

#### Interpret the results




#### Who made this?
[Graeme Fox](https://graemefox.github.io)
