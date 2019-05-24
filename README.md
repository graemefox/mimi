# 'Multi individual Microsatelite identification' - (MiMi)
### A tool to improve the design of novel microsatellite panels from genomic next-generation sequencing data.


#### What does it do?
MiMi is a Python script that attempts to build on the microsatellite markers design process [pal_finder](https://sourceforge.net/projects/palfinder/) (Castoe et al., 2012)
by increasing the rate at which markers amplify by PCR, allows the user to select polymorphic loci from the data and allows the avoidance of insertion/deletions in the flanking regions. It does this by
using genomic data from several individuals of the same species, rather than from a single indivdual which is more common in the microsatellite
design process.


#### What does it allow me to do?
MiMi allows you to visualise three important pieces of information which are not available when designing microsatellite markers from the genome of a single individual.

1) You can select primer pairs which show strong sequence conservation across several individuals. This gives a much higher rate of PCR success and should allow for a
reduction in the frequency of null alleles (in theory...). (Fig. 1).

**Fig 1**

A potential primer highlighted in red, showing 100% sequence conservaion in four individuals above.

![Figure1 - strong sequence conservation](/images/fig1.png)

2) By looking at the microsatellite locus itself in several individuals you can select loci which appear to be polymorphic in the amount of microsatellite repeats
and avoid loci where all individuals appear to have the same number of repeats. (Fig. 2).

**Fig. 2**

A "CATA" microsatellite motif in five individuals highlighted in red. The number of repeats varies between individuals with strong sequence conservation either side of the microsatellite.

![Figure2 - variable number of repeats](/images/fig3.png)

3) You can detect whether a potential microsatellite marker contains other fragment length altering polymorphisms, outside of those caused by the change in number
of motif repeats. Insertion/deletions in the flanking regions are clearly visible and can be avoided when designing your microsatellite panel. Mutations of this sort
are an important source of error and would otherwise be very difficult to detect in a panel designed in a single individual. (Fig. 3).

**Fig. 3**

A section of flanking sequence (in between a priming region and the microsatellite itsel). One individual has a significant deletion mutation (highlighted in red) in otherwise strong
sequence conservation.

![Figure3 - insertion/deletion mutation](/images/fig2.png)

### Versions
v0.2 (release date: ?????) - This is an updated version of the MiMi software. The core functionality of locus detection is unchanged. I have implemented some additional, automated filtering
of the output. Alignments are automatically aligned using MUSCLE (additional software requirements now apply), alignments are trimmed to the primer positions, low quality alignments are
automatically removed (based upon a metric of 'gappiness' of an alignment - it's crude but effective), loci are ranked by the conservation of primer sequences. A new log file is produced
detailing the proportions of putative markers which have been removed by the filters, and also retaining all the loci should the user want to access any markers which have been filtered out.

v0.1 (no longer available) - This is the initial version which should establish core functionality. Bugs are likely to be found and improvements/modifications will come in the future.
This was made by a humble biology PhD student, not a "programmer", so please be patient!

### Installation Instructions
This was made and tested on [Ubuntu Linux](https://www.ubuntu.com/) (currently 18.04) but *should* also work on OSX (your mileage may vary).

#### Dependencies
You need to have [Python](https://www.python.org/), [Biopython](https://biopython.org/), [MUSCLE] (https://www.drive5.com/muscle/downloads.htm) and [PANDAseq](https://github.com/neufeld/pandaseq) installed.
PANDAseq and MUSCLE must be in your $PATH.

#### Test PANDAseq installation:
Run the following from a terminal on your system.
```
pandaseq
```

If PANDAseq is correctly installed, you will see:

```
You must supply both forward and reverse reads.
Too confused to continue.
Try -h for help.
```

#### Test Biopython Installation
Run the following from a terminal on your system:
```
python
```

and then at the Python prompt, type:
```
import Bio
```

If no error message appears, Biopython is correctly installed.

Type:
```
exit()
```

to exit the Python prompt.

#### Test MUSCLE installation
Run the following from a terminal on your system:
```
muscle
```

You should see a lot of information relating to the muscle program, which begins:
(although the version number may vary)
```
MUSCLE v3.8.31 by Robert C. Edgar
```

#### Troubleshooting error messages
Upon running MiMi, if you receive the error message: "ImportError: No module named Bio" this indicates that Biopython is not correctly installed.
Please review the installation instructions at "https://biopython.org/DIST/docs/install/Installation.html" to ensure the module is correctly installed.

If you do not see the information relating to the MUSCLE program, this is indicative that the program is not properly installed.
Please refer to the installation instructions at (https://www.drive5.com/muscle/manual/install.html) to ensure the program is correctly installed.

#### Clone the MiMi repository to your disk
```
git clone https://github.com/graemefox/mimi.git
```

#### Set up the configuration file
There is a MiMi_config.txt file which contains the parameters used to control MiMi. Open it in a **plain text** editor (not a word processor) and change the following fields to accomodate your data:

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
# Prior to running MiMi, each individual should have been processed with pal_finder and the
# additional filtering and PANDAseq QC steps by Griffiths et al (2016)
# This is most easily achieved using the Galaxy hosted tool at the University of Manchester,
# here: https://palfinder.ls.manchester.ac.uk/
# with the path to the "filtered_microsatellites_(full_details)].tabular" output file given below in the following format:

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

Finally, MiMi also brings the pal_finder tool when you download which contains two files: "pal_finder_v0.02.04.pl" and "config.txt". Provide the paths for these two files. Do not alter the config file as this is accessed and modified
by the MiMi script. Please note there are two config files; one which is accessed by MiMi (MiMi_config.txt by default) and one which is accessed by pal_finder (config.txt by default).
These are different and both are required.

#### Run the script
You may need to give MiMi permission to run:
```
sudo chmod +x MiMi.py

```

#### Usage
```
./MiMi.py -c /path/to/config.txt
```

The usage is very simple. You only need to pass the MiMi_config file with the "-c" flag. All other settings are contained within the config file itself.

#### Testing with the demo data
When you clone the MiMi repository, it supplies some demo data and a pre-configured MiMi_config file to test your installation.

Run the script and pass the demo_config file:
```
./MiMi.py -c /demo_data/demo_config.txt
```

#### Interpret demo data output
If the script ran correctly you will see a directory named "MiMi_output" in the MiMi directory. Within this directory is a file "MiMi_output.txt" and a directory "Alignments".
The demo data consisted of small 'shotgun', paired-end sequencing datasets of four individuals provided in the "sequence_data" directory of the demo data.
Microsatellites had already been detected in these individuals and primers designed
using the Griffiths et al. (2016) workflow (available here: https://palfinder.ls.manchester.ac.uk/). The output files from the Griffiths workflow were provided in the "pal_filter_output"
directory in the demo data. MiMi detected that one of these primer pairs was found in three of the four individuals
and has extracted those reads and placed them into a FASTA file in the Alignments directory. The FASTA file is named with the forward primer sequence. Furthermore, in the "MiMi_output.txt"
file, the three alleles which were found at this locus are listed (Fig 4).

**Fig. 4**

![Figure4 - Demo output](/images/demo_output.png)
Figure 4. Showing the MiMi output for one microsatellite locus. The primer sequences are provided along with the number of alleles which have been detected, the number of individuals' datasets
in which this locus has been detected, the alleles present (numbers in brackets represent the number of repeats) and the size range between the smallest and largest allele.

### Run MiMi with real data

#### Next-generation sequencing of your samples.
MiMi requires paired-end, genomic sequence data in FASTQ format. In our lab we generate sequence data using the Illumina Nextera protocol and sequence using a MiSeq platform. Generally eight individuals are sequenced on a single MiSeq flowcell. Ideally, these would be from multiple sites to help MiMi counter any site-specific variation in the primer regions which may lead to null alleles.

#### Required bioinformatics prior to MiMi
You should have already detected microsatellite loci and designed primers using the workflow described in Griffiths et al. (2016.) This is most easily performed using the Galaxy version of the tool hosted at the University of Manchester (https://palfinder.ls.manchester.ac.uk/). This must be performed seperately for each of your individual samples. The Galaxy workflow will produce several output files for each dataset; it is the files containing "filtered_microsatellites_(full_details)].tabular" in the filename which are required by MiMi and are referred to as the "pal_filter output files" in the MiMi documentation.


Configure the MiMi_config.txt file to contain paths to each of your paired-end sequencing FASTQ files, each of your pal_filter output files and the pal_finder scripts. Finally, run the script and pass the MiMi_config file:
```
./MiMi.py -c /path/to/MiMi_config.txt
```
On my modest desktop machine, using data from a single Miseq run (approx 12-16Gb) the MiMi process runs in approximately four hours.

#### Interpret the results
Data interpretation is identical to that described in the demo_data section, however you will hopefully have many more results. The rows in the "MiMi_output.txt" are ranked by the
"Size Range" column as we propose that a large range in allele size is most likely to be indicative of a true polymorphic microsatellite as opposed to a sequencing error which may
result in smaller slippages producing an inflated number of alleles.

#### Cite MiMi
If you found this tool useful please cite MiMi:

Fox, G., Antwis, R., Preziosi, R.F. and Rowntree, J.K. (in progress) Multi individual Microsatellite identification (MiMi). A microsatellite design workflow incorporating multiple genomes.

(MiMi is currently a work in progress but will be submitted for peer review and publication very soon).

#### References
Castoe, T.A., Poole, A.W., Jason de Koning, A. P., Jones, K.L., Tomback, D.F., Oyler-McCance, S.J., Fike, J.A., Lance, S.L., Streicher, J.W., Smith, E.N. and Pollock, D.D. (2012) Rapid Microsatellite Identification from Illumina Paired-End Genomic Sequencing in Two Birds and a Snake. *PLoS ONE*. 7(2): e30953

Griffiths, S.M., Fox, G., Briggs, J., Donaldson, I.J., Hood, S., Richardson, P., Leaver, G.W., Truelove, N.K. and Preziosi, R.F. (2016) A Galaxy-based bioinformatics pipeline for optimised, streamlined microsatellite development from Illumina next-generation sequencing data. *Conservation Genetics Resources.* 8(4)pp. 481-486.

#### Who made this?
[Graeme Fox](https://graemefox.github.io)
