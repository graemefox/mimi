#!/usr/bin/python3 -tt
import configparser, subprocess, os, time, csv, Bio, argparse, re, shutil, time, tempfile, argparse, numpy, curses, sys
from Bio import SeqIO
from subprocess import Popen, PIPE
from os import walk
## update Aug 2021
## Updated script to run on Python3.x rather than require Python2.x due to
# Python2 reaching end of life.
#
## updates: 23/05/2018
## added additional filtering functions to the output files
## MiMi now automatically detects and filters/retains:
## a) polymorphic/monomorphic loci
## b) low quality alignments (developed a new metric to test the overall
#'quality' of an alignment)
## c) ranks the quality of the conservation of primer regions based
#on the following:

# - all forward primer regions should be perfectly conserved otherwise would #
# not pass MiMi
# - full length, perfectly conserved regions in >1 individuals are given h
#ighest priority
# - partial, perfectly conserved regions in > 1 individuals are given
#second priority
# - no additional primer info available are given third priority
# - primer regions where we can visualise mutations (SNP or INDEL etc.)
#are filtered out

## a new log file is procued giving statistics on how many putative markers
#have been removed,
## and under which filtering conditions. All putative markers are retained
#n the log
## file should the user wish to access them

###########################################################
# FUNCTION LIST
###########################################################
# Reverse complement a sequence
def ReverseComplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    return "".join([seq_dict[base] for base in reversed(seq)])

# write a config file suitable for pal_finder
def write_to_config(filename, old_string, new_string):
    s=open(filename).read()
    if old_string in s:
        s=s.replace(old_string, new_string)
        f=open(filename, 'w')
        f.write(s)
        f.flush()
        f.close()

# take a file, make a copy and name it something else
def copy_and_rename_file(old_file_name, new_file_name):
        in_file = open(old_file_name)
        out_file = open(new_file_name, "w")
        for line in in_file:
            out_file.write(line)
        in_file.close()
        out_file.close()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# get sequences containing a primer sequence from a fastq file

#input is a list of primers
#it searches FastQ files for sequences which contain those primer sequences
#and outputs ID and sequence information

def get_seqs(script, get_paired, sequencefile, paired_sequence_file, primer_seq):
    with tempfile.NamedTemporaryFile(mode='w') as scriptfile:
        scriptfile.write(script)
        scriptfile.flush()
        p = subprocess.Popen(['/bin/bash', scriptfile.name], stdout=subprocess.PIPE)
        out, err = p.communicate()
        # break up the output string and put it into a list
        number_of_seqs = int(len(out.decode().split("\n"))/5)
        sequence = 1
        output, seq_ID_line_numbers = [], []
        count = 0
        while sequence <= number_of_seqs:
            n = 0
            n = n + 1
            # get the line number to access the fastq information
            matchObj = re.findall( r'^[0-9]*', out.decode().split("\n")[n], re.M|re.I)
            seq_ID_line_number = (int(str(matchObj).lstrip("\[\'").rstrip("\'\]"))-2)
            # get sequence itself
            seq = out.decode().split("\n")[n].lstrip(str(matchObj))
            output.append(seq.lstrip(":"))
            # increment 3 as a fastQ file is in blocks of 5
            n = n + 3
            sequence = sequence + 1
            if get_paired == "1":
                get_fastq_format_from_ID_line_number(sequencefile, \
                                                    paired_sequence_file, \
                                                    seq_ID_line_number, \
                                                    primer_seq, \
                                                    count)
            count = count + 1
        return(output)

# get paired reads from fastQ files from the line number
#takes the line number (counting from zero) of a sequence ID, and get both fastq
#format sequences from the raw files (which count from 1)

def get_fastq_format_from_ID_line_number(forward_paired_file, \
                                        reverse_paired_file, ID_line_number, \
                                        primer_seq, count):
    with open("Forward_reads_for_assembly.fastq", 'a') as ffa, \
        open("Reverse_reads_for_assembly.fastq", 'a') as rfa:
        F_output = ""
        R_output = ""
        with open(forward_paired_file) as f, open(reverse_paired_file) as r:
            # yes, this is clunky. put a bit more thought into how to do this.
            # works in current form just janky
            for i, F_line in enumerate(f):
                if i == ID_line_number:  # ID line
                    for x in F_line.split(":")[0:9]:
                        F_output = F_output + x + ":"
                    F_output = F_output + primer_seq + "\n"
                if i == ID_line_number+1:  # sequence
                    F_output = F_output + F_line
                if i == ID_line_number+2:    # spacer
                    F_output = F_output + F_line
                if i == ID_line_number+3:   # quality info
                    F_output = F_output + F_line
            for i, R_line in enumerate(r):
                if i == ID_line_number:  # ID line
                    for x in R_line.split(":")[0:9]:
                        R_output = R_output + x + ":"
                    R_output = R_output + primer_seq + "\n"
                if i == ID_line_number+1:  # sequence
                    R_output = R_output + R_line
                if i == ID_line_number+2:    # spacer
                    R_output = R_output + R_line
                if i == ID_line_number+3:   # quality info
                    R_output = R_output + R_line
                count = count + 1
            ffa.write(F_output)
            rfa.write(R_output)

# count the occurrences of a primer sequence in a FastQ file
def count_primers(script, sequence_filecount, direction):
    with tempfile.NamedTemporaryFile(mode='w') as scriptfile:
        scriptfile.write(script)
        scriptfile.flush()
        p = subprocess.Popen(['/bin/bash', scriptfile.name], \
                                stdout=subprocess.PIPE)
        out = p.communicate()
        # break up the output string and put it into a list
        n = 0
        for x in out[0].decode():
            if x != "\n":
                if direction == "F":
                    if x != "0":
                        Fprimercountlist[sequence_filecount].append("1")
                    else:
                        Fprimercountlist[sequence_filecount].append("0")
                if direction == "R":
                    if x != "0":
                        Rprimercountlist[sequence_filecount].append("1")
                    else:
                        Rprimercountlist[sequence_filecount].append("0")
            n = n + 1

# count how many individuals data contains the primer sequence
#we have a list of lists *all_F_lists or *all_R_list, all the same length
#one list corresponds to each sequence file and each "row" is a locus
#zip over all the lists to determine which individuals have the primer
#at least x number of times (minimum individuals)

def how_many_files(group_of_lists, direction, \
                   minimum_individuals, how_many_individuals):
    position = 0
    for x in zip (*group_of_lists):
        count = 0
    # in how many individuals was it found at least once?
        for i in x:
            if i != "0":
                count = count + 1
        if int(count) > float(minimum_individuals):
            if direction == "F":
                wanted_F_primers.append(all_F_primers[position])
                how_many_individuals.append(count)
                file_contains_F_primer.append(x)
            else:
                wanted_R_primers.append(all_R_primers[position])
                file_contains_R_primer.append(x)
        position = position + 1

def get_reads(sequencefile, wanted_primers, file_contains_primer, n, \
              paired_sequence_file, output, containing_file):
    wanted_reads = []
    for x, y in zip(wanted_primers, file_contains_primer):
        # if y==1 then the primer appears in that file, so do the search
        if str(y[n:n+1]) == "('1',)":
            grep_script = "grep -nrh -B 1 -A 2 \"" + \
                            x.rstrip("\n") + "\" " + sequencefile
            result = get_seqs(grep_script, "1", sequencefile, \
                              paired_sequence_file, x)
            #if result is not None:
            if len(result) > 0:
                for sequence in result:
                    output.append(">" + x + "\n" + sequence.rstrip("\n"))
                    containing_file.append(sequencefile)
    return(output, containing_file)

def assemble_reads(forward_fastq, reverse_fastq, assembled_reads):
    # description of options:
    if not configParser.get('config_file', 'PANDAseq_exe'):
        pandaseq_command = 'pandaseq -f ' + forward_fastq + \
                            ' -r ' + reverse_fastq + ' -w ' + \
                             assembled_reads
    else:
        if configParser.get('config_file', 'PANDAseq_exe'):
            if os.path.isfile(configParser.get('config_file', 'PANDAseq_exe')):
                pandaseq_command = configParser.get('config_file', 'PANDAseq_exe') \
                    + ' -f ' + forward_fastq + ' -r ' + reverse_fastq + ' -w ' + \
                    assembled_reads
            else:
                print("ERROR............")
                print("The location of PANDAseq specified in the config file")
                print("is not correct.")
                print("Please double check and update the filepath.\n")
                quiit()

   ## redirect all pandaseqs spiel to dev/null/
    FNULL = open(os.devnull, 'w')
    subprocess.call(pandaseq_command, shell=True, stdout=FNULL, \
                    stderr=subprocess.STDOUT)

###########################################################
# MAIN PROGRAM
###########################################################
if __name__ == "__main__":
    wd = os.getcwd()

    # parse arguments
    parser = argparse.ArgumentParser(description='arguments for MiMi.py')
    parser.add_argument('-c','--config1', help='MiMi configuration file', \
                    required=True)
    args = parser.parse_args()

    # print obnoxious ascii text art
    print("\n          ~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("           __  __   _   __  __   _ ")
    print("          |  \/  | (_) |  \/  | (_)")
    print("          | \  / |  _  | \  / |  _ ")
    print("          | |\/| | | | | |\/| | | |")
    print("          | |  | | | | | |  | | | |")
    print("          |_|  |_| |_| |_|  |_| |_|")
    print("          ~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    print("Multi-Individual-Microsatellite-Identification\n\n\n")
    print("Reading config file......\n")
    print("MiMi_v0.03.py requires Python3")
    print("Currently running with the following version of Python:"+"\n")
    print(sys.version+"\n\n")
    print("Reading config file......\n")
    # Read and parse the config file
    configParser = configparser.RawConfigParser()
    configParser.read(args.config1)

    # Get number of samples
    number_of_samples = configParser.get('config_file', 'number_of_samples')
    print(number_of_samples + " samples to be analysed.\n")
    proportion_of_individuals = configParser.get('config_file', \
                                                 'proportion_of_individuals')
    pal_finder_path = configParser.get('config_file', 'pal_finder_path')
    pal_finder_config = configParser.get('config_file', 'pal_finder_config')

    ## check location of pal_finder script and config. Really don't want to
    ## find that these are missing at the end
    if not os.path.isfile(wd + "/" + pal_finder_path) or not \
           os.path.isfile(wd + "/" + pal_finder_config) == True:
        print("Fail:    Cannot find pal_finder script or pal_finder config file.")
        print("Please check paths given in the config file")
        quit()
    ## check if the user has chosen to provide absolute paths to muscle and PANDAseq
    if configParser.get('config_file', 'muscle_exe'):
        muscle_path = configParser.get('config_file', 'muscle_exe')
    ## if so, are they correct?
        if os.path.isfile(configParser.get('config_file', 'muscle_exe')):
            print("Location of muscle aligner specified in config file \
            and successfully found.")
        else:
            print("ERROR.........")
            print("The path to the MUSCLE aligner specified in the")
            print("config file is not correct.")
            print("Please double check and update the filepath.\n")
            quit()
    if configParser.get('config_file', 'PANDAseq_exe'):
        PANDAseq_path = configParser.get('config_file', 'PANDAseq_exe')
        if os.path.isfile(configParser.get('config_file', 'PANDAseq_exe')):
            print("Location of pandaseq specified in config file \
            and successfully found.")
        else:
            print("ERROR.........")
            print("The path to the PANDAseq aligner specified in the")
            print("config file is not correct.")
            print("Please double check and update the filepath.\n")
            quit()

    ##############################
    # need a method of checking that "sample number" matches the amount of
    # files entered in config"
    # if they don't match, then it just errors
    ###############################
    # Find the input files and check they exist (raw MiSeq fastQ files)
    n = 1
    check_for_duplicates = []
    # make lists of all the file names for easy access later
    list_of_R1_files = []
    list_of_R2_files = []
    list_of_pal_finder_output = []
    while n <= int(number_of_samples):
        R1input = "input" + str(n) + "_R1"
        R2input = "input" + str(n) + "_R2"
        pal_finder_output = "input" + str(n) + "_pal_finder"
        R1_file_url = configParser.get('config_file', R1input)
        list_of_R1_files.append(R1_file_url)
        R2_file_url = configParser.get('config_file', R2input)
        list_of_R2_files.append(R2_file_url)
        if R1_file_url in check_for_duplicates or \
           R2_file_url in check_for_duplicates:
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
                   ~~~~~~~~\n")
            print("Fail:    Duplicate entry in list of .fastq input files. \
                   Check config file.\n")
            print("ERROR is: Duplicate entry in list of R1 and R2 files. \
                   Every file must be unique.\n\nQuitting program....\n")
            quit()
        check_for_duplicates.append(wd + "/" + R1_file_url)
        check_for_duplicates.append(wd + "/" + R2_file_url)
        pal_finder = configParser.get('config_file', pal_finder_output)
        list_of_pal_finder_output.append(pal_finder)
    # troubleshooting file path(s):
        print("Checking sequencing files and pal_finder output exist for \
               sample " + R1input + ":\n")
        if os.path.isfile(R1_file_url) and os.path.isfile(R2_file_url) == True:
            print("Success: Found both sequencing files for sample \"" + \
                   R1input + "\"")
            if os.path.isfile(pal_finder) == True:
                print("Success: Found pal_finder output file for sample " + \
                       pal_finder_output + "\n")
            else:
                if os.path.isfile(pal_finder) == False:
                    print("\nFail: pal_finder output file(s) for samples \"" + \
                           R1input + "\" and \"" + R2input + "\" not found. \
                           Check filepath in config file.")
                    print("\nError is: Sample \"" + R1input + "\" and \"" + \
                           R2input + "\" pal_finder output file missing. \
                           Quitting program.\n")
                    quit()
        else:
            print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
                   ~~~~~\n")
            print("Fail: One or more sequence files missing for sample \"R" + \
                   str(n))
            print("\". Please check filepaths in config file, or change sample \
                   number\n\n")
            print("ERROR is: FILE MISSING.\n\nQuitting program....\n")
            quit()
        n  = n + 1
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("All input files present.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Parsing the pal_finder output files.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    # lists to hold all the primer sequences
    all_F_primers, all_R_primers = [],[]
    # lists to be generated for each input file
    Fprimerlist,Rprimerlist,seqIDs,motif = {},{},{},{}
    for x, y in zip (list_of_pal_finder_output, range (0,int(number_of_samples))):
        Fprimerlist[y],Rprimerlist[y],seqIDs[y],motif[y] = [],[],[],[]
        with open(x) as pal_finder_csv:
            pal = csv.reader(pal_finder_csv, delimiter = "\t")
            next(pal)
            for line in pal:
                seqIDs[y].append(line[0])
                motif[y].append(line[2])
                Fprimerlist[y].append(line[7])
                all_F_primers.append(line[7])
                Rprimerlist[y].append(ReverseComplement(line[9]))
                all_R_primers.append(ReverseComplement(line[9]))
    # write out combined primer lists to files
    filePath = str(os.getcwd())

    # generate a big list of all F primers and all R primers
    # also write these out to files as they are needed later
    for x in range(0, int(number_of_samples)):
        for F_primer, R_primer in zip(Fprimerlist[x], Rprimerlist[x]):
            all_F_primers.append(F_primer)
            all_R_primers.append(R_primer)

    with open(filePath + "/Fprimerlist.txt", 'w') as f, \
         open(filePath + "/Rprimerlist.txt", 'w') as r:
        for F_primer, R_primer in zip(all_F_primers, all_R_primers):
            f.write(F_primer + "\n")
            r.write(R_primer + "\n")
    f.close()
    r.close()
    # generate an empty list for each Individual to hold the primer count data:
    Fprimercountlist, Rprimercountlist = {}, {}
    all_F_lists, all_R_lists = [], []
    primer_counts = int(number_of_samples)
    for a in range (0,int(primer_counts)):
        Fprimercountlist[a] = []
        Rprimercountlist[a] = []
        all_F_lists.append(Fprimercountlist[a])
        all_R_lists.append(Rprimercountlist[a])

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Searching each individual FastQ for common primer sequences.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    # for every primer sequence, grep to count the occurrences of that sequence
    # in the fastQ files of each sample

    # Pair up every forward and reverse primer pair in a dict
    primer_pairs = dict(zip(all_F_primers, all_R_primers))

    # it selects different regions of the big Fprimerlist.txt to not bother
    # checking for primer sequences in the indiviual where they were detected
    Fprimerfile = "Fprimerlist.txt"
    Rprimerfile = "Rprimerlist.txt"
    sequencefile_count = 0
    sample_number = 1
    for Fsequencefile, Rsequencefile in zip(list_of_R1_files, list_of_R2_files):
        print((time.strftime("%H:%M:%S")))
        print("Currently searching sample: " + str(sample_number) + "\n")
        start_pos = 0
        primer_list = 0
        while primer_list < int(number_of_samples):
            end_pos = start_pos + len(Fprimerlist[primer_list]) - 1
            diff = end_pos - start_pos
            if int(primer_list) != int(sequencefile_count):
                # then it is a primer -> sequence file comparison that need to make
                Fgrep_script = "cat " + Fprimerfile + " 2>/dev/null | head -n " + \
                                str(end_pos+1) + " | tail -n " + str(diff+1) + \
                                " | while read line ; do grep -c $line " + \
                                Fsequencefile + "; done"
                count_primers(Fgrep_script, sequencefile_count, "F")
                start_pos = start_pos + len(Fprimerlist[primer_list])
                primer_list = int(primer_list) + 1
            else:

                # then we are searching for primer sequences in the individual
                # they were discovered
                # so it can be skipped and recorded as '1'
                # fill the primer count data with '1'
                count = start_pos
                n = 0
                while int(count) < int(end_pos+1):
                    Fprimercountlist[sequencefile_count].append("1")
                    Rprimercountlist[sequencefile_count].append("1")
                    count = count + 1
                    n = n + 1
                start_pos = start_pos + len(Fprimerlist[primer_list])
                primer_list = int(primer_list) + 1

        sequencefile_count = int(sequencefile_count ) + 1
        sample_number = sample_number + 1

    ## now that we have all the sequence counts, we can go through the lists
    # and get a list of just those primers which occur in > 50% of individuals
    # get a list of all the primers which appear in >50% of individuals
    wanted_F_primers, wanted_R_primers = [], []
    file_contains_F_primer, file_contains_R_primer = [], []
    position = 0
    minimum_individuals = float(number_of_samples) * \
                          float(proportion_of_individuals)
    how_many_individuals = []
    how_many_files(all_F_lists, "F", minimum_individuals, how_many_individuals)
    #how_many_files(all_R_lists, "R", minimum_individuals, how_many_individuals)

    with open("wanted_F_primers.txt", 'w') as wantedF, \
         open("wanted_R_primers.txt", 'w') as wantedR:
        for x in wanted_F_primers:
            wantedF.write(x + "\n")
        for y in wanted_R_primers:
            wantedR.write(y + "\n")
    # put the primer seqs and how many inds they were found in, into a dict
    which_individuals_F = dict(zip(wanted_F_primers, how_many_individuals))

    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("Merging all reads containing common primer sequences.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    # this produces an F fastq and a R fastq for assembly.
    # they are open as append so need to be removed first, if they already exist
    if os.path.exists("Forward_reads_for_assembly.fastq"):
        os.remove("Forward_reads_for_assembly.fastq")
    if os.path.exists("Reverse_reads_for_assembly.fastq"):
        os.remove("Reverse_reads_for_assembly.fastq")

    wanted_F_reads, wanted_R_reads, F_sequence_file = [], [], []
    R_sequence_file, containing_file = [], []
    n = 0
    for sequencefile, paired_sequence_file in \
                       zip(list_of_R1_files, list_of_R2_files):
        print(time.strftime("%H:%M:%S"))
        print("Extracting reads from file: " + sequencefile + "\n")
        get_reads(sequencefile, wanted_F_primers, file_contains_F_primer, n, \
                  paired_sequence_file, wanted_F_reads, containing_file)
        n = n + 1

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print (time.strftime("%H:%M:%S"))
    print("Finished scanning samples for common primer sequences")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print (time.strftime("%H:%M:%S"))
    print("Creating config file for pal_finder")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    # create a unique pal_finder config file for each
    # check for file of pal_finder configs, delete and re-create if necessary

    if os.path.isdir(wd + "/MiMi_output/pal_finder_files"):
        shutil.rmtree(wd + "/MiMi_output/pal_finder_files")
    if os.path.isdir(wd + "/MiMi_output"):
        shutil.rmtree(wd + "/MiMi_output")
    os.mkdir(wd + "/MiMi_output")
    os.mkdir(wd + "/MiMi_output/pal_finder_files")
    output_path = wd + "/MiMi_output/pal_finder_files/"
    assemble_reads("Forward_reads_for_assembly.fastq", \
                    "Reverse_reads_for_assembly.fastq", \
                    "Assembled_reads.fasta" )

    wd = os.getcwd()
    if os.path.isdir(wd + "/MiMi_output/Alignments"):
        shutil.rmtree(wd + "/MiMi_output/Alignments")
    os.mkdir(wd + "/MiMi_output/Alignments")
    output_path = wd + "/MiMi_output/Alignments/"

    fasta1 = SeqIO.parse("Assembled_reads.fasta",'fasta')
    forward_reads = SeqIO.parse("Forward_reads_for_assembly.fastq",'fastq')

    ## filter the list of source_files to just those reads which made it through
    ## the assembly process (this should be the vast majority of them)

    ## go through and get the reads IDs that were successfully assembled.
    list_of_assembled_reads = []
    with open("Assembled_reads.fasta", 'r') as assembled_reads_fasta:
        for line in assembled_reads_fasta:
            if line.startswith(">"):
                list_of_assembled_reads.append(line.rsplit(":", 1)[0].lstrip(">"))
    pre_assembled_reads = []

    ## need to get the IDs from the forward_for_assembly file and add them into this list
    for record in forward_reads:
        pre_assembled_reads.append(record.id)

    ## filter the list of containing files to just the entries that made it through assembly.
    wanted_containing_files = []
    for pre_assembled, source_file in zip(pre_assembled_reads, containing_file):
        if pre_assembled in list_of_assembled_reads:
            wanted_containing_files.append(source_file)

    list_of_assembled_sequences = []
    list_of_assembled_IDs = []
    step = 2
    with open("Assembled_reads.fasta", 'r') as assembled_reads_file:
        for number, line in enumerate(assembled_reads_file):
            if number % step != 0:
                list_of_assembled_sequences.append(line)
            else:
                list_of_assembled_IDs.append(line.rstrip("\n").split(";")[0])

    wanted = set()
    for ID, seq, containing_file in zip(list_of_assembled_IDs, list_of_assembled_sequences, wanted_containing_files):
        sequence_ID = (ID.split(":")[7])
        sequence = str(seq.rstrip("\n"))
        filePath = str(os.getcwd()) + '/MiMi_output/Alignments/%s.fasta' % \
                       (sequence_ID)
        if sequence_ID in wanted:
            if os.path.exists(filePath):
                with open(filePath, 'a') as f:
                    f.write(">" + containing_file.split("/")[len(containing_file.split("/"))-1] + "\n")
                    f.write(sequence + "\n")
        else:
            with open(filePath, 'w') as f:
                ## write the forward primer into the alignment
                f.write(">F_primer:_" + sequence_ID + "\n")
                f.write(sequence_ID + "\n")
                f.write(">" + containing_file.split("/")[len(containing_file.split("/"))-1] + "\n")
                f.write(sequence + "\n")
        wanted.add(sequence_ID)

    # even though some primer sequences are found in the data files enough times
    # to get through the process this far. it is possible that not every read
    # is sufficiently high quality to form an assembly.
    # this can result in the situation where alignments with only a single
    # sequence get through

    # check length of files in Alignments folder
    # remove anything with only one sequence
    for filename in os.listdir(wd + "/MiMi_output/Alignments"):
        if not file_len(wd + "/MiMi_output/Alignments/" + filename) > 2:
            os.remove(wd + "/MiMi_output/Alignments/" + filename)
    output_path = wd

    copy_and_rename_file(pal_finder_config, output_path + \
                         "pal_finder_config_file.txt")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "findPrimers 1", "findPrimers 0")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "platform Illumina", "platform 454")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "inputFormat fastq", "inputFormat fasta")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "pairedEnd  1","pairedEnd  0" )
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "input454reads  test/data/454_All_python.fna", \
                    "input454reads " + wd + "/Assembled_reads.fasta")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "MicrosatSumOut  test/output/test_microsat_summary.txt", \
                    "MicrosatSumOut " + output_path + \
                    "/MiMi_output/pal_finder_files/pal_finder_summary_out.txt")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "PALsummaryOut  test/output/test_PAL_summary.txt", \
                    "PALsummaryOut " + output_path + \
                    "/MiMi_output/pal_finder_files/pal_finder_PAL_summary.txt")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "2merMinReps 	6", "2merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "3merMinReps 	0", "3merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "4merMinReps 	0", "4merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "5merMinReps 	0", "5merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt", \
                    "6merMinReps 	0", "6merMinReps 	6")
    pal_finder_command = 'perl ' + pal_finder_path + " " + \
                          output_path + "pal_finder_config_file.txt"
    print("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print (time.strftime("%H:%M:%S"))
    print("Running pal_finder.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    subprocess.call(pal_finder_command, shell=True)

    print("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print (time.strftime("%H:%M:%S"))
    print("pal_finder complete.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

    #### parse the pal_finder output to find the variable loci
    # get a list of unique IDs (the forward primer sequences)
    with open(wd + "/MiMi_output/pal_finder_files/pal_finder_PAL_summary.txt") \
        as pf:
        unique_primers = set()
        allele_count = []
        unique_alleles = []
        next(pf)
        for line in pf:
            motifs = (line.split("\t")[3])
            if (len(motifs.split(" "))-1) == 1:
                ## it doesn't seem to be getting here
                ID = line.split("\t")[0]
                unique_primers.add(ID.split(":")[7].split(";")[0])
                allele_count.append(0)
                unique_alleles.append("Motifs: ")
    list_unique_primers = list(unique_primers)

    ## go through and get unique alleles associated with each primer sequence
    all_data = []
    with open(wd + "/MiMi_output/pal_finder_files/pal_finder_PAL_summary.txt") \
        as pf:
        for line in pf:
            motifs = (line.split("\t")[3])
            if (len(motifs.split(" "))-1) == 1:
                ID = line.split("\t")[0]
                for i, primer_line in enumerate(list_unique_primers):
                    if ID.split(":")[7].split(";")[0] == primer_line:
                        allele_count[i] = allele_count[i] + 1
                        unique_alleles[i] = unique_alleles[i] + motifs

    # merge multiple lists
    for x, y, z in zip(unique_primers, allele_count, unique_alleles):
        all_data.append("\t".join([x, str(y), z]))

    # remove any lines which have reported multiple motifs
    single_motif_only = []
    for row in all_data:
        unique = set()
        motifs = row.split("\t")[2].split(":")[1].lstrip(" ").rstrip(" ")
        for single_motif in motifs.split(" "):
            unique.add(single_motif.split("(")[0])
        if (len(unique) == 1 and len(motifs.split(" ")) > 1):
            single_motif_only.append(row.rstrip(" "))

    # calcualte difference in size between the biggest numof repeats and smallest
    diff_in_motif_size = []
    for row in single_motif_only:
        number_of_repeats = []
        for x in row.split("\t")[2].split(" ")[1:]:
            for result in (re.findall(r'\d+', x)):
                number_of_repeats.append(int(result))
        diff_in_motif_size.append((max(number_of_repeats) - \
                                   min(number_of_repeats)))
    # rank by difference in largest and samllest motifs
    ranked_output = []
    count = 0
    if len(diff_in_motif_size) > 0:
        while count <= max(diff_in_motif_size):
            for x, y in zip(single_motif_only, diff_in_motif_size):
                if y == count:
                    ranked_output.insert(0, x.replace("Motifs: ", " ") + \
                                         "\t" + str(y))
            count = count + 1
        # write out final output
        with open("MiMi_output/MiMi_output_all_loci.txt", 'w') as final_output:
            final_output.write("Forward_primer_seq(5'-3')\tReverse_primer_seq(5'-3')\t\
                               Number_of_alleles\tFound_in_individuals\t\
                               Alleles_present\tSize-Range\n")
            for x in ranked_output:
                final_output.write(x.split("\t")[0] + "\t" + \
                             ReverseComplement(primer_pairs[x.split("\t")[0]]) \
                             + "\t" + x.split("\t")[1] + "\t" + \
                             str(which_individuals_F[x.split("\t")[0]]) + \
                             "\t" + x.split("\t")[2] + "\t" + \
                             x.split("\t")[3] + "\n")
        final_output.close()
    else:
        print("Something went wrong.MiMi has not found any microsatellites" \
                " in the sequence data, which occur in multiple individuals.")
        print("Please check all your input files. failing that, please contact"\
        " the author.")

    # some tidying up of temporary files

    if os.path.exists("Forward_reads_for_assembly.fastq"):
        os.remove("Forward_reads_for_assembly.fastq")
    if os.path.exists("Reverse_reads_for_assembly.fastq"):
        os.remove("Reverse_reads_for_assembly.fastq")
    if os.path.exists("Fprimerlist.txt"):
        os.remove("Fprimerlist.txt")
    if os.path.exists("Rprimerlist.txt"):
        os.remove("Rprimerlist.txt")
    if os.path.exists("pal_finder_input_file.fasta"):
        os.remove("pal_finder_input_file.fasta")
    if os.path.exists("wanted_F_primers.txt"):
        os.remove("wanted_F_primers.txt")
    if os.path.exists("wanted_R_primers.txt"):
        os.remove("wanted_R_primers.txt")
    if os.path.exists("Assembled_reads.fasta"):
        os.remove("Assembled_reads.fasta")
    if os.path.isdir(wd + "/MiMi_output/pal_finder_files"):
        shutil.rmtree(wd + "/MiMi_output/pal_finder_files")
    if os.path.exists(output_path + "pal_finder_config_file.txt"):
        os.remove(output_path + "pal_finder_config_file.txt")

    ### additional filtering added as the suggestion by peer review

    ##### align all the fasta files in the Alignment directory using MUSCLE
    wd = os.getcwd()
    for (dirpath, dirnames, filenames) in walk(wd + "/MiMi_output/Alignments"):
        for MSA in filenames:
            if not configParser.get('config_file', 'muscle_exe'):
                align_command = "muscle -in " + wd + "/MiMi_output/Alignments/" \
                + MSA + " -out " + wd + "/MiMi_output/Alignments/" + MSA \
                + ".aln -quiet"
            else:
                if configParser.get('config_file', 'muscle_exe'):
                    if os.path.isfile(configParser.get('config_file', 'muscle_exe')):
                        align_command = configParser.get('config_file', 'muscle_exe') \
                        + " -in " + wd + "/MiMi_output/Alignments/" + MSA \
                        + " -out " + wd + "/MiMi_output/Alignments/" + MSA \
                        + ".aln -quiet"
                    else:
                        print("ERROR.........")
                        print("The path to the MUSCLE aligner specified in the")
                        print("config file is not correct.")
                        print("Please double check and update the filepath.\n")
                        quit()
            subprocess.call(align_command, shell=True)

    ## trim everything previous to the position of the forward primer in the alignment
    ### get the position where the left side trim should occur (ie. where  F primer sits)
    left_side_trim = []
    with open("MiMi_output/MiMi_output_all_loci.txt", 'r') as mimi_output:
        next(mimi_output)
        for line in mimi_output:
            record_dict = SeqIO.parse("MiMi_output/Alignments/" \
                    + line.split("\t")[0] + ".fasta.aln", "fasta")
            for record in record_dict:
                ## find the entry in the MSA which is just the F primer
                if not ".fastq" in record.id:
                    f_positional_counter = 0
                    for f_position in record.seq:
                        if f_position != "-":
                            left_side_trim.append(f_positional_counter)
                            break
                        else:
                            f_positional_counter = f_positional_counter + 1
    mimi_output.close()

    ## perform the trim
    with open("MiMi_output/MiMi_output_all_loci.txt", 'r') as mimi_output:
        next(mimi_output)
        for line, trim_pos in zip(mimi_output, left_side_trim):
            with open("MiMi_output/Alignments/" + line.split("\t")[0] \
            + ".trimmed", 'w') as trimmed_fasta:
                record_dict = SeqIO.parse("MiMi_output/Alignments/" \
                            + line.split("\t")[0] + ".fasta.aln", "fasta")
                for record in record_dict:
                    trimmed_fasta.write(">" + record.id + "\n")
                    trimmed_fasta.write(str(record.seq[trim_pos:]) + "\n")
    mimi_output.close()

    full_length_perfect = []
    conserved_reverse_primer = []
    rev_primer_missing = []
    partial_primer_match = []
    with open("MiMi_output/MiMi_output_all_loci.txt", 'r') as mimi_output:
        next(mimi_output)
        non_perfect = 0
        partial_match = 0
        for line in mimi_output:
            counter = 0
            record_dict = SeqIO.parse("MiMi_output/Alignments/" \
                            + line.split("\t")[0] + ".trimmed", "fasta")
            ## get the position of the reverse primer in the alignment
            for record in record_dict:
                if ReverseComplement(line.split("\t")[1]) in record.seq:
                    rev_position = str(record.seq).index(ReverseComplement(line.split("\t")[1]))

            ## handle non-perfect matches, and choose which loci to discard if have mis-matches etc.
            record_dict2 = SeqIO.parse("MiMi_output/Alignments/" \
                            + line.split("\t")[0] + ".trimmed", "fasta")

            for record in record_dict2:
                ### look at the section of alignment where the rev primer falls
                ## if string contains ACGT - ie. is not just an empty piece of alignment
                if bool(re.search('[ACGT]', str(record.seq[rev_position:rev_position \
                                            + len(line.split("\t")[1])]))):
                    ## if it is a full, length perfect match:
                    if ReverseComplement(line.split("\t")[1]) in record.seq[rev_position:rev_position + len(line.split("\t")[1])]:
                        #counter = counter + 1  # not sure that counter is doing anything here
                        non_perfect = 0  # this should already be zero, in fact
                    ## otherwise flag it as non-perfect
                    else:
                        non_perfect = 1
                    if "-" in str(record.seq[rev_position:rev_position + len(line.split("\t")[1])]):
                        if str(record.seq[rev_position:rev_position + len(line.split("\t")[1])].lstrip("-").rstrip("-")) in ReverseComplement(line.split("\t")[1]):
                            partial_match = 1

            ## build up lists depending on the quality of the alignment at the reverse primer
            if non_perfect == 1:
                full_length_perfect.append("0")
                non_perfect = 0
                if counter > 1:
                    conserved_reverse_primer.append("0")
            else:
                full_length_perfect.append("1")
                if counter > 1:
                    conserved_reverse_primer.append("1")
            if counter <= 1:
                conserved_reverse_primer.append("0")
            if partial_match == 1:
                partial_primer_match.append("1")
                partial_match = 0
            else:
                partial_primer_match.append("0")

    filtered_output = []
    monomorphic_loci = []
    low_qual_alignment_loci = []
    low_qual_to_remove = []
    priority2_output = []
    priority3_output = []
    mutations_in_rev_primer = []
    total_loci = float(0)
    retain_low_qual_alignments = []
    with open("MiMi_output/MiMi_output_all_loci.txt", 'r') as mimi_output:
        next(mimi_output)
        for line, rev_primer, full_length_cons, partial in zip(mimi_output, conserved_reverse_primer, full_length_perfect, partial_primer_match):
            total_loci = total_loci + 1
            record_dict = SeqIO.parse("MiMi_output/Alignments/" + line.split("\t")[0] + ".trimmed", "fasta")
            average = []
            for record in record_dict:
                ### measure the overall alignment "quality" - crude but effective to remove messy alignments
                if not line.split("\t")[0] in record.id:
                    ## remove trailing and leading spaces in alignment as not indicative of low qual
                    curr_seq = str(record.seq).rstrip("-").lstrip("-")
                    # count spaces remaining in alignment - measure of gappiness
                    average.append(curr_seq.count("-"))

            ## anything where this score is < user-definable quality score (default: 3)
            # 'is a "lower quality" alignment
            if float(sum(average))/float(len(average)) < int(configParser.get('config_file', 'overall_alignment_qual_score')):
                low_qual_alignment_loci.append("1")
            else:
                low_qual_alignment_loci.append("0")
    mimi_output.close()

    with open("MiMi_output/MiMi_output_filtered_loci.txt", 'w') as filtered_output:
        filtered_output.write("Forward_primer_seq(5'-3')\tReverse_primer_seq(5'-3')\tNumber_of_alleles\tFound_in_individuals\tAlleles_present\tSize_range\n")
        filtered_output.write("HIGH QUALITY LOCI\n")
        with open("MiMi_output/MiMi_output_all_loci.txt", 'r') as mimi_output:
            next(mimi_output)
            high_qual_counter = 0
            OK_qual_counter = 0
            for line, rev_primer, full_length_cons, partial, low_qual in zip(mimi_output, conserved_reverse_primer, full_length_perfect, partial_primer_match, low_qual_alignment_loci):
                ## remove monomorphic loci
                if line.split("\t")[5] == "0\n" or low_qual != "1":
                    if line.split("\t")[5] == "0\n":
                        monomorphic_loci.append(line)
                    ## remove any "low quality alignments"
                    elif low_qual != "1":
                        low_qual_to_remove.append(line)
                        retain_low_qual_alignments.append(line)
                ## loci with conserved rev primer written to file
                elif rev_primer == "1":
                    filtered_output.write(line)
                    high_qual_counter = high_qual_counter + 1
                ## anything else goes to priority2 output, for easy ranking
                elif rev_primer != "1":
                    OK_qual_counter = OK_qual_counter + 1
                    priority2_output.append(line.rstrip("\n") + "\t" \
                            + rev_primer + "\t" + full_length_cons \
                            + "\t" + partial + "\t" + low_qual + "\n")
        mimi_output.close()
        filtered_output.write("GOOD QUALITY LOCI\n")
        for priority2 in priority2_output:
            ## if it is a full length match, or a partial match
            if priority2.split("\t")[7] == "1" or priority2.split("\t")[8] == "1":
                if priority2.split("\t")[7] == "1":
                    filtered_output.write("\t".join(priority2.split("\t")[0:6])+"\n")
                elif priority2.split("\t")[8] == "1":
                    priority3_output.append(priority2.rstrip("\n") + "\t" \
                            + rev_primer + "\t" + full_length_cons + "\t" \
                            + partial + "\t" + low_qual + "\n")
            else:
                ## otherwise, we specifically have a mutation in the primer and it should be avoided
                mutations_in_rev_primer.append("\t".join(priority2.split("\t")[0:6])+"\n")
        for priority3 in priority3_output:
            filtered_output.write("\t".join(priority3.split("\t")[0:6])+"\n")
    OK_qual_counter = OK_qual_counter-len(mutations_in_rev_primer)

    ### write the log file with details of which loci were removed and why
    with open("MiMi_output/MiMi_loci_filter.log", 'w') as filter_log:
        filter_log.write(time.strftime("%H:%M:%S") + "\n")
        filter_log.write("Statistics relating to numbers of loci filtered out by automated MiMi quality filters\n\n")
        filter_log.write("Total number of putative microsatellite markers detected by MiMi: " + str(total_loci) + "\n")
        filter_log.write("Number of putative markers classified as HIGH quality: " \
                        + str(high_qual_counter)
                        + " (" + str("{:.1f}".format((high_qual_counter/total_loci)*100)) \
                            + "% passed)" + "\n")
        filter_log.write("Number of putative markers classified as GOOD quality: " \
                        + str(OK_qual_counter) \
                        + " (" + str("{:.1f}".format((OK_qual_counter/total_loci)*100)) \
                            + "% passed)" + "\n")
        filter_log.write("Markers determined to be low(er) quality alignments: " \
        + str(len(retain_low_qual_alignments)) + \
        " (" + str("{:.1f}".format((len(retain_low_qual_alignments)/total_loci)*100)) \
                            + "% filtered)" + "\n")
        filter_log.write("Markers determined to be monomorphic: " \
                            + str(len(monomorphic_loci)) + \
                " (" + str("{:.1f}".format((len(monomorphic_loci)/total_loci)*100)) \
                                    + "% filtered)" + "\n")
        filter_log.write("Mutations in reverse primer binding site: " \
        + str(len(mutations_in_rev_primer)) + \
        " (" + str("{:.1f}".format((len(mutations_in_rev_primer)/total_loci)*100)) \
                            + "% filtered)" + "\n")
        filter_log.write("\n\nThe filtered loci are printed below for your reference: \n\n")
        filter_log.write("LOW QUALITY ALIGNMENTS\n")
        filter_log.write("Foward_primer_seq\tReverse_primer_seq\tNumber_of_alleles\tFound_in_individuals\tAlleles_present\tSize_range\n")
        for low_qual_alignment in retain_low_qual_alignments:
            filter_log.write(low_qual_alignment)
        filter_log.write("\n\n")
        filter_log.write("MONOMORPHIC LOCI\n")
        filter_log.write("Foward_primer_seq\tReverse_primer_seq\tNumber_of_alleles\tFound_in_individuals\tAlleles_present\tSize_range\n")
        for monomorph in monomorphic_loci:
            filter_log.write(monomorph)
        filter_log.write("\n\n")
        filter_log.write("MUTATIONS IN REVERSE PRIMER SEQUENCE\n")
        filter_log.write("Foward_primer_seq\tReverse_primer_seq\tNumber_of_alleles\tFound_in_individuals\tAlleles_present\tSize_range\n")
        for missing_rev in mutations_in_rev_primer:
            filter_log.write(missing_rev)
    print("\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print(time.strftime("%H:%M:%S"))
    print("\n\nSuccessfully completed MiMi analysis.")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
