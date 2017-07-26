#!/usr/bin/python
import ConfigParser, subprocess, os, time, csv, Bio, argparse, re, shutil, time, tempfile
from Bio import SeqIO
####
# this is the version that appears to be working on 25_07_17
# check with several sets of data and assembly, filtering etc all appears to be working
#
#
# Additional features to build at some point:
#
# make it tell you which sequence file each read came from
#
#
# searching for reads containing the reverse primer will almost certainly improve the rate at
# which conserved reverse primer regions are discovered. There will be a performance hit
# but I think it is worth it.

# it would be really useful if the primer sequences were automatically output into
# the alignments rather than having to be added by hand later

# similarly, if the alignment were done automatically
# some sort of quality filter as to how good the primers score would also be nice.
# (no idea how this would work in practice though)

###
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
    with tempfile.NamedTemporaryFile() as scriptfile:
        scriptfile.write(script)
        scriptfile.flush()
        p = subprocess.Popen(['/bin/bash', scriptfile.name], stdout=subprocess.PIPE)
        out, err = p.communicate()
        # break up the output string and put it into a list
        number_of_seqs = int(len(out.split("\n")))/5
        sequence = 1
        output, seq_ID_line_numbers = [], []
        count = 0
        while sequence <= number_of_seqs:
            n = 0
            # get sequence ID - working but not currently in use
            #matchObj = re.findall( r'^[0-9]*-', out.split("\n")[n], re.M|re.I)
            n = n + 1
            # get the line number to access the fastq information
            matchObj = re.findall( r'^[0-9]*', out.split("\n")[n], re.M|re.I)
            seq_ID_line_number = (int(str(matchObj).lstrip("\[\'").rstrip("\'\]"))-2)
            # get sequence itself
            seq = out.split("\n")[n].lstrip(str(matchObj))
            output.append(seq.lstrip(":"))
            # increment 3 as a fastQ file is in blocks of 5
            n = n + 3
            sequence = sequence + 1
            if get_paired == "1":
                get_fastq_format_from_ID_line_number(sequencefile, paired_sequence_file, seq_ID_line_number, primer_seq, count)
            count = count + 1
        return(output)

# get paired reads from fastQ files from the line number

#takes the line number (counting from zero) of a sequence ID, and get both fastq
#format sequences from the raw files (which count from 1)

def get_fastq_format_from_ID_line_number(forward_paired_file, reverse_paired_file, ID_line_number, primer_seq, count):
    with open("Forward_reads_for_assembly.fastq", 'a') as ffa, open("Reverse_reads_for_assembly.fastq", 'a') as rfa:
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
    with tempfile.NamedTemporaryFile() as scriptfile:
        scriptfile.write(script)
        scriptfile.flush()
        p = subprocess.Popen(['/bin/bash', scriptfile.name], stdout=subprocess.PIPE)
        out, err = p.communicate()
        # break up the output string and put it into a list
        n = 0
        for x in out:
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

def how_many_files(group_of_lists, direction, minimum_individuals, how_many_individuals):
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

def get_reads(sequencefile, wanted_primers, file_contains_primer,n, paired_sequence_file, output):
    wanted_reads = []
    for x, y in zip(wanted_primers, file_contains_primer):
        # if y==1 then we know that the primer appears in that file, so do the search
        if str(y[n:n+1]) == "('1',)":
            grep_script = "grep -nr -B 1 -A 2 \"" + x.rstrip("\n") + "\" " + sequencefile
            result = get_seqs(grep_script, "1", sequencefile, paired_sequence_file, x)
            if result is not None:
                for sequence in result:
                    output.append(">" + x + "\n" + sequence.rstrip("\n"))
    return(output)

def assemble_reads(forward_fastq, reverse_fastq, assembled_reads):
    # description of options:
    pandaseq_command = 'pandaseq -f ' + forward_fastq + \
                        ' -r ' + reverse_fastq + ' -w ' + \
                         assembled_reads
    subprocess.call(pandaseq_command, shell=True)

###########################################################
# MAIN PROGRAM
###########################################################
if __name__ == "__main__":
    wd = os.getcwd()
    # print obnoxious ascii text art
    print "\n          ~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "           __  __   _   __  __   _ "
    print "          |  \/  | (_) |  \/  | (_)"
    print "          | \  / |  _  | \  / |  _ "
    print "          | |\/| | | | | |\/| | | |"
    print "          | |  | | | | | |  | | | |"
    print "          |_|  |_| |_| |_|  |_| |_|"
    print "          ~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    print "Multi-Individual-Microsatellite-Identification\n\n\n"
    time.sleep(1)
    print "Reading config file......\n"
    # Read and parse the config file
    time.sleep(1)
    configParser = ConfigParser.RawConfigParser()
    configFilePath = r'config.txt'
    configParser.read(configFilePath)
    # Get number of samples
    number_of_samples = configParser.get('config_file', 'number_of_samples')
    print number_of_samples + " samples to be analysed."
    proportion_of_individuals = configParser.get('config_file', 'proportion_of_individuals')
    pal_finder_script = configParser.get('config_file', 'pal_finder_path')
    pal_finder_config = configParser.get('config_file', 'pal_finder_config')

    if not os.path.isfile(pal_finder_script) and os.path.isfile(pal_finder_config) == True:
        print "Cannot find pal_finder script or pal_finder config file."
        print "Please check paths given in the config file"
        quit()
    ##############################
    # need a method of checking that "sample number" matches the amount of files entered in config"
    # if they don't match, then it just errors
    ###############################
    time.sleep(1)
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
        if R1_file_url in check_for_duplicates or R2_file_url in check_for_duplicates:
            print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
            print "Fail:    Duplicate entry in list of .fastq input files. Check config file.\n"
            print "ERROR is: Duplicate entry in list of R1 and R2 files. Every file must be unique.\n\nQuitting program....\n"
            quit()
        check_for_duplicates.append(R1_file_url)
        check_for_duplicates.append(R2_file_url)
        pal_finder = configParser.get('config_file', pal_finder_output)
        list_of_pal_finder_output.append(pal_finder)
    # troubleshooting file path(s):
        print "Checking sequencing files and pal_finder output exist for sample " + R1input + ":\n"
        if os.path.isfile(R1_file_url) and os.path.isfile(R2_file_url) == True:
            print "Success: Found both sequencing files for sample \"" + R1input + "\""
            if os.path.isfile(pal_finder) == True:
                print "Success: Found pal_finder output file for sample " + pal_finder_output + "\n"
            else:
                if os.path.isfile(pal_finder) == False:
                    print "\nFail: pal_finder output file(s) for samples \"" + R1input + "\" and \"" + R2input + "\" not found. Check filepath in config file."
                    print "\nError is: Sample \"" + R1input + "\" and \"" + R2input + "\" pal_finder output file missing. Quitting program.\n"#
                    quit()
        else:
            print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
            print "Fail:    One or more sequence files missing for sample \"R" + str(n)
            print "\". Please check filepaths in config file, or change sample number\n\n"
            print "ERROR is: FILE MISSING.\n\nQuitting program....\n"
            quit()
        n  = n + 1
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "All input files present."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Parsing the pal_finder output files."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"

    # lists to hold all the primer sequences
    all_F_primers, all_R_primers = [],[]

    # lists to be generated for each input file
    #Fprimerlist,Rprimerlist,Fprimercount,Rprimercount,seqIDs,motif = {},{},{},{},{},{}
    Fprimerlist,Rprimerlist,seqIDs,motif = {},{},{},{}
    for x, y in zip (list_of_pal_finder_output, range (0,int(number_of_samples))):
        #Fprimerlist[y],Rprimerlist[y],Fprimercount[y],Rprimercount[y],seqIDs[y],motif[y] = [],[],[],[],[],[]
        Fprimerlist[y],Rprimerlist[y],seqIDs[y],motif[y] = [],[],[],[]
        with open(x) as pal_finder_csv:
            pal = csv.reader(pal_finder_csv, delimiter = "\t")
            next(pal)
            for line in pal:
                seqIDs[y].append(line[0])
                motif[y].append(line[6])
                Fprimerlist[y].append(line[1])
                all_F_primers.append(line[1])
                Rprimerlist[y].append(ReverseComplement(line[3]))
                all_R_primers.append(ReverseComplement(line[3]))
    # write out combined primer lists to files
    filePath = str(os.getcwd())

    # generate a big list of all F primers and all R primers
    # also write these out to files as they are needed by Blast later
    with open(filePath + "/Fprimerlist.txt", 'w') as f, open(filePath + "/Rprimerlist.txt", 'w') as r:
        for x in range(0, int(number_of_samples)):
            for F_primer, R_primer in zip(Fprimerlist[x], Rprimerlist[x]):
                all_F_primers.append(F_primer)
                all_R_primers.append(R_primer)
                f.write(F_primer + "\n")
                r.write(R_primer + "\n")

    # generate an empty list for each Individual to hold the primer count data:
    Fprimercountlist, Rprimercountlist = {}, {}
    all_F_lists, all_R_lists = [], []
    primer_counts = int(number_of_samples)
    for a in range (0,int(primer_counts)):
        Fprimercountlist[a] = []
        Rprimercountlist[a] = []
        all_F_lists.append(Fprimercountlist[a])
        all_R_lists.append(Rprimercountlist[a])

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Searching each individual FastQ for common primer sequences."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    # for every primer sequence, use grep to count the occurrences of that sequence
    # in the fastQ files of each sample

    # Pair up every forward and reverse primer pair in a dict
    primer_pairs = dict(zip(all_F_primers, all_R_primers))

    # it selects different regions of the big Fprimerlist.txt file so as to not bother
    # checking for primer sequences in the indiviual where they were first detected
    Fprimerfile = "Fprimerlist.txt"
    Rprimerfile = "Rprimerlist.txt"
    sequencefile_count = 0
    sample_number = 1
    for Fsequencefile, Rsequencefile in zip(list_of_R1_files, list_of_R2_files):
        print (time.strftime("%H:%M:%S"))
        print "Currently searching sample: " + str(sample_number) + "\n"
        start_pos = 0
        primer_list = 0
        while primer_list < int(number_of_samples):
            end_pos = start_pos + len(Fprimerlist[primer_list]) - 1
            diff = end_pos - start_pos
            if int(primer_list) != int(sequencefile_count):
                # then it is a primer -> sequence file comparison that we need to make
                Fgrep_script = "cat " + Fprimerfile + " | head -n " + \
                                str(end_pos+1) + "| tail -n " + str(diff+1) + \
                                " | while read line ; do grep -c $line " + \
                                Fsequencefile + "\n" + "done"
                Rgrep_script = "cat " + Rprimerfile + " | head -n " + \
                                str(end_pos+1) + "| tail -n " + str(diff+1) + \
                                " | while read line ; do grep -c $line " + \
                                Rsequencefile + "\n" + "done"
                count_primers(Fgrep_script, sequencefile_count, "F")
                #count_primers(Rgrep_script, sequencefile_count, "R")
                start_pos = start_pos + len(Fprimerlist[primer_list])
                primer_list = int(primer_list) + 1
            else:
                # then we are searching for primer sequences in the individual they were discovered
                # so it can be skipped and recorded as '1'
                # fill the primer count data with '1'
                count = start_pos
                n = 0
                while int(count) < int(end_pos+1):
                    Fprimercountlist[sequencefile_count].append("1")
                    count = count + 1
                    n = n + 1
                start_pos = start_pos + len(Fprimerlist[primer_list])
                primer_list = int(primer_list) + 1
        sequencefile_count = int(sequencefile_count ) + 1
        sample_number = sample_number + 1

    ## now that we have all the sequence counts, we can go through the lists
    # and get a list of just those primers which occur in > 50% of individuals

    # get a list of all the primers which appear in >50% of individuals
    wanted_F_primers, wanted_R_primers, file_contains_F_primer, file_contains_R_primer = [], [], [], []
    position = 0
    minimum_individuals = float(number_of_samples) * float(proportion_of_individuals)
    how_many_individuals = []
    how_many_files(all_F_lists, "F", minimum_individuals, how_many_individuals)
    #how_many_files(all_R_lists, "R", minimum_individuals, how_many_individuals)
    with open("wanted_F_primers.txt", 'w') as wantedF:
        for x in wanted_F_primers:
            wantedF.write(x + "\n")

    # put the primer seqs and how many individuals they were found in, into a dict
    which_individuals = dict(zip(wanted_F_primers, how_many_individuals))


    print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Merging all reads containing common primer sequences."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"

    # this produces an F fastq and a R fastq for assembly.
    # they are open as append so need to be removed first, if they already exist
    if os.path.exists("Forward_reads_for_assembly.fastq"):
        os.remove("Forward_reads_for_assembly.fastq")
    if os.path.exists("Reverse_reads_for_assembly.fastq"):
        os.remove("Reverse_reads_for_assembly.fastq")

    wanted_F_reads, wanted_R_reads = [], []
    n = 0
    for sequencefile, paired_sequence_file in zip(list_of_R1_files, list_of_R2_files):
        print (time.strftime("%H:%M:%S"))
        print "Extracting reads from file: " + sequencefile + "\n"

        get_reads(sequencefile, wanted_F_primers, file_contains_F_primer, n, paired_sequence_file, wanted_F_reads)
        get_reads(sequencefile, wanted_R_primers, file_contains_R_primer, n, paired_sequence_file, wanted_R_reads)
        n = n + 1

    # create a file to go into pal_finder which enables us to count the repeats in each read
    filePath = "pal_finder_input_file.fasta"
    unique_number = 0
    with open(filePath, 'w') as f:
        for row in wanted_F_reads:
            f.write(row + "\n")
    f.close()

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    print "Finished scanning samples for common primer sequences"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"

    wd = os.getcwd()
    if os.path.isdir(wd + "/Alignments"):
        shutil.rmtree(wd + "/Alignments")
    os.mkdir(wd + "/Alignments")
    output_path = wd + "/Alignments/"

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Creating config file for pal_finder"
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    print (time.strftime("%H:%M:%S"))
    # create a unique pal_finder config file for each
    # check for file of pal_finder configs, delete and re-create if necessary

    if os.path.isdir(wd + "/pal_finder_files"):
        shutil.rmtree(wd + "/pal_finder_files")
    os.mkdir(wd + "/pal_finder_files")
    output_path = wd + "/pal_finder_files/"
    assemble_reads("Forward_reads_for_assembly.fastq", "Reverse_reads_for_assembly.fastq", "Assembled_reads.fasta" )


    # section below needs 1 indent to the right
        # create new alignment for each set of primer sequences - made of the assembled reads
    # get a list of unique IDs contained in the file
    wanted = set()
    fasta1 = SeqIO.parse("Assembled_reads.fasta",'fasta')
    for seq in fasta1:
        sequence_ID = (seq.id.split(":")[7].split("_")[0])
        sequence = str(seq.seq)
        filePath = str(os.getcwd()) + '/Alignments/%s.fasta' % (sequence_ID)
        if sequence_ID in wanted:
            if os.path.exists(filePath):
                with open(filePath, 'a') as f:
                    f.write(">" + sequence_ID + "\n")
                    f.write(sequence + "\n")
                    #SeqIO.write([seq], f, "fasta")
        else:
            with open(filePath, 'w') as f:
                #SeqIO.write([seq], f, "fasta")
                f.write(">" + sequence_ID + "\n")
                f.write(sequence + "\n")
        wanted.add(sequence_ID)

    # even though some primer sequences are found in the data files enough times
    # to get through the process this far. it is possible that not every read
    # is sufficiently high quality to form an assembly.
    # this can result in the situation where alignments with only a single
    # sequence get through

    # check length of files in Alignments folder
    # remove anything with only one sequence
    for filename in os.listdir(wd + "/Alignments"):
        if not file_len(wd + "/Alignments/" + filename) > 2:
            os.remove(wd + "/Alignments/" + filename)


    """
    # get set of unique sequence IDs
    for record in wanted_F_reads:
        for x in record.split("\n"):
            wanted.add(record.split("\n")[0].lstrip(">"))
    # create a file for each and write out the sequences into each files
    for record in wanted_F_reads:
        #for x in record.split("\n"):
        if record.split("\n")[0].lstrip(">") in wanted:
            filePath = str(os.getcwd()) + '/Alignments/%s.fasta' % (record.split("\n")[0].lstrip(">"))
            if os.path.exists(filePath):
                with open(filePath, 'a') as f:
                    f.write(record + "\n")
            else:
                with open(filePath, 'w') as f:
                    f.write(record + "\n")
    """
    copy_and_rename_file(pal_finder_config, output_path + "pal_finder_config_file.txt")
    write_to_config(output_path + "pal_finder_config_file.txt","findPrimers 1", "findPrimers 0")
    write_to_config(output_path + "pal_finder_config_file.txt","platform Illumina", "platform 454")
    write_to_config(output_path + "pal_finder_config_file.txt","inputFormat fastq", "inputFormat fasta")
    write_to_config(output_path + "pal_finder_config_file.txt","pairedEnd  1","pairedEnd  0" )
    write_to_config(output_path + "pal_finder_config_file.txt","input454reads  test/data/454_All_python.fna", "input454reads " + wd + "/Assembled_reads.fasta")
    write_to_config(output_path + "pal_finder_config_file.txt","MicrosatSumOut  test/output/test_microsat_summary.txt", "MicrosatSumOut " + output_path + "pal_finder_summary_out.txt")
    write_to_config(output_path + "pal_finder_config_file.txt","PALsummaryOut  test/output/test_PAL_summary.txt", "PALsummaryOut " + output_path + "pal_finder_PAL_summary.txt")
    write_to_config(output_path + "pal_finder_config_file.txt","2merMinReps 	6", "2merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt","3merMinReps 	0", "3merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt","4merMinReps 	0", "4merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt","5merMinReps 	0", "5merMinReps 	6")
    write_to_config(output_path + "pal_finder_config_file.txt","6merMinReps 	0", "6merMinReps 	6")
    pal_finder_command = 'perl ' + pal_finder_script + " " + output_path + "pal_finder_config_file.txt"
    print "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    print "Running pal_finder."
    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
    print (time.strftime("%H:%M:%S"))
    subprocess.call(pal_finder_command, shell=True)
    wanted = set()
    fasta1 = SeqIO.parse("Assembled_reads.fasta",'fasta')
    for seq in fasta1:
        sequence_ID = (seq.id.split(":")[7].split("_")[0])
        filePath = str(os.getcwd()) + '/Alignments/%s.fasta' % (sequence_ID)
        if sequence_ID in wanted:
            if os.path.exists(filePath):
                with open(filePath, 'a') as f:
                    SeqIO.write([seq], f, "fasta")
            else:
                with open(filePath, 'w') as f:
                    SeqIO.write([seq], f, "fasta")
                wanted.add(sequence_ID)
        wd = os.getcwd()

    #### parse the pal_finder output to find the variable loci
    # get a list of unique IDs (the forward primer sequences)
    with open(wd + "/pal_finder_files/pal_finder_PAL_summary.txt") as pf:
        unique_primers = set()
        allele_count = []
        unique_alleles = []
        for line in pf:
            motifs = (line.split("\t")[3])
            if (len(motifs.split(" "))-1) == 1:
                ID = line.split("\t")[0]
                unique_primers.add(ID.split(":")[7])
                allele_count.append(0)
                unique_alleles.append("Motifs: ")
    list_unique_primers = list(unique_primers)

    ## go through and get the unique alleles associated with each primer sequence
    all_data = []
    with open(wd + "/pal_finder_files/pal_finder_PAL_summary.txt") as pf:
        for line in pf:
            motifs = (line.split("\t")[3])
            if (len(motifs.split(" "))-1) == 1:
                ID = line.split("\t")[0]
                for i, primer_line in enumerate(list_unique_primers):
                    if ID.split(":")[7] == primer_line:
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

    # calcualte the difference in size between the biggest number of repeats and the smallest
    diff_in_motif_size = []
    for row in single_motif_only:
        number_of_repeats = []
        for x in row.split("\t")[2].split(" ")[1:]:
            for result in (re.findall(r'\d+', x)):
                number_of_repeats.append(int(result))
        diff_in_motif_size.append((max(number_of_repeats) - min(number_of_repeats)))

    # rank by difference in largest and samllest motifs
    ranked_output = []
    count = 0
    if len(diff_in_motif_size) > 0:
        while count < max(diff_in_motif_size):
            for x, y in zip(single_motif_only, diff_in_motif_size):
                if y == count:
                    ranked_output.insert(0, x.replace("Motifs: ", " ") + "\t" + str(y))
            count = count + 1
        # write out final output
        with open("Mimi_output.txt", 'w') as final_output:
            final_output.write("Foward_primer_seq\tReverse_primer_seq\tNumber_of_alleles\tFound_in_individuals\tAlleles_present\tSize-Range\n")
            for x in ranked_output:
                final_output.write(x.split("\t")[0] + "\t" + ReverseComplement(primer_pairs[x.split("\t")[0]]) + "\t" + x.split("\t")[1] + "\t" + str(which_individuals[x.split("\t")[0]]) + "\t" + x.split("\t")[2] + "\t" + x.split("\t")[3] + "\n")
    else:
        print("Some error message here. Need to work out what exactly causes it to get here")

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
    print("\n\nSuccessfully Complete")
    print (time.strftime("%H:%M:%S"))
