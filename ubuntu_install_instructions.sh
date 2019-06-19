## MiMi instructions on a fresh Ubuntu 19.04 installation
## this document could also be run as a bash script, and should install all
## software correctly
printf "\n\n\nMiMi - Multi-individual Microsatellite Identification - installer script\n\n\n"
sleep 1
printf "This script will attempt to install everything you need to run MiMi on a Ubuntu 19.04 installation.\n\n\n"
sleep 3
########################
## PANDASEQ INSTALLATION
########################

## install PANDAseq dependencies
printf "Installing PANDAseq dependencies....\n\n\n"
sleep 3
sudo apt-get install build-essential libtool automake zlib1g-dev libbz2-dev pkg-config

print "\n\n\nInstalling git (if required)....\n\n\n"
sleep 3
## install git
sudo apt install git

## clone PANDAseq repository
printf "Downloading PANDAseq....\n\n\n"
sleep 3
git clone https://github.com/neufeld/pandaseq.git

## move into PANDAseq direction
cd pandaseq/

## build the PANDAseq package

printf "Installing PANDAseq....\n\n\n"
sleep 3
/autogen.sh && ./configure && make && sudo make install

## then run

sudo ldconfig

## test PANDAseq
## the command:

printf "Testing Muscle installation.......\n\n\n"
sleep 5
printf "You are hopefully about to see some output which begins: You must supply both forward and reverse reads.\n\n\n"
sleep 5
pandaseq

## should produce the following output:

# You must supply both forward and reverse reads.
# Too confused to continue.
# Try -h for help.

#########################
## INSTALL PYTHON MODULES
#########################

## install pip
printf "\n\n\nInstalling pip.....\n\n\n"
sleep 5
sudo apt install python-pip

## install biopython

printf "\n\n\nInstalling biopython.....\n\n\n"
sleep 5
pip install biopython

#########################
## INSTALL MUSCLE ALIGNER
#########################

printf "\n\n\nInstalling the MUSCLE aligner.....\n\n\n"
sleep 5
## install muscle
sudo apt install muscle

## test muscle
printf "Testing Muscle installation.......\n\n\n"
sleep 5
printf "You should hopefully see some output which begins:\n MUSCLE v3.8.1551 by Robert C. Edgar\n\n\n"
sleep 3
muscle

printf "\n\n\nIf no error messages occurred, and you saw the expected output during tests, then MiMi should be ready to run.\n\n"
printf "Test MiMi with the following command: (check the version number and amend if necessary)\n\n"
printf "./MiMi_v0.04.py -c demo_data/demo_config.txt\n\n\n"
