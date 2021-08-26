## MiMi instructions on a fresh Ubuntu 19.04 installation
## this document could also be run as a bash script, and should install all
## software correctly
printf "\n\n\nMiMi - Multi-individual Microsatellite Identification - installer script\n\n\n"
sleep 1
printf "This script will attempt to install everything you need to run MiMi on a Ubuntu 21.04 installation.\n\n\n"
sleep 3
########################
## PANDASEQ INSTALLATION
########################

## install PANDAseq dependencies
printf "\n\n\nInstalling PANDAseq dependencies....\n\n\n"
sudo apt-get install build-essential libtool automake zlib1g-dev libbz2-dev pkg-config -y

printf "\n\n\nInstalling git (if required)....\n\n\n"
## install git
sudo apt install git

## clone PANDAseq repository
printf "\n\n\nDownloading PANDAseq....\n\n\n"
git clone https://github.com/neufeld/pandaseq.git

## move into PANDAseq direction
cd pandaseq/

## build the PANDAseq package

printf "\n\n\nInstalling PANDAseq....\n\n\n"
./autogen.sh && ./configure && make && sudo make install

## then run
sudo ldconfig
## install pip
printf "\n\n\nInstalling pip.....\n\n\n"
sudo apt install python3-pip -y

## install biopython
printf "\n\n\nInstalling biopython.....\n\n\n"
pip3 install biopython

#########################
## INSTALL MUSCLE ALIGNER
#########################

printf "\n\n\nInstalling the MUSCLE aligner.....\n\n\n"
## install muscle
sudo apt install muscle -y

printf "\n\n\nTest MiMi with the following command: (check the version number and amend if necessary)\n\n"
printf "./MiMi_v0.03.py -c demo_data/demo_config.txt\n\n\n"
