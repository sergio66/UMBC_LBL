## these files gets copied to /home/sergio/SPECTRA/HITRAN_Database_important_scripts/2024

echo "Copy main dir"
/bin/cp -a *.m *.sc Readme* size_xsec_zipped_files_from_hitran CDC-NIOSH_Pocket_Guide_to_ChemicalHazards-CASNumbers.pdf readme-1.txt        /home/sergio/SPECTRA/HITRAN_Database_important_scripts/2024/.

########################################################################

echo "Copy subdir ISOTOPES"
DIRECTORY_NAME="/home/sergio/SPECTRA/HITRAN_Database_important_scripts/2024/ISOTOPES/"

# Check if the directory exists
if [ ! -d "$DIRECTORY_NAME" ]; then
    echo "Directory $DIRECTORY_NAME does not exist. Creating it now..."
    # Create the directory
    mkdir "$DIRECTORY_NAME"
    echo "Directory created."
else
    echo "Directory $DIRECTORY_NAME already exists."
fi

/bin/cp -a ISOTOPES/* /home/sergio/SPECTRA/HITRAN_Database_important_scripts/2024/ISOTOPES/.

########################################################################

echo "Copy subdir TIPS"
DIRECTORY_NAME="/home/sergio/SPECTRA/HITRAN_Database_important_scripts/2024/QTIPS/"

# Check if the directory exists
if [ ! -d "$DIRECTORY_NAME" ]; then
    echo "Directory $DIRECTORY_NAME does not exist. Creating it now..."
    # Create the directory
    mkdir "$DIRECTORY_NAME"
    echo "Directory created."
else
    echo "Directory $DIRECTORY_NAME already exists."
fi

/bin/cp -a QTIPS/* /home/sergio/SPECTRA/HITRAN_Database_important_scripts/2024/QTIPS/.
########################################################################

