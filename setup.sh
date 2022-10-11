#!/bin/bash
#

### Doccumentation string for help.
if [ "$1" = "-h"  -o "$1" = "--help" ]     # Request help.
then
  echo; printf "This bash script initialises a target directory in which to record the raw experimental data.\nThe initialised directory should contain subdirectories for respective cameras and their world coordinates."; echo
  echo; printf "Usage: \n$ $0 run_number=<Run Number> set_date=<Folder Date>\n"; echo
  sed --silent -e '/DOCUMENTATIONXX$/,/^DOCUMENTATIONXX$/p' "$0" |
  sed -e '/DOCUMENTATIONXX$/d'; exit $DOC_REQUEST; fi


: <<DOCUMENTATIONXX
run_number: The run number of the experiment to be recorded. (Defaults to be the first run)
set_date: The date set to be used when creating the target directory. (Defualts to the execution date)
------------------------------------------------------------------------------------------------------
Please note: In the parent directory there should exist a path to 
"Project/TwoCamerasSuperTestUpdated210921.dfc" which is to be used when initialising both cameras. 
Additionally there should be two subdirectories "CamA/" and "CamB/" which each contain the following 
files:
CoordinateSystems.log
DigiFlow_Dialogs.dfs
DigiFlow_Status.dfs
wcs_first_cut.dfi
The above files should contain all relevent information to encapsulate the world co-ordinates for each
respective cameras.
DOCUMENTATIONXX

### We need to read in all values defined on the command line.
#loop over all arguments given on command line.
for ARGUMENT in "$@"
do

    #Identify input arguments and their values.
    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=)   

    #Match expected input arguments to that given on the command line.
    case "$KEY" in
            run_number)              run_number=${VALUE} ;;
    	    set_date)                set_date=${VALUE} ;;     
            *)   
    esac    
done

### Check if the user has entered predefined 'run_number' or 'set_date'
if [ -z ${run_number+x} ]; then printf "\n'run_number' is undefined by user."; run_number=1; fi
printf "\nrun_number is set to ${run_number}. \n"

if [ -z ${set_date+x} ]; then printf "\n'set_date' is undefined by user."; set_date=$(date +"%d%m%y"); setdate=$(date +"%d%m%Y"); fi
printf "\nset_date is set to '$set_date'. \n"

### Directory Structure inputted as fixed parameters
D_CamDir=$(pwd)
D_ProbeDir=$(pwd)/../Probe_Data
OD_CamDir="/c/Users/b5006861/OneDrive - Newcastle University/02_PhD_Project/04_Ice_Covered_Waters/02_Raw_data/01_CameraData/"
OD_ScriptsDir="/c/Users/b5006861/OneDrive - Newcastle University/02_PhD_Project/04_Ice_Covered_Waters/03_Code/01_LabPlotting/"

### Collate filenames
### Define target directories based on 'run_number' and 'set_date'.
if [[ $run_number -eq 1 ]]
then
	orig_cam_dir=$D_CamDir/$set_date
	orig_probe_dir=$D_ProbeDir/$setdate
	new_cam_dir=$OD_CamDir/$set_date
	D_CamSettings=$(pwd)/Settings

else
	orig_cam_dir=$D_CamDir/${run_number}_/$set_date
	orig_probe_dir=$D_ProbeDir/${run_number}_/$setdate
	new_cam_dir=$OD_CamDir/${run_number}_/$set_date
	D_CamSettings=$(pwd)/${run_number}_/Settings

fi
#CamA_dir=$orig_cam_dir/Setting/
#CamB_dir=$orig_cam_dir/CamB/


### If the target directory is pre existing the code is aborted
if [[ -d "$orig_cam_dir" ]]
then 
	printf "\n***ABORTED*** \nTarget directory already exists, please check you are not trying to override an existing directory. \n${orig_cam_dir}"
	exit 1
fi

### Create the target directory.
printf "\nCreating Directory: \n$orig_cam_dir"
mkdir -p $orig_cam_dir; 


### Copy premade Camera content.
printf "\nCopying premade camera folder content."
cp ${D_CamSettings}/* ${orig_cam_dir};

### Copy pre-made Probe Content
hostnm="$(hostname)"
if [[ ${hostnm} == "18S-STB-48039" ]]
then  
printf "\nCreating Directory: \n$orig_probe_dir"
mkdir -p $orig_probe_dir; 
printf "\n Copying premade Probe Folder content"
cp ${D_ProbeDir}/ProbeParse.m $orig_probe_dir

printf "\nMaking OneDrive Copy: \n$new_cam_dir"
mkdir -p "${new_cam_dir}/CamA"; mkdir -p "${new_cam_dir}/CamB"; mkdir -p "${new_cam_dir}/CamC"; 
printf "\nGetting RunParticle Displacmeents Script:"
cp "${OD_ScriptsDir}/run_particledisplacements.m" "$new_cam_dir"

fi


printf "\nComplete!"
