#!/bin/bash

#GLOBALS
SOURCE_FOLDER=""
TARGET_FOLDER=""
REF_FOLDER=""
REF_NAME=""
REF_DATABASE=""

INPUT_ARRAY=()
NUMBER_OF_INPUT_FILES=0


WORK_FILE_BASES=()
NUMBER_OF_WORKFILES=0



READS_TYPE="default"
OUTPUT_DIR="./output"
is_log=false
logfileBase="LogfileID"
logfile=""



function globalConfiguration() {
    echo $(date)
    echo "Welcome to aligment pipline Gruppe3Pipeline"
    echo "Select an Option:"
    echo "Please provide a filename in json-Format for configuration path. If you want to enter the parameters manually, type 'manual'"
    read -p "ConfigFile " config

    

    if [ "$config" == "manual" ]; then
        read -p "SourceFolder: " SOURCE_FOLDER
        read -p "TargetFolder: " TARGET_FOLDER
        read -p "ReferencegenomeFolder: " REF_FOLDER
    else
        SOURCE_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['sourceFolder'])")
        TARGET_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['targetFolder'])")
        REF_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['referenceGenome'])")
        REF_DATABASE=$(python3 -c "import json; print(json.load(open('$config'))['refDataBase'])")
        
    fi

    if [ -d $SOURCE_FOLDER ] &&  [ -d $SOURCE_FOLDER ] && [ -e $REF_FOLDER ] ; then
        echo "File exists!"
        else
        echo "File does not exist!"
    fi



    read -p "You have want to generate a logging file? y/n " log
    if [ "$log" == "y" ]; then
        is_log=true
        timestamp=$(date +%s)
        logfile="$logfileBase$timestamp"
        echo "Logfile will be created with the name: $logfile"
        echo "Variant Calling Pipeline by Group 3. Logging mode True."
        echo $(date) 
        echo "Unique Timestamp" 
    fi
   
} 


function fileConfiguration() {
  echo "File configuration"
    for file in "$SOURCE_FOLDER"/*.gz; do
        if [ -f "$file" ]; then
            ln -s "$file" "$(basename "$file")"
            echo "Created symlink for $(basename "$file")" >>log.txt 
            echo "BASE" $(basename "$file") >>log.txt 
            INPUT_ARRAY+=($(basename "$file")) 
        fi
    done
    NUMBER_OF_INPUTFILES=${#INPUT_ARRAY[@]}
    echo "INPUT:" "$INPUT_ARRAY"
    echo "NUMBER:" "$NUMBER_OF_INPUTFILES"

    echo "You will either work with single-or paired-end data. Enter read Mode:(single/paired)" 
    read num
    case $num in
            single)
                echo "You selected option single"
                READS_TYPE=$num
                reads_mode="y"
                ;;
            paired)
                echo "You selected option paird"
                READS_TYPE=$num
                reads_mode="y"
                ;;
        
            *)
                echo "Invalid selection"
                ;;
        esac


}


function printStartMessage () { 

    echo "_____________________________________________"
    echo "#####  ########    ###     #####    ########"
    echo "#         #       #   #    #    #       #  "
    echo "#####     #      #######   #####        # "
    echo "    #     #      #     #   #    #       #"      
    echo "#####     #      #     #   #     #      #"
   
    echo "Starting pipeline with follwowing configutation: "
    echo "Source Folder         ->" $SOURCE_FOLDER
    echo "Target Folder         ->" $TARGET_FOLDER
    echo "Reference File        ->" $REF_FOLDER
    echo "Reference DB Genome   ->" $REF_DATABASE
    echo "Read Mode             ->" $READS_TYPE
    echo "Input Files           ->  "${INPUT_ARRAY[@]}""
    echo "Number of Input Files -> $NUMBER_OF_INPUTFILES"
    echo "_____________________________________________"


}


function unit2_Alignment() {

  
  echo "Starting Unit 2...Performing Alignment via Sammtools"
  if [ "$READS_TYPE" == "paired" ]; then
    for ((i=0; i<$NUMBER_OF_INPUTFILES; i+=2))
        do
            echo "${INPUT_ARRAY[i]} + ${INPUT_ARRAY[i+1]}"
            basename=$(echo "${INPUT_ARRAY[i]}" | awk -F'.R1' '{print $1}')
            outputBase="$TARGET_FOLDER$basename"
            WORK_FILE_BASES+=("$outputBase")

            #outputFile="$OUTPUT_DIR$basename".alignment.bam 

            bwa mem  $REF_FOLDER "${INPUT_ARRAY[i]}" "${INPUT_ARRAY[i+1]}" > $outputBase.alignment.bam 
            samtools sort -o "$outputBase".alignment.sorted.bam $outputBase.alignment.bam 
            samtools index "$outputBase.alignment.sorted.bam" 
            samtools flagstat "$outputBase.alignment.sorted.bam" 
            #samtools tview -d C "$outputBase.alignment.sorted.bam" >>log.txt 
            #ALIGNMENT_FILES_ARRAY+=("$outputBase.alignment.sorted.bam")

        done     
    elif [ "$READS_TYPE" == "single" ]; then
        for ((i=0; i<$NUMBER_OF_INPUTFILES; i+=1))
            do
                echo "${INPUT_ARRAY[i]}"
                basename=$(echo "${INPUT_ARRAY[i]}" | awk -F'.fastq' '{print $1}')
                WORK_FILE_BASES+=("$OUTPUT_DIR$basename")
                #outputFile="$OUTPUT_DIR$basename".alignment.bam 
                bwa mem  $REF_FOLDER "${INPUT_ARRAY[i]}" > $outputFile >>log.txt 
                samtools sort -o "$outputBase".alignment.sorted.bam $outputBase.alignment.bam >>log.txt  
                samtools index "$outputBase.alignment.sorted.bam" >>log.txt 
                samtools flagstat "$outputBase.alignment.sorted.bam" >>log.txt 
                samtools tview -d C "$outputBase.alignment.sorted.bam" >>log.txt 
                #ALIGNMENT_FILES_ARRAY+=("$outputBase.alignment.sorted.bam")
                #ALIGNMENT_FILES_ARRAY+=("$OUTPUT_DIR$basename".alignment.sorted.bam)

            done    

    fi

    NUMBER_OF_WORKFILES="${#WORK_FILE_BASES[@]}" 
    #NUMBER_OF_ALIGNMENTFILES=${#ALIGNMENT_FILES_ARRAY[@]}
 

  }


function unit3_ReadgroupsAndDuplicateRemoval() {
     
    echo "Starting Unit 3...Performing Removal of Duplicats and adding of Readgroups"
    
    for file in "${WORK_FILE_BASES[@]}"; do

            

            java -jar /group/bin/picard.jar MarkDuplicates \
                I="$file.alignment.sorted.bam"\
                O="$file.marked_duplicates.bam" \
                M="$file.duplication_metrics.txt" \
                REMOVE_DUPLICATES=true
            echo "#######################################"
            java -jar /group/bin/picard.jar AddOrReplaceReadGroups \
                I="$file.marked_duplicates.bam" \
                O="$file.rg_added.bam" \
                RGID=1 \
                RGLB=lib1 \
                RGPL=illumina \
                RGPU=unit1 \
                RGSM=sample1
        #DUPLICTATE_REMOVED_FILES_ARRAY+=("$file.rg_added.bam")
        #NUMBER_OF_REMOVED_DUPLICATES_FILES=${#DUPLICTATE_REMOVED_FILES_ARRAY[@]}       
        done

    }

function unit4_BaseRecalibration() {
    for ((i=0; i<$NUMBER_OF_REMOVED_DUPLICATES_FILES; i+=1))
    do
        echo "Starting Unit 4...Performing Base Calibration and applying BaseQualityScoreRecalibration"
        basename=$(echo "${DUPLICTATE_REMOVED_FILES_ARRAY[i]}" | awk -F'.rg_added' '{print $1}')

        snp="/group/lectures/variantcalling/bioinfo23_data/reference/GATK_bundle_20220203/hg38_v0/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        gatk BaseRecalibrator -I "${DUPLICTATE_REMOVED_FILES_ARRAY[i]}" --known-sites $snp -O $basename.table --reference $REF_FOLDER
        gatk ApplyBQSR --input "${DUPLICTATE_REMOVED_FILES_ARRAY[i]}" --output $basename.recalibrated_reads.bam -bqsr $basename.table 
        CALIBRATED_FILES+=($basename.recalibrated_reads.bam)

    done 

    NUMBER_OF_CALIBRATED_FILES=${#CALIBRATED_FILES[@]}


}

function unit5_VariantDetection() {
    echo "Starting Unit 5...Performing Variant Detection and Filtering"
    #gatk Mutect2 -I $OUTPUT_DIR/recalibrated_reads.bam -O $OUTPUT_DIR/output_variants.vcf -R $REF_FOLDER
    gatk FilterMutectCalls -V $OUTPUT_DIR/output_variants.vcf -R  $REF_FOLDER -O $OUTPUT_DIR/filtered.vcf


}

function unit6_snpEffects() { 
    echo "Starting Unit 6...Performing SnpEffects"
    java -jar /group/opt/snpEff/snpEff.jar eff -csvStats summary.csv GRCh38.86 $OUTPUT_DIR/filtered.vcf >annotated.vcf
    
}




#BASIS KONFIGURATION (normal oder dev)
if [ "$1" == "dev" ]; then
    config="config.json"
    SOURCE_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['sourceFolder'])")
    TARGET_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['targetFolder'])")
    REF_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['referenceGenome'])")
    fileConfiguration
    printStartMessage
    unit2_Alignment
    echo "${WORK_FILE_BASES[@]}"
    #ALIGNMENT_FILES_ARRAY=("output/tumor.alignment.sorted.bam" "output/control.alignment.sorted.bam")
    #NUMBER_OF_ALIGNMENTFILES=2

    unit3_ReadgroupsAndDuplicateRemoval
    #unit4_BaseRecalibration
else
    runtime_paths_confirmed="n"
    runtime_files_confirmed="n"
    while [ $runtime_paths_confirmed != "y" ] ; do
        echo "Looping because var is false."
        # Set var to true or exit the loop after some condition to stop it
        globalConfiguration  
        echo "______________________________________"
        echo "Source Folder " $SOURCE_FOLDER
        echo "Target Folder " $TARGET_FOLDER
        echo "Reference File " $REF_FOLDER
        echo "Reference DB Genome "$REF_DATABASE
        read -p "Want to continue?(y/n)n will allow you to adjust the settings" runtime_paths_confirmed
    done

    while [ $runtime_files_confirmed != "y" ] ; do
        fileConfiguration
        echo "Following files have been found and will be analysed:"
        echo  "Number of files" "${INPUT_ARRAY[@]}"
        read -p "Want to continue? (y/n) 'n' will allow you to adjust the settings" runtime_files_confirmed
    done
fi



printStartMessage

#unit2_Alignment
#unit3_ReadgroupsAndDuplicateRemoval
#unit4_BaseRecalibration
#unit5_VariantDetection
#unit6_snpEffects





