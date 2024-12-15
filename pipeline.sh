#!/bin/bash

#GLOBALS
SOURCE_FOLDER=""
TARGET_FOLDER=""
REF_FOLDER=""
REF_NAME=""
REF_DATABASE=""
INPUT_ARRAY=()
NUMBER_OF_INPUT_FILES=0
PAIRED_DELIMITER=""
WORK_FILE_BASES=()
NUMBER_OF_WORKFILES=0
THREADS=""
DATA_FORMAT=""
KNOWN_SNPS=""
#visuals 
STEP_SEPERATOR_VISUAL="<<<------------------------------------------>>>"
checkmark="✔"
cross="✘"

#FUNCTION BLOCK
function unitInitReport (){
    local unit=$1
    local step=$2
    echo $STEP_SEPERATOR_VISUAL
    echo "Performing Step[ $unit / 7 ] -> $step "
    echo $STEP_SEPERATOR_VISUAL
}



function unitStatusReportSuccess () {
    local unit=$1
    local step=$2
    echo $STEP_SEPERATOR_VISUAL
    echo "Completed Step[ $unit / 7 ] -> $step succesfully performed $checkmark"
    echo $STEP_SEPERATOR_VISUAL
}

function unitStatusReportFailed () {
    local unit=$1
    local step=$2
    echo $STEP_SEPERATOR_VISUAL
    echo "Error at Step[ $unit / 7 ] -> $step failed $cross Exiting the pipeline now"
    echo $STEP_SEPERATOR_VISUAL
}




function globalConfiguration() {
    

# Create the folder if it doesn't exist
    

    # config 
    echo "Please provide a filename in json-Format for configuration path. If you want to enter the parameters manually, type 'manual'"
    read -p "ConfigFile " config

    
    # manual path configuration 
    if [ "$config" == "manual" ]; then
        read -p "Source Folder: " SOURCE_FOLDER
        read -p "Target Folder: " TARGET_FOLDER
        read -p "Referencegenome Folder: " REF_FOLDER
        read -p "Referencegenome Database: " REF_DATABASE
        read -p "Input dataformat ( gz, fastq, fasta...): " DATA_FORMAT
        read -p "Mode ( paired, single) " DATA_MODE
        read -p "Delimiter for forward and reverse reads (. , _ ) " PAIRED_DELIMITER
        read -p "Threads" THREADS
        read -p "Known Snps ( vcf file)" KNOWN_SNPS

    else
        SOURCE_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['sourceFolder'])")
        TARGET_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['targetFolder'])")
        REF_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['referenceGenome'])")
        REF_DATABASE=$(python3 -c "import json; print(json.load(open('$config'))['refDataBase'])")
        DATA_FORMAT=$(python3 -c "import json; print(json.load(open('$config'))['inputDataFormat'])")
        DATA_MODE=$(python3 -c "import json; print(json.load(open('$config'))['mode'])")
        PAIRED_DELIMITER=$(python3 -c "import json; print(json.load(open('$config'))['delimiter_pairedReads'])")
        THREADS=$(python3 -c "import json; print(json.load(open('$config'))['threads'])")
        KNOWN_SNPS=$(python3 -c "import json; print(json.load(open('$config'))['known_SNPS'])")

        
    fi

    if [ -d $SOURCE_FOLDER ] &&  [ -d $SOURCE_FOLDER ] && [ -e $REF_FOLDER ] ; then
        echo "Your choosen folders exist!"
        else
        echo "Attention! Folders do not exist!"
    fi
    mkdir -p "$TARGET_FOLDER"
    echo "Output Folder '$folder' has been created or already exists."

    # logging style -> NOT FUNCTIONAL YET
   
   
} 


function fileConfiguration() {

   # declaring input files 
  unitInitReport 1 "File configuration"
  for file in "$SOURCE_FOLDER"/*.$DATA_FORMAT; do
        #if [[ -f "$file" ]]; then
            INPUT_ARRAY+=($(basename "$file"))
            ln -s "$file" "$(basename "$file")"
            echo "Created symlink for $(basename "$file")"
            #echo "BASE" $(basename "$file") 
             
        #fi
    done
    NUMBER_OF_INPUTFILES=${#INPUT_ARRAY[@]}
 




    if [ $? -eq 0 ]; then
        unitStatusReportSuccess 1 "File configuration"
    else
        unitStatusReportFailed 1 "File configuration"
    fi



}


function printStartMessage () { 

    echo "_____________________________________________"
    echo "#####  ########    ###     #####    ########"
    echo "#         #       #   #    #    #       #  "
    echo "#####     #      #######   #####        # "
    echo "    #     #      #     #   #    #       #"      
    echo "#####     #      #     #   #     #      #"
    
    #summary of inital configuration
    echo "Starting pipeline with follwowing configutation: "
    echo "Source Folder         ->" $SOURCE_FOLDER
    echo "Target Folder         ->" $TARGET_FOLDER
    echo "Reference File        ->" $REF_FOLDER
    echo "Reference DB Genome   ->" $REF_DATABASE
    echo "Read Mode             ->" $DATA_MODE
    echo "Data Format           ->" $DATA_FORMAT
    echo "Input Files           -> "${INPUT_ARRAY[@]}""
    echo "Number of Input Files -> $NUMBER_OF_INPUTFILES"
    echo "Number of Threads     -> $THREADS"
    


}


function unit2_Alignment() {

  
  unitInitReport 2 "Alignment"
  if [ "$DATA_MODE" == "paired" ]; then
    for ((i=0; i<$NUMBER_OF_INPUTFILES; i+=2))
        do
            echo "${INPUT_ARRAY[i]} + ${INPUT_ARRAY[i+1]}"
            basename=$(echo "${INPUT_ARRAY[i]}" | awk -F $PAIRED_DELIMITER '{print $1}')
            outputBase="$TARGET_FOLDER$basename"
            WORK_FILE_BASES+=("$outputBase")
            bwa mem  -t $THREADS $REF_FOLDER "${INPUT_ARRAY[i]}" "${INPUT_ARRAY[i+1]}" > $outputBase.alignment.bam 
            samtools sort -@ $THREADS -o "$outputBase".alignment.sorted.bam $outputBase.alignment.bam &>>logfile.txt
            samtools index -@ $THREADS "$outputBase.alignment.sorted.bam" &>>logfile.txt
            samtools flagstat "$outputBase.alignment.sorted.bam" &>>logfile.txt
        done     
    elif [ "$DATA_MODE" == "single" ]; then
        for ((i=0; i<$NUMBER_OF_INPUTFILES; i+=1))
            do
                echo "${INPUT_ARRAY[i]}"
                basename=$(echo "${INPUT_ARRAY[i]}" | awk -F'.$DATA_FORMAT' '{print $1}')
                WORK_FILE_BASES+=("$OUTPUT_DIR$basename")
                bwa mem -t $THREADS $REF_FOLDER "${INPUT_ARRAY[i]}" > $outputFile &>>log.txt 
                samtools sort -o "$outputBase".alignment.sorted.bam $outputBase.alignment.bam &>>log.txt  
                samtools index "$outputBase.alignment.sorted.bam" &>>log.txt 
                samtools flagstat "$outputBase.alignment.sorted.bam" &>>log.txt 
                samtools tview -d C "$outputBase.alignment.sorted.bam" &>>log.txt 
            done    
    fi

    #number of files that will be generated for each step in the following workflow
    NUMBER_OF_WORKFILES="${#WORK_FILE_BASES[@]}"

    if [ $? -eq 0 ]; then
        unitStatusReportSuccess 2 "Alignment"
    else
        unitStatusReportFailed 2 "Alignment"
    fi
  }


function unit3_ReadgroupsAndDuplicateRemoval() {
     
    unitInitReport 3 "Readgroups and duplicare removal"
    for file in "${WORK_FILE_BASES[@]}"; do

            java -jar /group/bin/picard.jar MarkDuplicates \
                I="$file.alignment.sorted.bam"\
                O="$file.marked_duplicates.bam" \
                M="$file.duplication_metrics.txt" \
                REMOVE_DUPLICATES=true \
                &>>logfile.txt
          
            java -jar /group/bin/picard.jar AddOrReplaceReadGroups \
                I="$file.marked_duplicates.bam" \
                O="$file.rg_added.bam" \
                RGID=1 \
                RGLB=lib1 \
                RGPL=illumina \
                RGPU=unit1 \
                RGSM=sample1 \
                &>>logfile.txt
        done

    if [ $? -eq 0 ]; then
        unitStatusReportSuccess 3 "Readgroups and duplicare removal"
    else
        unitStatusReportFailed 3 "Readgroups and duplicare removal"
    fi


    }

function unit4_BaseRecalibration() {
    unitInitReport 4 "Base recalibration"
    for file in "${WORK_FILE_BASES[@]}";do 
        gatk BaseRecalibrator -I "$file.rg_added.bam" --known-sites $KNOWN_SNPS -O $file.table --reference $REF_FOLDER &>>logfile.txt
        gatk ApplyBQSR --input "$file.rg_added.bam" --output $file.recalibrated_reads.bam -bqsr $file.table &>>logfile.txt
    done 


     if [ $? -eq 0 ]; then
        unitStatusReportSuccess 4 "Base recalibration"
    else
        unitStatusReportFailed 4 "Base recalibration"
    fi


}

function unit5_VariantDetection() {
    unitInitReport 5 "Variant Detection"
    
    for file in "${WORK_FILE_BASES[@]}"; do
        gatk Mutect2 -I $file.recalibrated_reads.bam -O "$file.output_variants.vcf" -R $REF_FOLDER &>>logfile.txt
        gatk FilterMutectCalls -V "$file.output_variants.vcf" -O "$file.filtered.vcf" -R $REF_FOLDER &>>logfile.txt
    done


     if [ $? -eq 0 ]; then
        unitStatusReportSuccess 5 "Variant Detection"
    else
        unitStatusReportFailed 5 "Variant Detection"
    fi
}

function unit6_snpEffects() { 
  

    unitInitReport 6 "SNP Effect"
    for file in "${WORK_FILE_BASES[@]}"; do
        java -jar /group/opt/snpEff/snpEff.jar eff -htmlStats "$file.summary.html"  -csvStats "$file.summary.csv" GRCh38.86 "$file.filtered.vcf" > $file.annotated.vcf 
    done

     if [ $? -eq 0 ]; then
        unitStatusReportSuccess 6 "SNP Effect"
    else
        unitStatusReportFailed 6 "SNP Effect"
    fi
}


function unit7_GenerateHighImpactINDELSCSV() { 
    
    unitInitReport 7 "High Impact CSV"

    echo "Explorative file for high impact indels in the files "${WORK_FILE_BASES[@]}"" > $TARGET_FOLDER/high_impact.summary.csv 
    for file in "${WORK_FILE_BASES[@]}"; do
        echo "File: $file" >> $TARGET_FOLDER/high_impact.summary.csv 
        grep -w "HIGH" $file.annotated.vcf >> $TARGET_FOLDER/high_impact.summary.csv   
    done

    if [ $? -eq 0 ]; then
        unitStatusReportSuccess 7 "High Impact CSV"
    else
        unitStatusReportFailed 7 "High Impact CSV"
    fi
    

}



#BASIS KONFIGURATION (normal oder dev)
if [ "$1" == "dev" ]; then
    config="$2"
    SOURCE_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['sourceFolder'])")
    TARGET_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['targetFolder'])")
    REF_FOLDER=$(python3 -c "import json; print(json.load(open('$config'))['referenceGenome'])")
    REF_DATABASE=$(python3 -c "import json; print(json.load(open('$config'))['refDataBase'])")
    DATA_FORMAT=$(python3 -c "import json; print(json.load(open('$config'))['inputDataFormat'])")
    DATA_MODE=$(python3 -c "import json; print(json.load(open('$config'))['mode'])")
    PAIRED_DELIMITER=$(python3 -c "import json; print(json.load(open('$config'))['delimiter_pairedReads'])")
    THREADS=$(python3 -c "import json; print(json.load(open('$config'))['threads'])")
    KNOWN_SNPS=$(python3 -c "import json; print(json.load(open('$config'))['known_SNPS'])")
    mkdir -p "$TARGET_FOLDER"
    echo "Output Folder '$folder' has been created or already exists."


    fileConfiguration
    printStartMessage
    unit2_Alignment
    unit3_ReadgroupsAndDuplicateRemoval
    unit4_BaseRecalibration
    unit5_VariantDetection
    unit6_snpEffects
    unit7_GenerateHighImpactINDELSCSV


#MAIN PIPELINE
else
    runtime_paths_confirmed="n"
    runtime_files_confirmed="n"
    while [ $runtime_paths_confirmed != "y" ] ; do
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
        echo  "Input files" "${INPUT_ARRAY[@]}"
        read -p "Want to continue? (y/n) 'n' will allow you to adjust the settings" runtime_files_confirmed
    done
fi


printStartMessage
startPipeline="n"
read -p "Start pipeline?(y/n) 'n' will terminate the process"  startPipeline

if [ $startPipeline == "y" ];then

    unit2_Alignment
    unit3_ReadgroupsAndDuplicateRemoval
    unit4_BaseRecalibration
    unit5_VariantDetection
    unit6_snpEffects
    unit7_GenerateHighImpactINDELSCSV

else
 echo "Abort pipeline start. Goodbye"
 exit 0
fi





