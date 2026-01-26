#!/bin/bash

################# This script takes tab-/comma-delimited file with sample information, and executes GATK-SV pipeline
################# Tab-/comma-delimited file should have the following columns: Sample_ID, Path/to/normal.bam, Path/to/normal.bam.bai

## REQUIRED INPUTS (Change these paths as needed)
export MANIFEST="/projects/standard/venteicher_30050/balay011/gatk-sv/chordbai_gatksv.tsv" # <- adjust
export OUTDIR="/scratch.global/balay011/GermlineSV_outs" # <- adjust
export GATK_SV_DIR="/scratch.global/balay011/GermGatkSVonHPC" # <- adjust
export CROMWELL_CONF="$GATK_SV_DIR/cromwell.conf"
export CROMWELL_EXE_JAR="/users/1/balay011/helper_scripts/gatk_sv_cromwell_scripts/cromwell-91.jar" # <- adjust
export DEPS_ZIP="$GATK_SV_DIR/deps.zip" # <- adjust

# Check required variables
if [[ -z "$MANIFEST" || -z "$OUTDIR" || -z "$GATK_SV_DIR" ]]; then
    echo "Error: Missing required environment variables."
    exit 1
fi

# ---------------- GatherSampleEvidence ----------------
## Module parameters
export JSON_GATHERSAMPLEEVIDENCE_TEMPLATE="$GATK_SV_DIR/data/GatherSampleEvidence_inputs.json"

### Output directory
mkdir -p ${OUTDIR}/GatherSampleEvidence_out/{results,logs}

while IFS=$'\t' read -r SAMPLE_ID BAM_PATH BAI_PATH; do

    # Check if sample ID is valid
    [[ -z "$SAMPLE_ID" || "$SAMPLE_ID" == \#* ]] && continue

    # Check if BAM and index files exist
    if [[ ! -f "$BAM_PATH" || ! -f "$BAI_PATH" ]]; then
        echo "Error: Missing BAM/BAI for $SAMPLE_ID. Skipping..."
        continue
    fi

    # Sample-level output directory
    GSE_RESULTS_DIR="${OUTDIR}/GatherSampleEvidence_out/results"
    SAMPLE_WORK_DIR="$GSE_RESULTS_DIR/$SAMPLE_ID"
    mkdir -p "$SAMPLE_WORK_DIR"

    # Generate Sample-Specific JSON from Template
    sed -e "s|PLACEHOLDER_ID|$SAMPLE_ID|g" \
        -e "s|PLACEHOLDER_BAM|$BAM_PATH|g" \
        -e "s|PLACEHOLDER_BAI|$BAI_PATH|g" \
        "$JSON_GATHERSAMPLEEVIDENCE_TEMPLATE" > "$SAMPLE_WORK_DIR/${SAMPLE_ID}_inputs.json"

    # Submit the Master Cromwell Job to Slurm
    sbatch \
      --job-name="GSE_master_$SAMPLE_ID" \
      --output="${OUTDIR}/GatherSampleEvidence_out/logs/${SAMPLE_ID}_gse_master.out" \
      --error="${OUTDIR}/GatherSampleEvidence_out/logs/${SAMPLE_ID}_gse_master.err" \
      --account=aventeic \
      --partition=asvnode1,msibigmem \
      --time=48:00:00 \
      --cpus-per-task=4 \
      --mem=16G \
      --mail-type=END,FAIL \
      --mail-user=balay011@umn.edu \
      --wrap "module load java/openjdk-17.0.2 && cd $SAMPLE_WORK_DIR && java -Xmx12G -Dconfig.file=$CROMWELL_CONF -jar $CROMWELL_EXE_JAR run $GATHERSAMPLEEVIDENCE_WDL_PATH -i ${SAMPLE_ID}_inputs.json -p $DEPS_ZIP"

    echo "Submitted pipeline for sample: $SAMPLE_ID"

done < "$MANIFEST"