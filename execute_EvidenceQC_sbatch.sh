#!/bin/bash

################# This script takes GatherSampleEvidence outputs from execute_GatherSampleEvidence_sbatch.sh, and executes EvidenceQC module in GATK-SV pipeline
################# GatherSampleEvidence outputs MUST be located in the following directory structure: OUTDIR/GatherSampleEvidence_out/results/SAMPLE_ID/
################# and contain the following files: 
################# Read Counts TSV: OUTDIR/GatherSampleEvidence_out/results/SAMPLE_ID/*/call-CountsMetrics/inputs/*/SAMPLE_ID.counts.tsv.gz
################# Standardized Manta VCF (optional): OUTDIR/GatherSampleEvidence_out/results/SAMPLE_ID/*/call-Manta_Std/execution/SAMPLE_ID.manta.std.vcf.gz
################# Standardized Wham VCF (optional): OUTDIR/GatherSampleEvidence_out/results/SAMPLE_ID/*/call-Wham_Std/execution/SAMPLE_ID.wham.std.vcf.gz
################# Scramble VCF (optional): OUTDIR/GatherSampleEvidence_out/results/SAMPLE_ID/*/call-MakeScrambleVcf/execution/SAMPLE_ID.scramble.vcf.gz

## REQUIRED INPUTS (Change these paths as needed)
export OUTDIR="/scratch.global/balay011/GermlineSV_outs" # <- adjust
export GATK_SV_DIR="/scratch.global/balay011/GermGatkSVonHPC" # <- adjust
export CROMWELL_CONF="$GATK_SV_DIR/cromwell.conf"
export CROMWELL_EXE_JAR="/users/1/balay011/helper_scripts/gatk_sv_cromwell_scripts/cromwell-91.jar" # <- adjust
export DEPS_ZIP="$GATK_SV_DIR/deps.zip" # <- adjust
export EVIDENCEQC_WDL_PATH="/projects/standard/venteicher_30050/balay011/gatk-sv/wdl/EvidenceQC.wdl" # <- adjust

# Check required variables
if [[ -z "$OUTDIR" || -z "$GATK_SV_DIR" ]]; then
    echo "Error: Missing required environment variables."
    exit 1
fi

# ---------------- EvidenceQC ----------------
## Module parameters
export JSON_EVIDENCEQC_TEMPLATE="$GATK_SV_DIR/data/EvidenceQC_inputs.json"
export GATHERSAMPLEEVIDENCE_OUTDIR="$OUTDIR/GatherSampleEvidence_out/results"
export EVIDENCEQC_OUTDIR="$OUTDIR/EvidenceQC_out"

SAMPLES=""
COUNTS=""
MANTA=""
WHAM=""
SCRAMBLE=""
echo "Gathering paths from $GATHERSAMPLEEVIDENCE_OUTDIR..."
for SAMPLE_PATH in $(ls -d ${GATHERSAMPLEEVIDENCE_OUTDIR}/*/ | sort); do
    SAMPLE_ID=$(basename "$SAMPLE_PATH")
    
    # Locate files
    C_FILE=$(find "$SAMPLE_PATH" -name "${SAMPLE_ID}.counts.tsv.gz" | grep "call-CountsMetrics" | head -n 1)
    M_FILE=$(find "$SAMPLE_PATH" -name "${SAMPLE_ID}.manta.std.vcf.gz" | grep "execution" | head -n 1)
    W_FILE=$(find "$SAMPLE_PATH" -name "${SAMPLE_ID}.wham.std.vcf.gz" | grep "execution" | head -n 1)
    S_FILE=$(find "$SAMPLE_PATH" -name "${SAMPLE_ID}.scramble.vcf.gz" | grep "execution" | head -n 1)

    if [[ -n "$C_FILE" ]]; then
        SAMPLES+="\"$SAMPLE_ID\", "
        COUNTS+="\"$C_FILE\", "
        MANTA+="\"$M_FILE\", "
        WHAM+="\"$W_FILE\", "
        SCRAMBLE+="\"$S_FILE\", "
    fi
done

SAMPLES=${SAMPLES%, }
COUNTS=${COUNTS%, }
MANTA=${MANTA%, }
WHAM=${WHAM%, }
SCRAMBLE=${SCRAMBLE%, }

## Replace placeholders in JSON template by adding arrays for samples, counts and if found, Manta, Wham and Scramble VCFs
echo "Setting up EvidenceQC JSON input file..."
sed -i "/\"SAMPLE_ID1\"/,/\"SAMPLE_ID2\"/c\  $SAMPLES" "$JSON_EVIDENCEQC_TEMPLATE"
sed -i "/\"\/path\/to\/results\/SAMPLE_ID1\/SAMPLE_ID1.counts.tsv.gz\"/,/\"\/path\/to\/results\/SAMPLE_ID2\/SAMPLE_ID2.counts.tsv.gz\"/c\  $COUNTS" "$JSON_EVIDENCEQC_TEMPLATE"
sed -i "/\"\/path\/to\/results\/SAMPLE_ID1\/SAMPLE_ID1.manta.std.vcf.gz\"/,/\"\/path\/to\/results\/SAMPLE_ID2\/SAMPLE_ID2.manta.std.vcf.gz\"/c\ $MANTA" "$JSON_EVIDENCEQC_TEMPLATE"
sed -i "/\"\/path\/to\/results\/SAMPLE_ID1\/SAMPLE_ID1.wham.std.vcf.gz\"/,/\"\/path\/to\/results\/SAMPLE_ID2\/SAMPLE_ID2.wham.std.vcf.gz\"/c\ $WHAM" "$JSON_EVIDENCEQC_TEMPLATE"
sed -i "/\"\/path\/to\/results\/SAMPLE_ID1\/SAMPLE_ID1.scramble.vcf.gz\"/,/\"\/path\/to\/results\/SAMPLE_ID2\/SAMPLE_ID2.scramble.vcf.gz\"/c\ $SCRAMBLE" "$JSON_EVIDENCEQC_TEMPLATE"
echo "Setup complete: $JSON_EVIDENCEQC_TEMPLATE"

## Verify JSON file
python3 -c "import json; json.load(open('$JSON_EVIDENCEQC_TEMPLATE'))" && echo "JSON is valid" || exit 1

## Create output directories
mkdir -p "$EVIDENCEQC_OUTDIR/results"
mkdir -p "$EVIDENCEQC_OUTDIR/logs"

sbatch \
  --job-name="EvidenceQC_Cohort" \
  --output="$EVIDENCEQC_OUTDIR/logs/evidence_qc.out" \
  --error="$EVIDENCEQC_OUTDIR/logs/evidence_qc.err" \
  --account=aventeic \
  --partition=asvnode1,msibigmem \
  --time=96:00:00 \
  --cpus-per-task=8 \
  --mem=64G \
  --mail-type=END,FAIL \
  --mail-user=balay011@umn.edu \
  --wrap "module load java/openjdk-17.0.2 && cd $EVIDENCEQC_OUTDIR/results && java -Xmx48G -Dconfig.file=$CROMWELL_CONF -jar $CROMWELL_EXE_JAR run $EVIDENCEQC_WDL_PATH -i $JSON_EVIDENCEQC_TEMPLATE -p $DEPS_ZIP"

echo "Cohort EvidenceQC submitted."