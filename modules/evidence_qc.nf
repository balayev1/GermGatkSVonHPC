#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process EVIDENCE_QC {

    input:
    path "samples/*"

    output:
    path "results/*", emit: qc_results

    script:

    def template_path = file(params.evidqc_template ?: "${params.gatk_sv_dir}/data/EvidenceQC_inputs.json").toAbsolutePath()
    
    """
    # Use Python to load the template and inject dynamic file lists
    python3 - <<EOF
    import json
    import os
    import glob

    # Load the static JSON template
    with open("${template_path}", 'r') as f:
        data = json.load(f)

    samples = []
    counts = []
    manta = []
    wham = []
    scramble = []

    sample_dirs = sorted([d for d in os.listdir('samples') if os.path.isdir(os.path.join('samples', d))])

    for sample_id in sample_dirs:
        sample_path = os.path.join('samples', sample_id)
        
        def find_file(pattern, sub_string):
            found = glob.glob(f"{sample_path}/**/{pattern}", recursive=True)
            filtered = [os.path.abspath(f) for f in found if sub_string in f]
            return filtered[0] if filtered else None

        c_file = find_file(f"{sample_id}.counts.tsv.gz", "call-CountsMetrics")
        m_file = find_file(f"{sample_id}.manta.std.vcf.gz", "execution")
        w_file = find_file(f"{sample_id}.wham.std.vcf.gz", "execution")
        s_file = find_file(f"{sample_id}.scramble.vcf.gz", "execution")

        if c_file:
            samples.append(sample_id)
            counts.append(c_file)
            manta.append(m_file)
            wham.append(w_file)
            scramble.append(s_file)

    # Overwrite the dynamic arrays while keeping other static keys
    data["EvidenceQC.samples"] = samples
    data["EvidenceQC.counts"] = counts
    data["EvidenceQC.manta_vcfs"] = manta
    data["EvidenceQC.wham_vcfs"] = wham
    data["EvidenceQC.scramble_vcfs"] = scramble

    # Write the final merged JSON
    with open('evidence_qc_inputs.json', 'w') as f:
        json.dump(data, f, indent=4)
    EOF

    # Run Evidence_QC using Cromwell
    java -Xmx48G -Dconfig.file=${params.cromwell_conf} -jar ${params.cromwell_jar} \
        run ${params.evidqc_wdl} \
        -i evidence_qc_inputs.json \
        -p ${params.deps_zip}

    # Move files to output directory
    mkdir -p results
    mv cromwell-executions/EvidenceQC/*/call-* results/
    """
}