#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process EVIDENCE_QC {

    input:
    val sample_list
    path publish_dir_base

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

    for sample_id in ${sample_list.inspect()}:
        # Construct the path to the published directory
        sample_path = os.path.join("${publish_dir_base}", sample_id)
        
        if not os.path.exists(sample_path):
            print(f"WARNING: Published directory for {sample_id} not found at {sample_path}")
            continue

        def find_file(pattern, sub_string):
            # Search in the published results folder
            found = glob.glob(f"{sample_path}/**/{pattern}", recursive=True)
            filtered = [f for f in found if sub_string in f]
            if filtered:
                return os.path.realpath(filtered[0]) 
            return None

        # Gather files from the Outdir
        c_file = find_file("*.counts.tsv.gz", "call-CountsMetrics")
        m_file = find_file("*.manta.std.vcf.gz", "execution")
        w_file = find_file("*.wham.std.vcf.gz", "execution")
        s_file = find_file("*.scramble.vcf.gz", "execution")

        if c_file and m_file and w_file and s_file:
            samples.append(sample_id)
            counts.append(c_file)
            manta.append(m_file)
            wham.append(w_file)
            scramble.append(s_file)
        else:
            print(f"SKIPPING: {sample_id} - One or more files missing in {sample_path}")

    data["EvidenceQC.samples"] = samples
    data["EvidenceQC.counts"] = counts
    data["EvidenceQC.manta_vcfs"] = manta
    data["EvidenceQC.wham_vcfs"] = wham
    data["EvidenceQC.scramble_vcfs"] = scramble

    with open('evidence_qc_inputs.json', 'w') as f:
        json.dump(data, f, indent=4)
    EOF

    # Execute Cromwell
    java -Xmx48G -Dconfig.file=${params.cromwell_conf} -jar ${params.cromwell_jar} \
        run ${params.evidqc_wdl} \
        -i evidence_qc_inputs.json \
        -p ${params.deps_zip}

    mkdir -p results
    mv cromwell-executions/EvidenceQC/*/call-* results/
    """
}