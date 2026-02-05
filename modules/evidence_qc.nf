#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process EVIDENCE_QC {

    input:
    val sample_list
    path "all_samples/*" 

    output:
    path "evidence_qc_results/*", emit: evidence_qc_results

    script:
    def template_path = file(params.evidqc_template ?: "${params.gatk_sv_dir}/data/EvidenceQC_inputs.json").toAbsolutePath()
    
    """
    python3 - <<EOF
    import json
    import os
    import glob
    import sys

    with open("${template_path}", 'r') as f:
        data = json.load(f)

    samples, counts, manta, wham, scramble = [], [], [], [], []

    search_base = "all_samples"

    for sample_id in ${sample_list.inspect()}:
        sample_path = os.path.join(search_base, sample_id)
        
        if not os.path.exists(sample_path):
             continue

        def find_file(pattern, sub_string):
            found = glob.glob(f"{sample_path}/**/{pattern}", recursive=True)
            filtered = [f for f in found if sub_string in f]
            return os.path.realpath(filtered[0]) if filtered else None

        c_file = find_file("*.counts.tsv.gz", "call-CollectCounts")
        m_file = find_file("*.manta.std.vcf.gz", "execution")
        w_file = find_file("*.wham.std.vcf.gz", "execution")
        s_file = find_file("*.scramble.vcf.gz", "execution")

        if c_file and m_file and w_file and s_file:
            samples.append(sample_id)
            counts.append(c_file)
            manta.append(m_file)
            wham.append(w_file)
            scramble.append(s_file)

    if not samples:
        print("ERROR: Zero samples found in staged directory!")
        sys.exit(1)

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

    mkdir -p evidence_qc_results
    cp evidence_qc_inputs.json evidence_qc_results/
    find cromwell-executions/EvidenceQC/ -name "call-*" -type d -exec cp -r {} evidence_qc_results/ \\;
    """
}