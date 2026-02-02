process EVIDENCE_QC {
    input:
    val sample_list
    // This tells Nextflow to collect all folders into a local directory called 'all_samples'
    path "all_samples/*" 

    output:
    path "qc_results/*", emit: qc_results

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

    # CHANGE 2: Point the search_base to the local staged folder
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

    mkdir -p qc_results
    # Using find to avoid 'argument list too long' errors
    find cromwell-executions/EvidenceQC/ -name "call-*" -type d -maxdepth 2 -exec mv -t qc_results/ {} +
    """
}