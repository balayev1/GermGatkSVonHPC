#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRAIN_GCNV {

    input:
    val sample_list
    path "all_samples/*"

    output:
    path "train_gcnv_results", emit: train_gcnv_results

    script:
    def template_path = file(params.traingcnv_template ?: "${params.gatk_sv_dir}/data/TrainGCNV_inputs.json").toAbsolutePath()

    """
    python3 - <<EOF
    import json
    import os
    import glob
    import sys

    with open("${template_path}", 'r') as f:
        data = json.load(f)

    samples, counts = [], []

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
        
        if c_file:
            samples.append(sample_id)
            counts.append(c_file)
        
    if not samples:
        print("ERROR: Zero samples found in staged directory!")
        sys.exit(1)

    data["TrainGCNV.samples"] = samples
    data["TrainGCNV.count_files"] = counts

    with open('train_gcnv_inputs.json', 'w') as f:
        json.dump(data, f, indent=4)
    EOF

    # Execute Cromwell
    java -Xmx48G -Dconfig.file=${params.cromwell_conf} -jar ${params.cromwell_jar} \
        run ${params.traingcnv_wdl} \
        -i train_gcnv_inputs.json \
        -p ${params.deps_zip}
    
    mkdir -p train_gcnv_results
    cp train_gcnv_inputs.json train_gcnv_results/
    find cromwell-executions/TrainGCNV/ -name "*" -type d -exec cp -r {} train_gcnv_results/ \\;
    """
}
