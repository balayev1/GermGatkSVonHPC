//
// SAMPLE_PROCESSING
//

params.options = [:]
include { PICARD_COLLECTINSERTSIZEMETRICS } from '../../modules/nf-core/picard/collectinsertsizemetrics/main.nf' addParams( options: params.options )

workflow SAMPLE_PROCESSING {
    take:
    
}