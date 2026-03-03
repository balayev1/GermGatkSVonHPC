//
// This file holds several functions specific to the main.nf workflow in the mskilab-org/nf-jabba pipeline
//

import nextflow.Nextflow

class WorkflowMain {
    //
    // Get attribute from genome config file
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }

    //
    // Populate params.* keys from the selected genome block when unset
    //
    public static void setGenomeDefaults(params, List<String> attributes) {
        attributes.each { attribute ->
            if (params[attribute] == null) {
                def value = getGenomeAttribute(params, attribute)
                if (value != null) {
                    params[attribute] = value
                }
            }
        }
    }
}
