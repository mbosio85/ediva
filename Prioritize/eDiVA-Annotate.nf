#!/usr/bin/env nextflow

params.NAME = 'test'
params.VCF = 'combined.variants.supplement.vcf'
params.OUTF = '/users/so/mbosio/scratch/ediva_nextflow/'
params.CPU = '2'

VCF_f = file(params.VCF)

/*
 * Annnotate VCF
 */
process annotate {
    publishDir params.OUTF
    
    input:
    file VCF_f 

    output:
    file '*sorted.annotated.csv' into annotated_out
    
    shell:
    '''
     
    # annotation
    $PYTHON $EDIVA/Annotate/annotate.py --input !{VCF_f} --sampleGenotypeMode complete -f

    '''
}

/*
 * Rank
 */
process rank {
    publishDir params.OUTF, mode: 'move'
    
    input:
    file annotated_in from annotated_out

    output:
    file '*.ranked.csv' into ranked_out
    
    shell:
    '''
    OUTFILE=$(echo !{annotated_in}  | sed -e 's/sorted.annotated.csv/ranked.csv/g')
    # rank the variants given
    Rscript $EDIVA/Prioritize/wrapper_call.R $EDIVA/Prioritize/ediva_score.rds !{annotated_in}  $OUTFILE
 
    '''
}


