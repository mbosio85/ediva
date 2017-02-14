#!/usr/bin/env nextflow

params.NAME = 'test'
params.FAMILY_INFO = 'family_config.txt'
params.OUTF = '/users/so/mbosio/scratch/ediva_nextflow/'
params.CPU = '2'

FAM_INFO = file(params.FAMILY_INFO)

/*
 * Merge VCF
 */
process mergeVCF {
    publishDir params.OUTF, mode: 'move'
    
    input:
    file FAM_INFO

    output:
    file 'combined.variants.supplement.vcf' into merge_vcf
    file 'pedigree.tree' into ped_file
    
    shell:
    '''
    
    grep -v -P 'ID\tstatus' !{FAM_INFO}  | awk -F '\t' '{ print $4}'> bam.list
    echo "sample affected" > pedigree.tree
    grep -v -P 'ID\tstatus' !{FAM_INFO}  | awk -F '\t' '{ print $1,$2}' >> pedigree.tree
    sed -i 's/ /\t/g'  pedigree.tree
    
    # merge vcf files
    VARLINE=$(echo $(grep -v -P 'ID\tstatus' !{FAM_INFO}  |  awk -F '\t' '{print "--variant",$3}'))
    java -jar $GATK -T CombineVariants -R $REF $VARLINE -o combined.variants.vcf --unsafe LENIENT_VCF_PROCESSING
            
     # do Genotyping in all family members
     java -jar $GATK -T UnifiedGenotyper -R $REF -I bam.list --dbsnp $DBSNP -o combined.variants.supplement.temp.vcf -alleles combined.variants.vcf --output_mode EMIT_ALL_SITES --genotyping_mode GENOTYPE_GIVEN_ALLELES -glm BOTH
        
      $PYTHON $EDIVA/Prioritize/vcf_filter.py --infile combined.variants.supplement.temp.vcf --outfile combined.variants.supplement.vcf
       
      rm combined.variants.supplement.temp.vcf
    
    '''

    
    
}