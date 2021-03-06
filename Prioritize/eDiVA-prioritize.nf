#!/usr/bin/env nextflow

params.ANNOTATED = 'combined.variants.ranked.csv'
params.OUTF = '/users/so/mbosio/scratch/ediva_nextflow/'
params.MODE = 'standard'
params.PEDIGREE='pedigree.tree'
params.FAMILY_TYPE = 'trio'
params.HPO_list='/nfs/users2/GD/tools/ediva/Resource/empty_white_list.txt'
params.EXCLUSIONLIST='/nfs/users2/GD/tools/ediva/Resource/gene_exclusion_list.txt'
params.INHERITANCE='none'
ANNOTATED_f = file(params.ANNOTATED)
PED_f = file(params.PEDIGREE)
EXCLUDED = file(params.EXCLUSIONLIST)
HPO = file(params.HPO_list)
XLS = file(params.OUTF+"/" +'/variant_prioritization_report.xlsx' )
/*
 * Prioritize ranked
 */
process prioritize {
    publishDir params.OUTF+"/" + params.INHERITANCE, mode :'copy'
    
    input:
    file ANNOTATED_f
    file PED_f
    file EXCLUDED
    file HPO
     
    file XLS
    
    val INHERITANCE from params.INHERITANCE
    val type from params.FAMILY_TYPE
    val mode from params.MODE
    
    output:
    file "*${INHERITANCE}*" into prioritzed
    file '../variant_prioritization_report.xlsx' into xls2
    
    
    shell:
    if (['dominant_inherited','recessive','dominant_denovo','Xlinked','compound'].contains(INHERITANCE)){
        if (['trio','family'].contains(type)){    

                
                if (['standard','strict'].contains(mode)){
                    if (mode =='standard'){
                        '''
                        # run inheritance inheritance: dominant_inherited        
                        OUT1=$(echo !{ANNOTATED_f} | sed -e "s/ranked.csv/!{INHERITANCE}.csv/g")
                        OUT2=$(echo !{ANNOTATED_f} | sed -e "s/ranked.csv/filtered.!{INHERITANCE}.csv/g")
                        $PYTHON $EDIVA/Prioritize/familySNP.py --infile !{ANNOTATED_f} \
                                --outfile $OUT1 \
                                --filteredoutfile $OUT2  \
                                --family !{PED_f} \
                                --inheritance !{INHERITANCE}\
                                --familytype trio \
                                --geneexclusion !{EXCLUDED} \
                                --white_list !{HPO}
                        cp variant_prioritization_report.xlsx ../variant_prioritization_report.xlsx
                        '''
                    }
                    else{
                                                '''
                        # run inheritance inheritance: dominant_inherited        
                        OUT1=$(echo !{ANNOTATED_f} | sed -e "s/ranked.csv/!{INHERITANCE}.csv/g")
                        OUT2=$(echo !{ANNOTATED_f} | sed -e "s/ranked.csv/filtered.!{INHERITANCE}.csv/g")
                        $PYTHON $EDIVA/Prioritize/familySNP_2.0.py --infile !{ANNOTATED_f} \
                                --outfile $OUT1 \
                                --filteredoutfile $OUT2  \
                                --family !{PED_f} \
                                --inheritance !{INHERITANCE}\
                                --familytype trio \
                                --geneexclusion !{EXCLUDED} \
                                --white_list !{HPO}
                      cp variant_prioritization_report.xlsx ../variant_prioritization_report.xlsx 
                        '''
                    }
                }
                
                else{       
                error "Mode ${mode} is not available in eDiVA-prioritize, please choose one from:standard,strict"
                }
        }
        else{       
        error "FAMILY_TYPE ${type} is not available in eDiVA-prioritize, please choose one from:trio,family"
        }
    }
    else{
        error "Inheritance ${INHERITANCE} is not available in eDiVA-prioritize, please choose one from :dominant_inherited,recessive,dominant_denovo,Xlinked,compound"    
    }
}
         


