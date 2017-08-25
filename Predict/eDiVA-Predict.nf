#!/usr/bin/env nextflow

params.NAME = 'test'
params.READ1 = '/users/so/nrostan/PHD_exome_course_2015/TheExomeCourse/Data/Case_1_FHHt/CD2224/CD2224.read1.fastq.gz'
params.READ2 = '/users/so/nrostan/PHD_exome_course_2015/TheExomeCourse/Data/Case_1_FHHt/CD2224/CD2224.read2.fastq.gz'
params.OUTF = '/users/so/mbosio/scratch/ediva_nextflow/'
params.CPU = '2'
params.AFFECTED = '0'

READ1 = file(params.READ1)
READ2 = file(params.READ2)
File AA = new File(params.OUTF)
String OUTF_FULL_PATH = AA.absolutePath
/*
 * Align fastq files
 */
process alignReads {
    //publishDir params.OUTF+"/Intermediate_files"
	cpus params.CPU
    input:
    file READ1 
    file READ2
    val NAME from params.NAME
    //val CPU from params.CPU
    

    output:
    file "${NAME}.noChimeric.bam" into chimeric_ch
    
    shell:
    '''
    $BWA mem -M -t !{task.cpus} -R \"@RG\\tID:!{NAME}\\tSM:!{NAME}\" $REF !{READ1} !{READ2} |  $SAMTOOLS view -h -b -S -F 0x900 -  > !{NAME}.noChimeric.bam 
    ### check for Quality encoding and transform to 33 if 64 encoding is encountered
    OFFSET=$($SAMTOOLS view !{NAME}.noChimeric.bam  | $PYTHON $EDIVA/Predict/whichQuality_bam.py)
    if [[ $OFFSET == 64 ]];
    then
        echo 'fixing 64 quality encoding'
        $SAMTOOLS view -h !{NAME}.noChimeric.bam  | $PYTHON $EDIVA/Predict/bam_rescale_quals.py - | $SAMTOOLS view -bS - > !{NAME}.transformed.bam
        rm !{NAME}.noChimeric.bam
        mv !{NAME}.transformed.bam !{NAME}.noChimeric.bam 
    fi

    '''
}

/*
 * Sort Bam file
 */
process sortBam {
    //publishDir params.OUTF+"/Intermediate_files"
	cpus params.CPU
    input:
    val NAME from params.NAME
    val OUTF from params.OUTF
    
    file bam from chimeric_ch

    output:
    file "${NAME}.sort.bam" into sorted_bam
    file "${NAME}.sort.bam.bai" into sorted_idx

    shell:
    '''
    echo Sort BAM
    # NOVOSORT --threads !{task.cpus} --tmpdir !{OUTF} --forcesort --output !{NAME}.sort.bam -i -m 40G !{NAME}.noChimeric.bam
    $SAMTOOLS sort -@  !{task.cpus} -T !{NAME}.tmp_sorting   --reference  $REF  -o !{NAME}.sort.bam  !{NAME}.noChimeric.bam
    $SAMTOOLS index !{NAME}.sort.bam 
    '''
}

/*
 *  Local Re-alignment
 */
process local_realignment {
    //publishDir params.OUTF+"/Intermediate_files"
        cpus params.CPU
        input:
        val NAME from params.NAME
        val OUTF from params.OUTF
        //val CPU from params.CPU
        file sorted_bam_in from sorted_bam
        file sorted_idx_in from sorted_idx
        
        output:
        file "${NAME}.realigned.bam" into realigned_bam
        file "${NAME}.realigned.bai" into realigned_idx
    
        shell:
        '''
               echo Local Re-alignment
               java -jar $GATK -nt  !{task.cpus} -T RealignerTargetCreator -R $REF -I !{NAME}.sort.bam -o !{NAME}.intervals -known $DBINDEL --minReadsAtLocus 6 --maxIntervalSize 200 --downsampling_type NONE # -L $EXOME
               java -jar $GATK  -T IndelRealigner -R $REF -I !{NAME}.sort.bam -targetIntervals !{NAME}.intervals -o !{NAME}.realigned.bam -known $DBINDEL --maxReadsForRealignment 10000 --consensusDeterminationModel USE_SW --downsampling_type NONE # -L $EXOME
    '''
    }
    
/*
 *  Duplicate marking
 */
process duplicate_marking {
    //publishDir params.OUTF+"/Intermediate_files"
        cpus params.CPU
        input:
        val NAME from params.NAME
        val OUTF from params.OUTF
        file realign_bam from realigned_bam
        file realign_idx from realigned_idx
        
        output:
        file "${NAME}.realigned.dm.bam" into dm_bam_out
        file "${NAME}.realigned.dm.bai" into dm_idx_out
    
        shell:
        '''
               echo Duplicate marking
               java -jar $PICARD/MarkDuplicates.jar INPUT=!{NAME}.realigned.bam OUTPUT=!{NAME}.realigned.dm.bam METRICS_FILE=duplication_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
        '''
}          
            
/*
 *  Base quality recalibration
 */ 
process bam_recalibrate {
    publishDir params.OUTF, mode:'copy'
        cpus params.CPU
        input:
        val NAME from params.NAME
        val OUTF from params.OUTF
        file dm_bam_in from dm_bam_out
        file dm_idx_in from dm_idx_out
        
        output:
        file "${NAME}.realigned.dm.recalibrated.bam" into realigned_bam_out,realigned_bam_out2,realigned_bam_out3
        file "${NAME}.realigned.dm.recalibrated.bai" into realigned_idx_out,realigned_idx_out2,realigned_idx_out3
    
        shell:
        '''
               echo -e " \n #### Base quality recalibration \n "
               java -jar $GATK -T BaseRecalibrator -nct  !{task.cpus} --default_platform illumina -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R $REF -I !{NAME}.realigned.dm.bam -knownSites $DBSNP --downsampling_type NONE -o !{NAME}.recal_data.grp # -L $EXOME
               java -jar $GATK -T PrintReads -R $REF -I !{NAME}.realigned.dm.bam -BQSR !{NAME}.recal_data.grp -o !{NAME}.realigned.dm.recalibrated.bam # -L $EXOME
        '''

}             

            
/*
 *  GATK: Call SNPs and Indels with the GATK Haplotype Caller
 */

process GATK_call {
    //publishDir params.OUTF
    
        cpus params.CPU
        input:
        val NAME from params.NAME
        val OUTF from params.OUTF
        file realigned_bam_in from realigned_bam_out
        file realigned_idx_in from realigned_idx_out
        
        output:
        file "GATK.snps.raw.vcf" into gatk_raw_snps
        file "GATK.indel.raw.vcf" into gatk_raw_indel,gatk_raw_indel2
    
        shell:
        '''
               echo -e "\n #### GATK: Call SNPs and Indels with the GATK Unified Genotyper \n"
               java -jar $GATK -T HaplotypeCaller -nct 1 -R $REF --dbsnp $DBSNP -I !{NAME}.realigned.dm.recalibrated.bam -o GATK.both.raw.vcf -L $EXOME
            
               echo -e "\n #### GATK: Split SNPs and Indels \n"
               java  -jar $GATK -T SelectVariants -R $REF --variant GATK.both.raw.vcf -o GATK.snps.raw.vcf  -selectType SNP
               java  -jar $GATK -T SelectVariants -R $REF --variant GATK.both.raw.vcf -o GATK.indel.raw.vcf -selectType INDEL
        '''
}


/*
 * Filter and compare SNP calls from 2 different pipelines
 */
process GATK_filtering {
    //publishDir params.OUTF+"/SNP_Intersection"
        cpus params.CPU
        input:
        val NAME from params.NAME
        val OUTF from params.OUTF
        file raw_snps from gatk_raw_snps
        file raw_imdel from gatk_raw_indel2
        
        output:
        file "GATK.snps.filtered.cleaned.vcf" into gatk_clean_snps
        file "report.snps.txt" into snps_report
        file "snps.enriched.vcf" into snps_enriched
        file "report.snps.enriched.txt" into report_enriched_snps
    
        shell:
        '''
            # Filtering           
            java -jar $GATK -T VariantFiltration -R $REF -o GATK.snps.filtered.vcf --variant GATK.snps.raw.vcf --mask GATK.indel.raw.vcf --clusterWindowSize 10 --filterExpression \"MQ < 30.0 || QUAL < 25.0 \" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > 400 || GQ < 15\" --genotypeFilterName CRGg
            STATUS=\"${?}\"
            
            if [ \"$STATUS\" -gt 0 ];
            then
                echo Error in GATK Filtering VariantFiltration >&2
                exit 1
            fi
            
            ## isolate PASSed variants
             grep -E '^#|PASS' GATK.snps.filtered.vcf | grep -v CRGg > GATK.snps.filtered.cleaned.vcf

            # Evaluation
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBSNP  -o report.snps.txt --eval GATK.snps.filtered.cleaned.vcf -l INFO # -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"GATK\"' -selectName GATK
            STATUS=\"${?}\"
            
            if [ \"$STATUS\" -gt 0 ];
            then
                echo Error in GATK VariantEval >&2
                exit 1
            fi

            ## Annotate Enrichment
            $BEDTOOLS/intersectBed -a GATK.snps.filtered.cleaned.vcf -b $EXOME > merged.all.vcf
            
            # borrow header from GATK vcf file
            grep '^#' GATK.snps.filtered.cleaned.vcf > snps.enriched.vcf
            cat merged.all.vcf >> snps.enriched.vcf
            
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBSNP -o report.snps.enriched.txt --eval snps.enriched.vcf -l INFO # -select 'set==\"Intersection\"' -selectName Intersection -select 'set==\"GATK\"' -selectName GATK
            
            STATUS=\"${?}\"
            echo $STATUS
            if [ \"$STATUS\" -gt 0 ];
            then
                echo Error in Annotate VariantEval GATK  >&2
                exit 1
            fi
            '''
            
}



/*
 * GATK Indel preparation
 */
process GATK_filtering_indel {
//    publishDir params.OUTF+"/Indel_Intersection"
    cpus params.CPU
    input:
    val NAME from params.NAME
    val OUTF from params.OUTF
    file raw_indels from gatk_raw_indel
        
        output:
        file "GATK.indel.filtered.cleaned.vcf" into gatk_clean_indel
        file "indel.enriched.vcf" into indel_enriched
        file "report.indel.enriched.txt" into report_enriched_indel
    
        shell:
        '''
            
            java -jar -Xmx4g $GATK -T VariantFiltration -R $REF -o GATK.indel.filtered.vcf --variant GATK.indel.raw.vcf --filterExpression \"QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20\" --filterName CRG --genotypeFilterExpression \"DP < 5 || DP > 400 || GQ < 15\" --genotypeFilterName LOWQ
            STATUS=\"${?}\"
            if [ \"$STATUS\" -gt 0 ];
            then
                echo Error in GATK VariantFiltration  >&2
                exit 1
            fi
            # select PASSed variants
            grep -v \"CRG\" GATK.indel.filtered.vcf > GATK.indel.filtered.cleaned.vcf
            
            
            # filter Indels for enriched regions
            $BEDTOOLS/intersectBed -a GATK.indel.filtered.cleaned.vcf -b $EXOME > merged.all.vcf
            
            grep '^#' GATK.indel.filtered.cleaned.vcf > indel.enriched.vcf
            cat merged.all.vcf >> indel.enriched.vcf
            
            java -Xmx5g -jar $GATK -T VariantEval -R $REF --dbsnp $DBINDEL -select 'set==\"GATK\"' -selectName GATK -o report.indel.enriched.txt --eval indel.enriched.vcf -l INFO
            STATUS=\"${?}\"
            if [ \"$STATUS\" -gt 0 ];
            then
                echo Error in filter Indels for enriched regions >&2
                exit 1
            fi
        '''
}


/*
 * Fuse variants
 */
process Fuse_variants {
    publishDir params.OUTF,mode :'copy'
        cpus params.CPU
        input:
        val NAME from params.NAME
        val OUTF from params.OUTF
        file merge_indel from indel_enriched
        file merge_snps from snps_enriched
        
        output:
        file "all_variants.vcf" into final_vcf
    
        shell:
        '''
            #FUSEVARIANTS    
            # fuse indel and snp calls into one file
            java -Xmx5g -jar $GATK -T CombineVariants -R $REF --variant snps.enriched.vcf --variant indel.enriched.vcf -o all_variants.vcf -U LENIENT_VCF_PROCESSING --genotypemergeoption UNIQUIFY
            STATUS=\"${?}\"
            if [ \"$STATUS\" -gt 0 ];
            then
                echo Error in fusevariants >&2
                exit 1
            fi
        '''  
}

/*
 * Quality control
 */

 process Quality_control{
    publishDir params.OUTF+"/QualityControl"
        cpus params.CPU
        input:
	file READ1
	file READ2
        val NAME from params.NAME
        val OUTF from params.OUTF
        file real_bam_qc from realigned_bam_out2
        file real_bai_qc from realigned_idx_out2
        file snp_report from report_enriched_snps

        output:
        file "*fastqc.zip"  into reads_fastqc
        file 'coverage.pdf'             into coverage_plot
        file 'flagstat_kit.txt'         into flagstat
        file 'flagstat_report.txt'      into flagstat_rep
        file 'gatk_report.txt'          into gatk_rep
    
        shell:
        '''
            ## QUALITY CONTROL
            $PYTHON   $EDIVA/Predict/quality_control.py -bamfile !{NAME}.realigned.dm.recalibrated.bam   -bedfile $EXOME   --output ./   --samtools $SAMTOOLS   --bedtools $BEDTOOLS   --fastqc $FASTQC   --read1 !{READ1}   --read2 !{READ2}   --gatk  report.snps.enriched.txt
        '''
 }

 
 
 
 /*
 * Cleanup
 */

process Cleanup{
    publishDir params.OUTF, mode: 'copy'
    input:
        val NAME from params.NAME
        val OUTF from params.OUTF
        val CPU from params.CPU
        val AFFECTED from params.AFFECTED
        file qc_rep from gatk_rep
        val OUTF_FULL_PATH from OUTF_FULL_PATH
        //file real_bam_qc from realigned_bam_out3
       // file real_bai_qc from realigned_idx_out3
        //file allvar from final_vcf
   output:
        //file "all_variants.vcf" into final_vcf2
        //file "${NAME}.realigned.dm.recalibrated.bam" into realigned_bam_final
        //file "${NAME}.realigned.dm.recalibrated.bai" into realigned_idx_final
        file "sample_info.txt" into sample_info
	file "sample_info_docker.txt" into sample_info_docker
   shell :
   '''
   
   echo \"!{NAME} !{AFFECTED} !{OUTF_FULL_PATH}/all_variants.vcf !{OUTF_FULL_PATH}/!{NAME}.realigned.dm.recalibrated.bam "|sed -e 's/ /\t/g' > sample_info.txt 
   echo \"!{NAME} !{AFFECTED} /samples/!{NAME}/all_variants.vcf /samples/!{NAME}/!{NAME}.realigned.dm.recalibrated.bam "|sed -e 's/ /\t/g' > sample_info_docker.txt 
   
   '''
 
 
 }

