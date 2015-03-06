#!/usr/bin/env python
import pipeline_element
import subprocess
import pickle


pipe = list()
#Build a pipeline
text = ("""
### Align reads with bwa
$BWA mem -M -t 4  -R "@RG\\tID:$NAME\\tSM:$NAME" $REF $READ1 $READ2 | time $SAMTOOLS view -h -b -S -F 0x900 -  > $TMPDIR/$NAME.noChimeric.bam
### check for Quality encoding and transform to 33 if 64 encoding is encountered
OFFSET=$($SAMTOOLS view $TMPDIR/$NAME.noChimeric.bam | python $EDIVA/Predict/whichQuality_bam.py)
if [[ $OFFSET == 64 ]];
then
    echo "fixing 64 quality encoding"
    $SAMTOOLS view -h $TMPDIR/$NAME.noChimeric.bam | python $EDIVA/Predict/bam_rescale_quals.py - | $SAMTOOLS view -bS - > $TMPDIR/$NAME.transformed.bam
    rm $TMPDIR/$NAME.noChimeric.bam
    mv $TMPDIR/$NAME.transformed.bam $TMPDIR/$NAME.noChimeric.bam
fi
""")
p_element = pipeline_element.pipeline_element(text,"BWA alignment")
p_element.set_error("Error in changing encoding. Please refer to SGE job error file")
p_element.status = -1
pipe.append(p_element)


text =("""### Sort BAM file
if [ -s $TMPDIR/$NAME.noChimeric.bam ];
then
    echo Sort BAM
    $NOVOSORT --threads 4 --tmpdir $TMPDIR --forcesort --output $TMPDIR/$NAME.sort.bam -i -m 20G $TMPDIR/$NAME.noChimeric.bam
    cp $TMPDIR/$NAME.sort.bam* $OUTF/
    # clean up
    rm $TMPDIR/$NAME.noChimeric.bam
else
    echo $TMPDIR/$NAME.noChimeric.bam not found
    exit
fi""")
p_element = pipeline_element.pipeline_element(text,"BWA sorting")
p_element.set_error("$TMPDIR/$NAME.noChimeric.bam not found Please refer to SGE job error file")
pipe.append(p_element)

logfile = "log.log"
#Save it
pipeline_element.save_pipeline(pipe,logfile+'.pipe')


#Run it
with open("log.log",'a') as logfile:
    for p_el in pipe:
        logfile.write("Run " + p_el.alias +'\n')
        if p_el.status >=0:
            p_el.run()
            if p_el.status <0:
                logfile.write("Success for " + p_el.alias +'\n')
            else:
                logfile.write(p_el.err_msg +'\n')
                break
        else:
            logfile.write("Skip this step "+ p_el.alias +" because already executed successfully \n")
            print("Skip this step "+ p_el.alias +" because already executed successfully \n")

spipeline_element.save_pipeline(pipe,logfile+'.pipe')
        
        
     



#Done
