########################### SCRIPT ##############################

#!/usr/bin/env python

infodist = { 
	'P479_101_ThruPlex_Streptococcus_Pneumoniae'           : [ 'index', 'XXXXXX', '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_101_ThruPlex_Streptococcus_Pneumoniae/P479_101_ThruPlex_Streptococcus_Pneumoniae_ATCACGAT-TCTTTCCC_L001_R1_001.cutadapt.ar.fastq','/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_101_ThruPlex_Streptococcus_Pneumoniae/P479_101_ThruPlex_Streptococcus_Pneumoniae_ATCACGAT-TCTTTCCC_L001_R2_001.cutadapt.ar.fastq' ],
    	'P479_101_XTrobot_Streptococcus_Pneumoniae'           : [ 'index', 'XXXXXX', '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_101_XTrobot_Streptococcus_Pneumoniae/P479_101_XTrobot_Streptococcus_Pneumoniae_GGACTCCT-AGAGTAGA_L001_R1_001.cutadapt.ar.fastq','/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_101_XTrobot_Streptococcus_Pneumoniae/P479_101_XTrobot_Streptococcus_Pneumoniae_GGACTCCT-AGAGTAGA_L001_R2_001.cutadapt.ar.fastq' ],
    	'P479_102_ThruPlex_Streptococcus_Pneumoniae'           : [ 'index', 'XXXXXX', '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_102_ThruPlex_Streptococcus_Pneumoniae/P479_102_ThruPlex_Streptococcus_Pneumoniae_CGATGTAT-TCTTTCCC_L001_R1_001.cutadapt.ar.fastq','/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_102_ThruPlex_Streptococcus_Pneumoniae/P479_102_ThruPlex_Streptococcus_Pneumoniae_CGATGTAT-TCTTTCCC_L001_R2_001.cutadapt.ar.fastq' ],
    	'P479_102_XTrobot_Streptococcus_Pneumoniae'           : [ 'index', 'XXXXXX', '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_102_XTrobot_Streptococcus_Pneumoniae/P479_102_XTrobot_Streptococcus_Pneumoniae_CGTACTAG-AGAGTAGA_L001_R1_001.cutadapt.ar.fastq','/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Sample_P479_102_XTrobot_Streptococcus_Pneumoniae/P479_102_XTrobot_Streptococcus_Pneumoniae_CGTACTAG-AGAGTAGA_L001_R2_001.cutadapt.ar.fastq' ],

#    'Thruplex_i9'           : [ 'index', 'XXXXXX', '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/index_9_read1_cutadapt.fastq','/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/index_9_read2_cutadapt.fastq' ],
# 
 #  'Klebsiella_16_ThruPlex'           : [ 'index', 'XXXXXX', '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/BWA_test/Sample_16_ThruPlex_Klebsiella_ACAGTGAT-TCTTTCCC_L001_R1_001.cutadapt.ar.fastq','/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/BWA_test/Sample_16_ThruPlex_Klebsiella_ACAGTGAT-TCTTTCCC_L001_R2_001.cutadapt.ar.fastq' ],
 #   'bm55.66'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.66.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.66.r2.fastq.gz' ],
 #   'bm55.68'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.68.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.68.r2.fastq.gz' ],
 #   'bm55.71'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.71.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.71.r2.fastq.gz' ],
 #   'bm55.neg'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.neg.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.neg.r2.fastq.gz' ],
 #   'bm55.11'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.11.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.11.r2.fastq.gz' ],
 #   'bm55.24'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.24.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.24.r2.fastq.gz' ],
 #   'bm55.32'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.32.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.32.r2.fastq.gz' ],
 #   'bm55.33'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.33.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.33.r2.fastq.gz' ],
 #   'bm55.36'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.36.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.36.r2.fastq.gz' ],
 #   'bm55.39'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.39.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.39.r2.fastq.gz' ],
 #   'bm55.45'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.45.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.45.r2.fastq.gz' ],
 #   'bm55.49'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.49.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.49.r2.fastq.gz' ],
 #   'bm55.55'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.55.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.55.r2.fastq.gz' ],
 #   'bm55.62'           : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.62.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/bm55.62.r2.fastq.gz' ],
 #   'bm55.Donor'        : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/DON.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/DON.r2.fastq.gz' ],
 #   'bm55.Patient'      : [ 'index', 'XXXXXX', '/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/PAT.r1.fastq.gz','/proj/b2011011/fatCellTnt-EB/data/bm55.Concat/PAT.r2.fastq.gz' ]
    }

def main():
    import sys
    samples = []
    for sampleid in infodist.keys():
	s = Sample(sampleid,sys.argv[1])
	s.createDirs()
	s.makeScriptMapping()
	#s.makeScriptReAlNreCal()
	#s.makeScriptQC()
	samples.append(s)

    c = collection(samples,sys.argv[2])
    c.createDirs()
    c.makeScript_haplotypeCaller()
    c.makeScript_vqsr_indels()
    c.makeScript_vqsr_snps()
    
class Sample(object):

	def __init__(self,sampleid,basepath):
	    import os
	    self.sampleid = sampleid
	    self.index = [infodist[self.sampleid][0]]
	    self.indexseq = [infodist[self.sampleid][1]]
	    self.r1files = [infodist[self.sampleid][2]]
	    self.r2files = [infodist[self.sampleid][3]]
	    self.basepath = os.path.abspath(basepath)
	    self.path = os.path.abspath(basepath)+'/'+sampleid
	    ######################################### HaR e dina grejjer
            self.project = 'b2013064'
	    self.reference = '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Streptococcus_samples/Reference/Streptococcus_ref.fasta'
	    self.gatkreference = '/proj/b2013064/britta/NXT_deep/140305_M01320_0084_000000000-A8FKA/Data/Intensities/Unaligned/Project_NXT/Streptococcus_samples/Reference/Streptococcus_ref.fasta'
            ############################################################
	    self.scriptpath =	self.path+'/script'
	    self.sam =		self.path+'/'+self.sampleid+'.bwa.sam'
	    self.bam =		self.path+'/'+self.sampleid+'.bwa.bam'
	    self.sortedbam =	self.path+'/'+self.sampleid+'.bwa.sorted.bam'
	    self.mdmetrix =	self.path+'/'+self.sampleid+'.MarkDupsMetrix'
	    self.markedbam =	self.path+'/'+self.sampleid+'.bwa.marked.bam'
	    self.rginfobam =	self.path+'/'+self.sampleid+'.bwa.rgInfoFixed.bam'
	    self.realtargets =	self.path+'/'+self.sampleid+'.reAlignemntTargetIntervals.bed'	    
	    self.realignedbam = self.path+'/'+self.sampleid+'.realigned.bam'
	    self.bqsr =		self.path+'/'+self.sampleid+'.BQSR.grp'
	    self.recalbam =	self.path+'/'+self.sampleid+'.recal.final.bam'
	    self.realncalscript= self.scriptpath +'/'+self.sampleid+'.realNrecal.sh'
	    self.mappingscript = self.scriptpath +'/'+self.sampleid+'.mapping.sh'
	    self.qcscript = self.scriptpath +'/'+self.sampleid+'.qc.sh'
	    self.extraqcscript = self.scriptpath +'/'+self.sampleid+'.extra_qc.sh'

	def createDirs(self):
	    import os
	    try: os.makedirs(self.path)
	    except OSError:pass
	    try: os.makedirs(self.scriptpath)
	    except OSError:pass
	
	def makeScriptQC(self):
	    output = ''
	    output += '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+self.project+'\n'
	    output += '#SBATCH -n 16 -p node'+'\n'
	    output += '#SBATCH -t 72:00:00'+'\n'
	    output += '#SBATCH -J '+self.sampleid+'.QC'+'\n'
	    output += '#SBATCH -e '+self.path+'/stderr.QC.txt'+'\n'
	    output += '#SBATCH -o '+self.path+'/stdout.QC.txt'+'\n'
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'cd '+self.path+'\n'
	    output += 'picard=/proj/b2011168/private/bin/picard-tools-1.63'+'\n'
	    output += 'GATK=/proj/b2011011/fatCellTnt-EB/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar'+'\n'
	    output += 'GATKbundle=/proj/b2011011/fatCellTnt-EB/references/GATKbundle/'+'\n'
	    output += 'file='+self.recalbam+'\n'
	    output += 'echo "-----"'+'\n'

	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'echo -e "-> CallableLoci <-"'+'\n'
	    output += 'java -Xmx60g -jar $GATK -T CallableLoci -I $file -summary $file.summary.txt -o $file.callableLoci.bed -R '+self.gatkreference+' &'+'\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
	    output += 'echo'+'\n'
	    output += 'echo "-----"'+'\n'

	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'echo -e "-> flagstat <-"'+'\n'
	    output += 'samtools flagstat $file > $file.flagstat.txt &'+'\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
	    output += 'echo'+'\n'
	    output += 'echo "-----"'+'\n'

	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'echo -e "-> Pauls qacompute <-"'+'\n'
	    output += '/proj/b2010052/scripts/qaCompute -d -q 10 -m $file $file.qacompute.out > $file.qacompute.stdout.txt 2> $file.qacompute.stderr.txt &'+'\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
	    
	    output += 'java -Xmx60g -jar $picard/CalculateHsMetrics.jar BAIT_INTERVALS=/proj/b2011011/fatCellTnt-EB/references/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem TARGET_INTERVALS=/proj/b2011011/fatCellTnt-EB/references/truseq_exome_targeted_regions.hg19.bed.chr.columnReOrdered.withHeader.chrRem INPUT=$file OUTPUT=$file.hs_metrics.summary.txt PER_TARGET_COVERAGE=$file.hs_metrics.perTargetCoverage.txt REFERENCE_SEQUENCE='+self.reference+'\n'

	    output += 'echo'+'\n'
	    output += 'wait'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'
	    with open(self.qcscript,'w') as outfile: outfile.write(output)
	
	def makeScriptMapping(self):
	    output = '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+self.project+'\n'
	    output += '#SBATCH -n 16 -p node'+'\n'
	    output += '#SBATCH -t 5:00:00'+'\n'
	    output += '#SBATCH -J '+self.sampleid+'.map.mark'+'\n'
	    output += '#SBATCH -e '+self.path+'/stderr.mapNmark.txt'+'\n'
	    output += '#SBATCH -o '+self.path+'/stdout.mapNmark.txt'+'\n'
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'cd '+self.path+'\n'
	    output += 'echo'+'\n'
            #output += 'module load bioinfo-tools samtools picard/1.92 bowtie2/2.1.0'+'\n'
	    output += 'module load bioinfo-tools samtools picard/1.92 bwa'+'\n'
	    output += 'picard=/sw/apps/bioinfo/picard/1.92/milou'+'\n'
	    
            #output += 'bowtie2 -1 '+self.r1files[0]+' -2 '+self.r2files[0]+' -p16 -x '+self.reference+' > '+self.sam+'\n'
	    output += 'bwa mem '+self.reference+' '+self.r1files[0]+' '+self.r2files[0]+'  > '+self.sam+'\n'
	    output += 'echo -e "mapping Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    
	    output += 'java -Xmx60g -jar $picard/SamFormatConverter.jar \\'+'\n'
	    output += 'MAX_RECORDS_IN_RAM=2500000 INPUT='+self.sam+' OUTPUT='+self.bam+'\n'
	    output += 'echo -e "sam2bam Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+self.sam+''+'\n'
	    
	    #output += 'samtools view -bS -t '+self.reference+' '+self.sam+' > '+self.bam+'\n'
	    #output += 'echo -e "sam2bam Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    #output += 'rm -v '+self.sam+''+'\n'
	    
	    output += 'java -Xmx60g -jar $picard/SortSam.jar \\'+'\n'
	    output += 'MAX_RECORDS_IN_RAM=2500000 SORT_ORDER=coordinate INPUT='+self.bam+' '
	    output +='OUTPUT='+self.sortedbam+' CREATE_INDEX=true'+'\n'
	    output += 'echo -e "bam2sort Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+self.bam+''+'\n'
	    
	    #output += 'samtools sort '+self.bam+' '+self.sortedbam+'\n'
	    #output += 'echo -e "bam2sort Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    #output += 'rm -v '+self.bam+''+'\n'
	    
	    output += 'java -Xmx60g -jar $picard/MarkDuplicates.jar MAX_RECORDS_IN_RAM=2500000 VALIDATION_STRINGENCY=LENIENT INPUT='+self.sortedbam+' OUTPUT='+self.markedbam+' METRICS_FILE='+self.mdmetrix+'\n'
	    output += 'echo -e "mark Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+self.sortedbam+''+'\n'
	    
	    output += 'samtools flagstat '+self.markedbam+' \n'
	    output += 'echo "flagstat Done. $(date) Running on: $(hostname)" 1>&2'+'\n'

	    output += 'java -Xmx60g -jar $picard/CollectInsertSizeMetrics.jar MAX_RECORDS_IN_RAM=2500000 VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE='+self.markedbam+'.histogram.isize.pdf OUTPUT='+self.markedbam+'.isize.out INPUT='+self.markedbam+'\n'
	    output += 'echo -e "mark Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+self.sortedbam+''+'\n'
	    
	    output += 'java -Xmx60g -jar $picard/CollectGcBiasMetrics.jar MAX_RECORDS_IN_RAM=2500000 VALIDATION_STRINGENCY=LENIENT REFERENCE_SEQUENCE='+self.reference+' OUTPUT='+self.markedbam+'.GCbias.out INPUT='+self.markedbam+' CHART_OUTPUT='+self.markedbam+'.GCbias.pdf \n'
	    output += 'echo -e "GCbias Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+self.sortedbam+''+'\n'
	    
	    output += 'java -Xmx60g -jar $picard/AddOrReplaceReadGroups.jar \\'+'\n'
	    output += 'MAX_RECORDS_IN_RAM=2500000 \\'+'\n'
	    output += 'INPUT='+self.markedbam+' \\'+'\n'
	    output += 'OUTPUT='+self.rginfobam+' \\'+'\n'
	    output += 'CREATE_INDEX=true RGID='+self.sampleid+' RGLB='+self.sampleid+' RGPL=illumina RGSM='+self.sampleid+' RGCN="SciLifeLab" RGPU="barcode"'+'\n'
	    output += 'echo "addorreplace Done. $(date) Running on: $(hostname)" 1>&2'+'\n'
	    output += 'rm -v '+self.markedbam+''+'\n'
	    
	    output += 'wait'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'
	    with open(self.mappingscript,'w') as outfile: outfile.write(output)
	
	def makeScriptReAlNreCal(self):
	    output = '#! /bin/bash -l\n'
	    output += '#SBATCH -A '+self.project+'\n'
	    output += '#SBATCH -n 16 -p node\n'
	    output += '#SBATCH -t 72:00:00\n'
	    output += '#SBATCH -J '+self.sampleid+'.realNrecal\n'
	    output += '#SBATCH -e '+self.path+'/stderr.realNrecal.txt\n'
	    output += '#SBATCH -o '+self.path+'/stdout.realNrecal.txt\n'
	    output += '#SBATCH --mail-type=All\n'
	    output += '#SBATCH --mail-user=erik.borgstrom@scilifelab.se\n'
	    output += 'echo "$(date) Running on: $(hostname)"\n'

	    output += 'cd $workpath\n'
	    output += 'picard=/proj/b2011168/private/bin/picard-tools-1.63\n'
	    output += 'GATK=/proj/b2011011/fatCellTnt-EB/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar\n'
	    output += 'GATKbundle=/proj/b2011011/fatCellTnt-EB/references/GATKbundle/\n'
	    output += 'echo "$(date) Running on: $(hostname)"\n'
	    
	    output += 'echo -e "-> RealignerTargetCreator <-"\n'
	    output += 'java -Xmx72g -jar $GATK -T RealignerTargetCreator -nt 16 -I '+self.rginfobam+' -R '+self.gatkreference+' -o '+self.realtargets
	    output += ' -known $GATKbundle/Mills_and_1000G_gold_standard.indels.b37.vcf'
	    output += ' -known $GATKbundle/1000G_phase1.indels.b37.vcf;'
	    output += '\n'
	    
	    output += 'echo -e "-> IndelRealigner <-"\n'
	    output += 'java -Xmx72g -jar $GATK -T IndelRealigner -I '+self.rginfobam+' -R '+self.gatkreference+' -targetIntervals '+self.realtargets
	    output += ' -o '+self.realignedbam
	    output += ' -known $GATKbundle/Mills_and_1000G_gold_standard.indels.b37.vcf'
	    output += ' -known $GATKbundle/1000G_phase1.indels.b37.vcf;'
	    output += '\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"\n'
	    output += 'echo "$(date) Running on: $(hostname)"\n'
	    output += 'rm -v '+self.rginfobam+''+'\n'

	    output += 'echo -e "-> BaseRecalibrator <-"\n'
	    output += 'java -Xmx72g -jar $GATK -T BaseRecalibrator -I '+self.realignedbam+' -R '+self.gatkreference+' -o '+self.bqsr
	    output += ' -knownSites $GATKbundle/dbsnp_138.b37.vcf;'
	    output += '\n'

	    output += 'echo -e "-> PrintReads <-"\n'
	    output += 'java -Xmx72g -jar $GATK -T PrintReads -I '+self.realignedbam+' -R '+self.gatkreference+' -BQSR '+self.bqsr+' -o '+self.recalbam+';\n'
	    output += 'rm -v '+self.realignedbam+''+'\n'

	    output += 'echo "Done. $(date) Running on: $(hostname)"\n'
	    output += 'wait\n'
	    output += 'echo "$(date) AllDone"\n'
	    with open(self.realncalscript,'w') as outfile: outfile.write(output)

	def makeScriptExtraQC(self):
	    output = ''
	    output += '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+self.project+'\n'
	    output += '#SBATCH -n 16 -p node'+'\n'
	    output += '#SBATCH -t 2:00:00'+'\n'
	    output += '#SBATCH -J '+self.sampleid+'.ExtraQC'+'\n'
	    output += '#SBATCH -e '+self.path+'/stderr.ExtraQC.txt'+'\n'
	    output += '#SBATCH -o '+self.path+'/stdout.ExtraQC.txt'+'\n'
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'cd '+self.path+'\n'
	    output += 'picard=/proj/b2011168/private/bin/picard-tools-1.63'+'\n'
	    output += 'GATK=/proj/b2011011/fatCellTnt-EB/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar'+'\n'
	    output += 'GATKbundle=/proj/b2011011/polyG_seqcap/GATKBundle/'+'\n'
	    output += 'file='+self.recalbam+'\n'
	    output += 'echo "-----"'+'\n'

	    output += 'java -Xmx72g -jar /proj/b2011168/private/bin/picard-tools-1.63/CalculateHsMetrics.jar '+'BAIT_INTERVALS=/proj/b2011011/polyG_seqcap/targetBeds/agilentBaitIntervalls.bed '+'TARGET_INTERVALS=/proj/b2011011/polyG_seqcap/targetBeds/wholeGenomePolyGsPlusMinus100bp.bed '+'INPUT=$file '+'OUTPUT=$file.hs_metrics.summary.wholeGenome.txt PER_TARGET_COVERAGE=$file.hs_metrics.perTargetCoverage.wholeGenome.txt '+'REFERENCE_SEQUENCE='+self.reference+'\n'

	    output += 'echo'+'\n'
	    output += 'wait'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'
	    with open(self.extraqcscript,'w') as outfile: outfile.write(output)

class collection(object):

	def __init__(self,samples,basepath):
	    import os
	    import time

	    self.samples = samples
	    self.project = self.samples[0].project
	    self.reference = self.samples[0].reference
	    self.bamfiles = ' -I '.join([sample.recalbam for sample in self.samples])

	    self.path = os.path.abspath(basepath)
	    self.scriptpath = self.path+'/script'

	    self.hc		= self.scriptpath +'/haplotypeCaller.sh'
	    self.callindels	= self.scriptpath +'/call.indels.sh'
	    self.callsnps	= self.scriptpath +'/call.snps.sh'
	    self.vqsrindels	= self.scriptpath +'/vqsr.indels.sh'
	    self.vqsrsnps	= self.scriptpath +'/vqsr.snps.sh'

	def createDirs(self):
	    import os
	    try: os.makedirs(self.path)
	    except OSError:pass
	    try: os.makedirs(self.scriptpath)
	    except OSError:pass
	
	def makeScript_haplotypeCaller(self):
	    for part in ['xaa','xab','xac','xad','xae','xaf','xag','xah','xai','xaj','xak','xal','xam','xan','xao','xap','xaq','xar','xas','xat','xau','xav','xaw','xax','xay','xaz','xba','xbb','xbc','xbd','xbe','xbf','xbg','xbh','xbi','xbj','xbk','xbl','xbm','xbn','xbo','xbp','xbq','xbr','xbs','xbt','xbu','xbv','xbw','xbx','xby','xbz','xca','xcb','xcc','xcd','xce','xcf','xcg','xch','xci','xcj','xck','xcl','xcm','xcn','xco','xcp','xcq','xcr','xcs','xct','xcu','xcv','xcw','xcx','xcy','xcz','xda','xdb','xdc','xdd','xde','xdf','xdg','xdh','xdi','xdj','xdk','xdl','xdm','xdn','xdo','xdp','xdq','xdr','xds','xdt','xdu','xdv','xdw','xdx','xdy','xdz','xea','xeb','xec','xed','xee','xef','xeg','xeh','xei','xej','xek','xel','xem','xen','xeo','xep','xeq','xer','xes','xet','xeu','xev','xew','xex','xey','xez','xfa','xfb','xfc','xfd','xfe','xff','xfg','xfh','xfi','xfj','xfk','xfl','xfm','xfn','xfo','xfp','xfq','xfr','xfs','xft','xfu','xfv','xfw','xfx','xfy','xfz','xga','xgb','xgc','xgd','xge','xgf','xgg','xgh','xgi','xgj','xgk','xgl','xgm','xgn','xgo','xgp','xgq','xgr','xgs','xgt','xgu','xgv','xgw','xgx','xgy','xgz','xha','xhb','xhc','xhd','xhe','xhf','xhg','xhh','xhi','xhj','xhk','xhl','xhm','xhn','xho','xhp','xhq','xhr','xhs','xht']:
		output = '#!/bin/bash -l'+'\n'
		output += '#SBATCH -A '+self.project+'\n'
		output += '#SBATCH -n 16 -p node'+'\n'
		output += '#SBATCH -t 5:00:00'+'\n'
		output += '#SBATCH -J '+part+'.haplotypecaller'+'\n'
		output += '#SBATCH -e '+self.path+'/stderr.haplotypecaller.'+part+'.txt'+'\n'
		output += '#SBATCH -o '+self.path+'/stdout.haplotypecaller.'+part+'.txt'+'\n'
		output += '#SBATCH --mail-type=All'+'\n'
		output += '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'
		output += 'echo "$(date) Running on: $(hostname)"'+'\n'
		output += 'echo "$(date) Running on: $(hostname)" >&2'+'\n'
		output += 'cd '+self.path+'\n'
		output += 'picard=/proj/b2011168/private/bin/picard-tools-1.63'+'\n'
		output += 'GATK=/proj/b2011011/fatCellTnt-EB/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar'+'\n'
		output += 'GATKbundle=/proj/b2011011/fatCellTnt-EB/references/GATKbundle/'+'\n'
		output += 'echo "$(date) Running on: $(hostname)"'+'\n'
		
		if part == 'xaa': output += 'python /proj/b2011011/fatCellTnt-EB/script/mappingStatsForExcel.py '+self.samples[0].basepath+'/ > '+self.samples[0].basepath+'/mappingStats.tsv'+'\n'
		output += 'echo "HC" '+'\n'
		
		output += 'java -Xmx120g -jar $GATK -nct 16 '
		output += '-T HaplotypeCaller '
		output += '-R '+self.reference+' '
		output += '-I '+self.bamfiles+' '
		output += '--genotyping_mode DISCOVERY '
		output += '-stand_emit_conf 10 '
		output += '-stand_call_conf 30 '
		output += '-L /proj/b2011011/fatCellTnt-EB/references/truseq_exome_targeted_regions.part'+part+'.bed '
		output += '--dbsnp $GATKbundle/dbsnp_138.b37.vcf '
		output += '-o raw_variants.part.'+part+'.vcf &'+'\n'
		output += 'wait'+'\n'
		output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
		output += 'wait'+'\n'
		output += 'echo "$(date) AllDone"'+'\n'
		output += 'echo "$(date) AllDone" >&2'+'\n\n\n'
		with open(self.hc+part,'w') as outfile: outfile.write(output)

	def makeScript_vqsr_indels(self):
	    output = '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+self.project+'\n'
	    output += '#SBATCH -n 16 -p node'+'\n'
	    output += '#SBATCH -t 72:00:00'+'\n'
	    output += '#SBATCH -J VQSR.INDELs'+'\n'
	    output += '#SBATCH -e '+self.path+'/stderr.vqsr.INDELs.txt'+'\n'
	    output += '#SBATCH -o '+self.path+'/stdout.vqsr.INDELs.txt'+'\n'
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)" >&2'+'\n'
	    output += 'cd '+self.path+'\n'
	    output += 'reference='+self.reference+'\n'
	    output += 'picard=/proj/b2011168/private/bin/picard-tools-1.63'+'\n'
	    output += 'GATK=/proj/b2011011/fatCellTnt-EB/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar'+'\n'
	    output += 'GATKbundle=/proj/b2011011/fatCellTnt-EB/references/GATKbundle/'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'echo "VQSR" '+'\n'
	    output += 'java -Xmx24g -jar $GATK \\'+'\n'
	    output += '   -T VariantRecalibrator \\'+'\n'
	    output += '   -nt 16 \\'+'\n'
	    output += '   -R $reference \\'+'\n'
	    output += '   -input recalibrated_snps_raw_indels.vcf \\'+'\n'
	    output += '   -recalFile indels.raw.recal \\'+'\n'
	    output += '   -tranchesFile indels.raw.tranches \\'+'\n'
	    output += '   -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATKbundle/Mills_and_1000G_gold_standard.indels.b37.vcf \\'+'\n'
	    output += '   -an DP -an FS -an MQRankSum -an ReadPosRankSum \\'+'\n'
	    output += '   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\'+'\n'
	    output += '   --maxGaussians 4 \\'+'\n'
	    output += '   -rscriptFile recalibrate_INDEL_plots.R \\'+'\n'
	    output += '   -mode INDEL'+'\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
	    output += 'wait'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'echo "---> apply recal <---"'+'\n'
	    output += 'java -Xmx72g -jar $GATK \\'+'\n'
	    output += '   -nt 16 \\'+'\n'
	    output += '   -T ApplyRecalibration \\'+'\n'
	    output += '   -R $reference \\'+'\n'
	    output += '   -input recalibrated_snps_raw_indels.vcf \\'+'\n'
	    output += '   -tranchesFile indels.raw.tranches \\'+'\n'
	    output += '   -recalFile indels.raw.recal \\'+'\n'
	    output += '   -o indels.recalibrated.vcf \\'+'\n'
	    output += '   --ts_filter_level 99.0 \\'+'\n'
	    output += '   -mode INDEL'+'\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
	    output += 'wait'+'\n'
	    output += 'echo $file; done'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'
	    with open(self.vqsrindels,'w') as outfile: outfile.write(output)

	def makeScript_vqsr_snps(self):
	    output = '#! /bin/bash -l'+'\n'
	    output += '#SBATCH -A '+self.project+'\n'
	    output += '#SBATCH -n 16 -p node'+'\n'
	    output += '#SBATCH -t 72:00:00'+'\n'
	    output += '#SBATCH -J VQSR.SNPs'+'\n'
	    output += '#SBATCH -e '+self.path+'/stderr.vqsr.SNPs.txt'+'\n'
	    output += '#SBATCH -o '+self.path+'/stdout.vqsr.SNPs.txt'+'\n'
	    output += '#SBATCH --mail-type=All'+'\n'
	    output += '#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)" >&2'+'\n'
	    output += 'cd '+self.path+'\n'
	    output += 'reference='+self.reference+'\n'
	    output += 'picard=/proj/b2011168/private/bin/picard-tools-1.63'+'\n'
	    output += 'GATK=/proj/b2011011/fatCellTnt-EB/bin/GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar'+'\n'
	    output += 'GATKbundle=/proj/b2011011/fatCellTnt-EB/references/GATKbundle/'+'\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output += '\n'
	    output += 'echo "MERGING VCF PARTS" '+'\n'
	    output += 'grep -P "^#" raw_variants.part.xaa.vcf > raw_variants.vcf;'+'\n'
	    output += 'for file in raw_variants.part.x*.vcf;'
	    output += 'do grep -vP "^#" $file >> raw_variants.vcf;'
	    output += 'grep -P "#" raw_variants.vcf > raw_variants.SORTED.vcf;'+'\n'
	    output += 'grep -vP "#" raw_variants.vcf | sort -Vk1,2 >> raw_variants.SORTED.vcf'+'\n'
	    output += 'mv -v raw_variants.SORTED.vcf raw_variants.vcf'+'\n'
	    output += '\n'
	    output += 'echo "---> VQSR <---" '+'\n'
	    output += 'echo "---> VQSR <---" 1>&2'+'\n'
	    output += 'java -Xmx60g -jar $GATK -T VariantRecalibrator '
	    #output += '-nt 16'
	    output += '-R $reference -input raw_variants.vcf'
	    output += ' -resource:hapmap,known=false,training=true,truth=true,prior=15.0'+' $GATKbundle/hapmap_3.3.b37.vcf'
	    output += ' -resource:omni,known=false,training=true,truth=false,prior=12.0' +' $GATKbundle/1000G_omni2.5.b37.vcf'
	    output += ' -resource:dbsnp,known=true,training=false,truth=false,prior=6.0' +' $GATKbundle/dbsnp_138.b37.vcf'
	    output += ' -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -mode SNP'
	    output += ' -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0'
	    output += ' -minNumBad 1000'
	    output += ' -recalFile recalibrate_SNP.recal '
	    output += ' -tranchesFile recalibrate_SNP.tranches '
	    output += ' -rscriptFile recalibrate_SNP_plots.R '
	    output += '\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
	    output += 'wait'+'\n'
	    output += 'echo "-----"'+'\n'
	    output +='\n'
	    output += 'echo "$(date) Running on: $(hostname)"'+'\n'
	    output +='\n'
	    output +='\n'
	    output += 'echo "---> apply recal <---" '+'\n'
	    output += 'echo "---> apply recal <---" 1>&2'+'\n'
	    output +='\n'
	    output += 'java -Xmx60g -jar $GATK '
	    output += '-T ApplyRecalibration -R $reference '
	    output += ' -input raw_variants.vcf '
	    output += ' -recalFile recalibrate_SNP.recal '
	    output += ' -tranchesFile recalibrate_SNP.tranches '
	    output += ' -o recalibrated_snps_raw_indels.vcf '
	    output += ' --ts_filter_level 99.0 -mode SNP'+'\n'
	    output +='\n'
	    output +='python /proj/b2011011/fatCellTnt-EB/script/myHardFilter.py 0 3 recalibrated_snps_raw_indels.vcf > topVariationsBySample.tsv'+'\n'
	    output +='python /proj/b2011011/fatCellTnt-EB/script/myHardFilter.py 0 1 recalibrated_snps_raw_indels.vcf > variationsTable.tsv'+'\n'
	    output +='\n'
	    output += 'echo "Done. $(date) Running on: $(hostname)"'+'\n'
	    output += 'wait'+'\n'
	    output += 'echo "$(date) AllDone"'+'\n'
	    output += 'echo "$(date) AllDone" >&2'+'\n'
	    with open(self.vqsrsnps,'w') as outfile: outfile.write(output)

if __name__ == "__main__":
    main()

######################## END OF SCRIPT ##########################
