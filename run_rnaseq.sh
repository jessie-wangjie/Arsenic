# Configure
export PATH=$PATH:~/data/software/RSEM-1.2.23/:~/data/software/STAR-STAR_2.4.2a/bin/Linux_x86_64/

# parameters
nthread=8

# GTF
# zcat ~/data/annotation/gencode.v19.annotation.gtf.gz | awk -F "\t" '{OFS="\t"; if(!/^#/) print }' > gencode.v19.annotation.tRNA.gtf
# zcat ~/data/annotation/gencode.v19.tRNAs.gtf.gz | awk -F "\t" '{OFS="\t"; if(!/^#/) { $3="exon"; print } }' >> gencode.v19.annotation.tRNA.gtf

# STAR genome
# mkdir STARgenome
# STAR --runMode genomeGenerate --genomeDir STARgenome --genomeFastaFiles /data/db/genome_sequences/hg19/hg19.fa --sjdbGTFfile gencode.v19.annotation.tRNA.gtf --sjdbOverhang 50 --outFileNamePrefix STARgenome --runThreadN $nthread

# RSEM genome
# rsem-prepare-reference --gtf gencode.v19.annotation.tRNA.gtf /data/db/genome_sequences/hg19/hg19.fa RSEMgenome/RSEMref

# STAR
# for i in ../raw_data/Ren.Bladder.Arsenic.*.fastq.gz
for i in ../raw_data/Ren.Fry.*.0.fastq.gz ../raw_data/Ren.Fry.*.1.fastq.gz
# for i in ../raw_data/test.*.fastq.gz
do
	name=`basename $i .fastq.gz`
echo """
	# Align
	mkdir $name 
	STAR --runThreadN $nthread --genomeDir STARgenome --readFilesIn $i --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 --outSAMstrandField intronMotif --outSAMunmapped Within --outWigStrand Unstranded --outFilterMismatchNoverReadLmax 0.04 --outFileNamePrefix $name/ 

	# Signal
	STAR --runMode inputAlignmentsFromBAM --inputBAMfile $name/Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigStrand Unstranded --outFileNamePrefix $name/

	# bigWig
	bedGraphToBigWig $name/Signal.Unique.str1.out.bg STARgenome/chrNameLength.txt $name/$name.star.unique.bw
	bedGraphToBigWig $name/Signal.UniqueMultiple.str1.out.bg STARgenome/chrNameLength.txt $name/$name.star.all.bw

	# sort bam
	samtools sort -n $name/Aligned.toTranscriptome.out.bam -o $name/$name.toTranscriptome.bam -@ $nthread -T $name/$name.toTranscriptome 

	# RSEM
	rsem-calculate-expression --bam --estimate-rspd --calc-ci --seed 12345 -p $nthread --no-bam-output --ci-memory 30000 $name/$name.toTranscriptome.bam RSEMgenome/RSEMref $name/$name\_rsem
	rsem-plot-model $name/$name\_rsem $name/$name\_rsem.pdf

	# clean
	mv $name/Log.final.out $name/$name.star.stats
	rm $name/*.bg $name/Log.out
	
""" > job.$name.sh
	nohup bash job.$name.sh &> job.$name.log &
done

