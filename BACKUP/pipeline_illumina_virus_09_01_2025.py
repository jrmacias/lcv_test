import sys
import subprocess
import argparse
import os
import pandas as pd

def createSubdirectories(output):
	parent_dir=output
	
	primary_dir=os.path.join(parent_dir, "primario")
	secondary_dir=os.path.join(parent_dir, "secundario")
	terciary_dir=os.path.join(parent_dir, "terciario")
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")
	assembly_dir=os.path.join(parent_dir, "assembling")
	annotation_dir=os.path.join(parent_dir, "annotation")
	
	if not os.path.exists(primary_dir):
		os.mkdir(primary_dir)
	if not os.path.exists(secondary_dir):
		os.mkdir(secondary_dir)
	if not os.path.exists(terciary_dir):
		os.mkdir(terciary_dir)
	if not os.path.exists(taxonomy_dir):
		os.mkdir(taxonomy_dir)
	if not os.path.exists(assembly_dir):
		os.mkdir(assembly_dir)
	if not os.path.exists(annotation_dir):
		os.mkdir(annotation_dir)
	
def fastqc(dir, sample, threads, reads1, reads2):
	parent_dir=dir
	out_dir=os.path.join(parent_dir, "primario")
	#out_dir=os.path.join(parent_dir, "trimming")
	cmd="fastqc -o " + out_dir + " -t " + threads + " -f fastq " + reads1 + " " + reads2
	print("Start FastQC")
	print("------------")
	print(cmd)
	os.system(cmd)

def trimmomatic(dir, sample, reads1, reads2, threads="4", phred="phred33", adapters="TruSeq3-PE.fa:2:30:10:2:keepBothReads", leading="30", trailing="30", slidingwindow="10:25", minlen="40"):
	print("Start trimmomatic")
	print("-----------------")

	#java -jar trimmomatic SE -threads 4 -phred33 file_input file_output ILLUMINACLIP:adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	parent_dir=dir
	primary_dir=os.path.join(parent_dir, "primario")
	out_trimmed_reads1_paired=primary_dir+"/"+sample+"_trimmed_1P.fastq.gz"
	out_trimmed_reads2_paired=primary_dir+"/"+sample+"_trimmed_2P.fastq.gz"
	out_trimmed_reads1_unpaired=primary_dir+"/"+sample+"_trimmed_1U.fastq.gz"
	out_trimmed_reads2_unpaired=primary_dir+"/"+sample+"_trimmed_2U.fastq.gz"
	s=" "
	cmd1=" java -jar /usr/share/java/trimmomatic.jar PE -threads " + threads + s + "-phred33 " + reads1 + s + reads2 + s 
	cmd2=out_trimmed_reads1_paired + s + out_trimmed_reads1_unpaired + s + out_trimmed_reads2_paired + s + out_trimmed_reads2_unpaired + s
	#cmd3="ILLUMINACLIP:" + adapters + s + "LEADING:" + leading + " TRAILING:" + trailing + " SLIDINGWINDOW:" + slidingwindow + " MINLEN:" + minlen
	cmd3="ILLUMINACLIP:" + adapters + s + "LEADING:" + leading + " TRAILING:" + trailing + " MINLEN:" + minlen
	cmd=cmd1 + cmd2 + cmd3
	print(cmd)
	os.system(cmd)
	return(out_trimmed_reads1_paired, out_trimmed_reads2_paired, out_trimmed_reads1_unpaired, out_trimmed_reads2_unpaired)


def kraken(dir, sample, paired_reads1, paired_reads2, threads=4):
	print("Start kraken")
	print("------------")
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")
	#/home/vserver1/software/kraken2-2.0.8-beta/kraken2 --db /home/vserver1/software/kraken2-2.0.8-beta/taxoDB --threads 12 --output ./ --gzip-compressed --paired ./AV-01_S32_L001_R1_001.fastq.gz ./AV-01_S32_L001_R2_001.fastq.gz  
	#virus
	
	cmd="/home/admingenomica/software/kraken2/kraken2 --db /labgenomica/kraken_viral_db_v12_01_2024/ --threads " + threads + " --output " + taxonomy_dir + "/" + sample + "_kraken_viral_db.txt --gzip-compressed --paired " + paired_reads1 + " " + paired_reads2 + " --report " + taxonomy_dir + "/" + sample + "_kraken_viral_db_report.txt"
	print(cmd)
	os.system(cmd)

		#standar
	cmd2="/home/admingenomica/software/kraken2/kraken2 --db /labgenomica/kraken_standard_db_v12_01_2024/ --threads " + threads + " --output " + taxonomy_dir + "/" + sample + "_kraken_standard_db.txt --gzip-compressed --paired " + paired_reads1 + " " + paired_reads2 + " --report " + taxonomy_dir + "/" + sample + "_kraken_standard_db_report.txt"

	print(cmd2)
	os.system(cmd2)
	

def mash_screen (dir, db, read1, read2):
	print("Start Mash screen")
	print("------------")
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")

	cmd="mash screen " + db + " " + read1 + " " + read2 + " > " + taxonomy_dir +"/mash_out.tab"
	print(cmd)
	os.system(cmd)
	cmd1="sort -gr " + taxonomy_dir+"/mash_out.tab > " + taxonomy_dir + "/mash_out_ordered.tab"
	print(cmd1)
	os.system(cmd1)

	
def bwa(dir, sample, reads1, reads2, reference, threads=4):
	print("Start bwa")
	print("---------")
	#bwa mem -t 14 ./data/Reference.fna ./data/illumina/Illumina_R1.fastq.gz ./data/illumina/Illumina_R2.fastq.gz > ./mapping/illumina.bwa.sam
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	#secondary_dir=os.path.join(parent_dir, "secundario_ref_pOX_21_1088-10")
	s=" "
	
	cmd1="bwa mem -R \"@RG\\tID:"+ sample +"\\tSM:"+ sample +"\\tPL:ILLUMINA\" -t " + threads + s + reference + s + reads1 + s + reads2 
	cmd2= " | samtools view -Sbh - > "+ secondary_dir + "/" + sample + "_bwa.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	cmd="samtools sort -@" + threads + s + secondary_dir + "/" + sample + "_bwa.bam -o " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	cmd="samtools index " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	##Para amplicones comentar esta orden:
	cmd="java -jar /home/admingenomica/software/picard-2-27-1/picard.jar MarkDuplicates -I " + secondary_dir + "/" + sample + "_bwa_sort.bam -O " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam -M picard_metrics.txt"
	os.system(cmd)
	cmd="samtools index " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam"
	
	
def bowtie(dir, sample, reads1, reads2, reference, threads=4):
	print("Start bowtie")
	print("------------")
	#/home/vserver1/software/bowtie2-2.4.2-linux-x86_64/bowtie2
	#crear indice de referencia: bowtie2-build ./referencia_mpox/ON563414_monkey_pox_virus.fasta ON563414_monkey_pox_virus
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	s=" "
	cmd1="bowtie2 -p " + threads  + s + "--met-file" + s + secondary_dir + "/" + sample + "_bwt_metrics.txt" +  s + "--rg-id " + sample + " --rg SM:" + sample + " --rg PL:ILLUMINA  -x" + s + reference[:-6] + s + "-1" + s + reads1 + s + "-2" + s + reads2
	cmd2= " | samtools view -Sbh - > "+ secondary_dir + "/" + sample + "_bowtie2.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	#samtools sort MS3825.bam -o MS3825.sorted.bam
	cmd="samtools sort -@" + threads + s + secondary_dir + "/" + sample + "_bowtie2.bam -o " + secondary_dir + "/" + sample + "_bowtie2_sort.bam"
	os.system(cmd)
	#samtools index MS3825.sorted.bam
	cmd="samtools index " + secondary_dir + "/" + sample + "_bowtie2_sort.bam"
	os.system(cmd)
	cmd="java -jar /home/admingenomica/software/picard-2-27-1/picard.jar MarkDuplicates -I " + secondary_dir + "/" + sample + "_bowtie2_sort.bam -O " + secondary_dir + "/" + sample + "_bowtie2_sort_markedDup.bam -M picard_metrics.txt"
	os.system(cmd)
	cmd="samtools index " + secondary_dir + "/" + sample + "_bowtie2_sort_markedDup.bam"

def qualimap(dir, sample, alignment_tool):
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	
	#secondary_dir=os.path.join(parent_dir, "secundario_ref_pOX_21_1088-10")
	cmd="qualimap bamqc -sd -bam " + secondary_dir + "/" + sample + "_"+alignment_tool+"_sort_markedDup.bam"
	###para amplicones, sin eliminar duplicados:
	#cmd="qualimap bamqc -sd -bam " + secondary_dir + "/" + sample + "_"+alignment_tool+"_sort.bam"
	###
	print(cmd)
	os.system(cmd)
	

def samtools_coverage(dir, sample, coveragedir, alignment_tool, regions):
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	cmd="samtools coverage " + secondary_dir + "/" + sample + "_" + alignment_tool +"_sort_markedDup.bam > " + secondary_dir + "/" + sample + "_samtools_coverage.tab"
	###para amplicones, sin eliminar duplicados:
	#cmd="samtools coverage " + secondary_dir + "/" + sample + "_" + alignment_tool +"_sort.bam > " + secondary_dir + "/" + sample + "_samtools_coverage.tab"
	###

	print(cmd)
	os.system(cmd)
	dict_cov_by_seg={}
	for line in open(secondary_dir + "/" + sample + "_samtools_coverage.tab"):
		if not line.startswith("#"):
			line=line.rstrip()
			fields=line.split("\t")
			seg=fields[0]
			depth=fields[6]
			cov_hor=fields[5]
			dict_cov_by_seg[seg]=[depth,cov_hor]

	list_depth_by_seg=[]
	list_cov_by_seg=[]

	for chromosome in regions:
		start=regions[chromosome][0]
		end=regions[chromosome][1]
		list_depth_by_seg.append(dict_cov_by_seg[chromosome][0])
		list_cov_by_seg.append(dict_cov_by_seg[chromosome][1])
	'''
	for line in open(regions, "r"):
		line=line.rstrip()
		fields=line.split("\t")
		chromosome=fields[0]
		start=fields[1]
		end=fields[2]
		list_depth_by_seg.append(dict_cov_by_seg[chromosome][0])
		list_cov_by_seg.append(dict_cov_by_seg[chromosome][1])
	'''

	if os.path.isfile(coveragedir+ "/coverage_by_segments.txt"):	
		file2=open(coveragedir + "/coverage_by_segments.txt", "a")
		file2.write(sample + "\t" + "\t".join(list_depth_by_seg)+ "\t" + "\t".join(list_cov_by_seg)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:	
		file2=open(coveragedir + "/coverage_by_segments.txt", "a")
		file2.write("sample\tdepth_1\tdepth_2\tdepth_3\tdepth_4\tdepth_5\tdepth_6\tdepth_7\tdepth_8\tdepth_9\tdepth_10\tcov_1\tcov_2\tcov_3\tcov_4\tcov_5\tcov_6\tcov_7\tcov_8\tcov_9\tcov_10\n")
		file2.write(sample + "\t" + "\t".join(list_depth_by_seg)+ "\t" + "\t".join(list_cov_by_seg)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")


def coverage_alignment(dir, sample, coveragedir, alignment_tool, regions):
	print("Start merging coverage alignment results")
	print("-----------------------------------------")
	parent_dir=dir
	secundario_dir=os.path.join(parent_dir, "secundario")
	ref_length=0
	ass_length=0
	n50=0
	hay_dup="no"
	

	list_depth_by_seg=[]
	for chromosome in regions:
		start=regions[chromosome][0]
		end=regions[chromosome][1]
		list_depth_by_seg.append([chromosome,0])

	'''	
	for line in open(regions, "r"):
		line=line.rstrip()
		fields=line.split("\t")
		chromosome=fields[0]
		start=fields[1]
		end=fields[2]
		list_depth_by_seg.append([chromosome,0])
	'''
		
	###para amplicones, sin eliminar duplicados:
	#for line in open(secundario_dir+"/"+sample+"_"+alignment_tool+"_sort_stats/genome_results.txt", "r"):
	###
	for line in open(secundario_dir+"/"+sample+"_"+alignment_tool+"_sort_markedDup_stats/genome_results.txt", "r"):
		line=line.rstrip()
		if "number of reads = " in line:
			fields_list=line.split(" ")
			nreads=fields_list[-1]
		elif "number of mapped reads = " in line:
			fields_list=line.split(" ")
			nmappedreads=fields_list[-2]
			percentage_mapped_reads=fields_list[-1][1:-1]
		elif "number of duplicated reads (flagged) = " in line:
			fields_list=line.split(" ")
			dupreads=fields_list[-1]
			hay_dup="yes"
		elif "mean coverageData = " in line:
			fields_list=line.split(" ")
			average_depth=fields_list[-1]
		elif "of reference with a coverageData >= 20X" in line:
			fields_list=line.split(" ")
			coverage=fields_list[-8]

		for item in list_depth_by_seg:
			chromosome=item[0]
			if chromosome in line:
				fields_list=line.split("\t")
				depth=fields_list[4][0:6]
				item[1]=depth
		'''	
	#ESTOO ES LO NUEVO
	cmd="samtools depth -a " + secondary_dir + "/" + sample + "_" + alignment_tool +"_sort.bam "  + " | awk \'{sum+=$3} END { printf \"%.2f\\n\", sum/NR}\' > " + secondary_dir + "/" + sample + "samtools_depth_result.tsv"
	print(cmd)
	os.system(cmd)
	for line in open(assembling_dir + "/" + sample + "samtools_depth_result.tsv", "r"):
		samtools_depth=line.rstrip()
	
	cmd="samtools depth -a " + secondary_dir + "/" + sample + "_" + alignment_tool +"_sort.bam " + " | awk \'{if ($3 >= 20) {sum+=1}} END { printf \"%.2f\\n\", (sum * 100)/NR}\' > " + assembling_dir + "/" + sample + "samtools_cobertura_horizontal_result.tsv"
	print(cmd)
	os.system(cmd)
	for line in open(assembling_dir + "/" + sample + "samtools_cobertura_horizontal_result.tsv", "r"):
		samtools_coverage=line.rstrip()
	#HASTA AQUIIII
	'''
	
	if hay_dup=="no":
		dupreads="0"	
	if os.path.isfile(coveragedir+ "/alignment_coverage_stats.txt"):	
		file2=open(coveragedir + "/alignment_coverage_stats.txt", "a")
		#file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(dupreads) + "\t" + str(coverage) + "\t" + str(depth)+"\n")#"\t" +str(depth_1) +"\t"+str(depth_2)+"\t"+str(depth_3)+"\t"+str(depth_4)+"\t"+str(depth_5)+"\t"+str(depth_6)+"\t"+str(depth_7)+"\t"+str(depth_8)+"\t"+str(depth_9)+"\t"+str(depth_10)+"\n")
		file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(dupreads) + "\t" + str(coverage) + "\t" + str(average_depth)+"\t" + str(list_depth_by_seg[0][1]) +"\t"+str(list_depth_by_seg[1][1])+"\t"+str(list_depth_by_seg[2][1])+"\t"+str(list_depth_by_seg[3][1])+"\t"+str(list_depth_by_seg[4][1])+"\t"+str(list_depth_by_seg[5][1])+"\t"+str(list_depth_by_seg[6][1])+"\t"+str(list_depth_by_seg[7][1])+"\t"+str(list_depth_by_seg[8][1])+"\t"+str(list_depth_by_seg[9][1])+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:
		file2=open(coveragedir + "/alignment_coverage_stats.txt", "w")
		#file2.write("sample\ttotal_reads\tmapped_reads\tpercent_mapped_reads\tdup_reads\thorizontal_coverage\tdepth\n")#\tdepth_1\tdepth_2\tdepth_3\tdepth_4\tdepth_5\tdepth_6\tdepth_7\tdepth_8\tdepth_9\tdepth_10\n")
		file2.write("sample\ttotal_reads\tmapped_reads\tpercent_mapped_reads\tdup_reads\thorizontal_coverage\tdepth\tdepth_1\tdepth_2\tdepth_3\tdepth_4\tdepth_5\tdepth_6\tdepth_7\tdepth_8\tdepth_9\tdepth_10\n")
		file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(dupreads) + "\t" + str(coverage) + "\t" + str(average_depth) + "\t" + str(list_depth_by_seg[0][1]) +"\t"+str(list_depth_by_seg[1][1])+"\t"+str(list_depth_by_seg[2][1])+"\t"+str(list_depth_by_seg[3][1])+"\t"+str(list_depth_by_seg[4][1])+"\t"+str(list_depth_by_seg[5][1])+"\t"+str(list_depth_by_seg[6][1])+"\t"+str(list_depth_by_seg[7][1])+"\t"+str(list_depth_by_seg[8][1])+"\t"+str(list_depth_by_seg[9][1])+"\n")
		#file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(dupreads) + "\t" + str(coverage) + "\t" + str(depth)+"\t"+str(depth_1) +"\t"+str(depth_2)+"\t"+str(depth_3)+"\t"+str(depth_4)+"\t"+str(depth_5)+"\t"+str(depth_6)+"\t"+str(depth_7)+"\t"+str(depth_8)+"\t"+str(depth_9)+"\t"+str(depth_10)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")

def bcftools_varcalling(dir, sample, bam, reference, threads="4"):
	print("Start bcftools_variantcalling")
	print("-----------------------------")
	#bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
	parent_dir=dir
	terciary_dir=os.path.join(parent_dir, "terciario")

	s=" "
	
	#cmd1="bcftools mpileup --min-MQ 50 --min-BQ 30 -I --threads " + threads + " -s " + sample + s + " --annotate INFO/AD -f " + reference + s + bam + s + "| bcftools call -mv --ploidy 1 --threads " + threads + " -Ov -o " + terciary_dir + "/" + sample + "_raw.vcf"
	cmd1="bcftools mpileup --min-MQ 40 --min-BQ 30 -I --threads " + threads + " -s " + sample + s + " --annotate INFO/AD -f " + reference + s + bam + s + "| bcftools call -mv --ploidy 1 --threads " + threads + " -Ov -o " + terciary_dir + "/" + sample + "_raw.vcf"

	print(cmd1)
	os.system(cmd1)
	
	#cmd2="bcftools norm -Ov -f " + reference + " -d all -o " + terciary_dir + "/" + sample + "_norm_raw.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#print(cmd2)
	#os.system(cmd2)
	
	cmd3="bcftools filter -Ov -e 'QUAL<100||PL[0:0]<150||DP<5' --threads " + threads + " -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#cmd3="bcftools filter -Ov -e 'QUAL<30 || DP<10' -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_norm_raw.vcf"
	print(cmd3)
	os.system(cmd3)
	
	#cmd4="bcftools view -O v -V indels -o " + terciary_dir + "/" + sample + "_filter.vcf " + parent_dir + "/" + sample + "_AE017334_q30.vcf"
	#print(cmd4)
	#os.system(cmd4)
	

def generate_consensus_sequence(dir, sample, vcf, reference, bed_file, phages_bed="nan"):
	print("Start bcftools - generate consensus sequence")
	print("--------------------------------------------")
	#bcftools consensus -i -s NA001 -f in.fa in.vcf.gz > out.fa
	parent_dir=dir

	#Sacar las posiciones con profundidad < 10 (para susutituirlas por N)
	secundary_dir=os.path.join(parent_dir, "secundario")
	cmd="bedtools genomecov -bga -ibam " + secundary_dir+ "/" + sample + "_bwa_sort_markedDup.bam | awk '$4<10' > " + bed_file
	###para amplicones, sin eliminar duplicados:
	#cmd="bedtools genomecov -bga -ibam " + secundary_dir+ "/" + sample + "_bwa_sort.bam | awk '$4<10' > " + bed_file
	print(cmd)
	os.system(cmd)
	terciary_dir=os.path.join(parent_dir, "terciario")
	cmd="bgzip -c " + vcf + " > " + terciary_dir + "/" + sample + "_filter.vcf.gz"
	#cmd="bgzip -c " + vcf + " > " + terciary_dir + "/" + sample + "_filter_Nophages.vcf.gz"
	print(cmd)
	os.system(cmd)
	cmd1="tabix -p vcf " + terciary_dir + "/" + sample + "_filter.vcf.gz"
	#cmd1="tabix -p vcf " + terciary_dir + "/" + sample + "_filter_Nophages.vcf.gz"
	print(cmd1)
	os.system(cmd1)

	if phages_bed != "nan":
		cmd2="bcftools consensus -p " + sample + " -m " + phages_bed + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"
		print(cmd2)
		os.system(cmd2)
	else:
		cmd2="bcftools consensus -p " + sample + " -f " + reference + " -m" + bed_file+ " "+ terciary_dir + "/" + sample + "_filter.vcf.gz >" + terciary_dir + "/" + sample + "_with_N_consensus.fasta"
		#cmd2="bcftools consensus -p " + sample + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"
		
		print(cmd2)
		os.system(cmd2)	

def extract_reads(dir, sample, reference, alignment_tool):
	print("Start extracting reads")
	print("----------------------")

	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	s=" "
	
	#cmd1="samtools view -u -F 4 " + secondary_dir + "/" + sample + "_bowtie2_sort_markedDup.bam > " + secondary_dir + "/" + sample + "_nodup_mapped.bam"
	#cmd1="samtools view -u -F 4 " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam > " + secondary_dir + "/" + sample + "_nodup_mapped.bam"
	cmd1="samtools view -u -F 4 " + secondary_dir + "/" + sample + "_" + alignment_tool + "_sort_markedDup.bam >"   + secondary_dir + "/" + sample + "_nodup_mapped.bam"
	###para amplicones, sin eliminar duplicados:
	#cmd1="samtools view -u -F 4 " + secondary_dir + "/" + sample + "_" + alignment_tool + "_sort.bam >"   + secondary_dir + "/" + sample + "_nodup_mapped.bam"
	print(cmd1)
	os.system(cmd1)

	cmd2= "samtools sort -@ 10 -n " + secondary_dir + "/" + sample + "_nodup_mapped.bam" + s + "-o " + secondary_dir + "/" + sample + "_nodup_sort_mapped.bam"
	print(cmd2)
	os.system(cmd2)
	
	cmd3="bamToFastq -i " + secondary_dir + "/" + sample + "_nodup_sort_mapped.bam" + s + "-fq " + parent_dir + "/" + sample + "_ref_R1.fastq -fq2 " + parent_dir  + "/" + sample + "_ref_R2.fastq"
	print(cmd3)
	os.system(cmd3)
	
	cmd4="gzip " + parent_dir + "/" + sample + "_ref_*.fastq"
	print(cmd4)
	os.system(cmd4)


def spades(dir, sample, paired_reads1, paired_reads2, threads=20):
	print("Start SPADES")
	print("------------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	#para correr plasmidspades
	#assembling_dir=os.path.join(parent_dir, "plasmidspades")
	cmd1="/home/admingenomica/software/SPAdes-3.15.4-Linux/bin/spades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2
	
	#para correr plasmidspades:
	#cmd1="/home/vserver1/software/SPAdes-3.14.1/bin/plasmidspades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2 + " --s1 " + unpaired_reads1
	#cmd2=" -o " + assembling_dir

	#cmd2=" --cov-cutoff auto" + " -o " + assembling_dir + " --careful --isolate"
	cmd2=" --cov-cutoff auto" + " -o " + assembling_dir + " --careful"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)


def quast(dir, sample, scaffold, paired_reads1, paired_reads2, reference="nan", threads=4):
	print("Start QUAST")
	print("-----------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	
	s=" "
	
	cmd="bwa index " + scaffold
	print(cmd)
	os.system(cmd)

	cmd1="bwa mem -R \"@RG\\tID:"+ sample +"\\tSM:"+ sample +"\\tPL:ILLUMINA\" -t " + threads + s + scaffold + s + paired_reads1 + s + paired_reads2 
	cmd2= " | samtools view -Sbh - > "+ assembling_dir + "/" + sample + "_assembly_bwa.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	cmd="samtools sort -@" + threads + s + assembling_dir + "/" + sample + "_assembly_bwa.bam -o " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam"
	os.system(cmd)
	cmd="samtools index " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam"
	os.system(cmd)
	
	cmd="pilon -Xmx200G --genome " + scaffold + " --frags " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam --output " + sample + " --outdir " + assembling_dir
	print(cmd)
	os.system(cmd)

	cmd="python /labgenomica/scripts/cambiar_contigs_names.py " + assembling_dir + "/" + sample + ".fasta"
	print(cmd)
	os.system(cmd)
	
	if not reference == "nan":
		cmd="quast.py " + " -o " + assembling_dir + " -t " +threads + " -r " + reference + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " " + assembling_dir + "/" + sample + ".fasta"
		#cmd="quast.py " + " -o " + assembling_dir + " -t " +threads + " -r " + reference + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " " + assembling_dir + "/" + sample + ".fasta"
		#cmd="/home/vserver1/software/quast-5.0.2/quast.py " + " -o " + assembling_dir + " -t " +threads + " -r " + reference + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " --single " + unpaired_reads1 + " " + scaffold
	else:
		cmd="quast.py " + " -o " + assembling_dir + " -t " +threads + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " " + assembling_dir + "/" + sample + ".fasta"
		#cmd="quast.py " + " -o " + assembling_dir + " -t " +threads + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " --single " + " " + assembling_dir + "/" + sample + ".fasta"
	print(cmd)
	os.system(cmd)

def coverage_assembly(dir, sample, coveragedir):
	print("Start merging coverage assemblies results")
	print("-----------------------------------------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	#assembling_dir=os.path.join(parent_dir, "hybrid_assembling")
	ref_length=0
	ass_length=0
	n50=0

	for line in open(assembling_dir+"/report.tsv", "r"):
		field_list=line.split("\t")
		field=field_list[0]
		value=field_list[1]
		if field=="Reference length":
			ref_length=int(value)
		elif field=="Total length":
			ass_length=int(value)
		elif field=="# contigs (>= 0 bp)":
			n_contigs_total=int(value)
		elif field=="# contigs":
			n_contigs=int(value)		
		elif field=="# contigs (>= 10000 bp)":
			n_contigs_10000=int(value)
		elif field=="# contigs (>= 25000 bp)":
			n_contigs_25000=int(value)	
		elif field=="# contigs (>= 50000 bp)":
			n_contigs_50000=int(value)
		elif field=="Largest contig":
			largest_contig=int(value)
		elif field=="GC (%)":
			gc_content=float(value)
		elif field=="N50":
			n50=int(value)
		elif field=="L50":
			l50=int(value)
		elif field=="# total reads":
			n_total_reads=int(value)
		elif field=="Mapped (%)":
			percentage_mapped_reads_with_assembly=float(value)
		elif field=="Avg. coverage depth":
			assembly_depth=int(value)

	'''		
	cmd="samtools depth -a " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam" + " | awk \'{sum+=$3} END { printf \"%.2f\\n\", sum/NR}\' > " + assembling_dir + "/" + sample + "samtools_depth_result.tsv"
	print(cmd)
	os.system(cmd)
	
	for line in open(assembling_dir + "/" + sample + "samtools_depth_result.tsv", "r"):
		samtools_depth=line.rstrip()
	cmd="samtools depth -a " + assembling_dir + "/" + sample +"_assembly_bwa_sort.bam" + " | awk \'{if ($3 >= 20) {sum+=1}} END { printf \"%.2f\\n\", (sum * 100)/NR}\' > " + assembling_dir + "/" + sample + "samtools_cobertura_horizontal_result.tsv"
	print(cmd)
	os.system(cmd)
	
	for line in open(assembling_dir + "/" + sample + "samtools_cobertura_horizontal_result.tsv", "r"):
		samtools_coverage=line.rstrip()
	'''
	if os.path.isfile(coveragedir+ "/assembling_coverage_analysis_results.txt"):	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "a")
		#estaafile2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\t" + str(samtools_coverage)+ "\t" + str(samtools_depth)+"\n")
		file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) +"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "w")
		#file2.write("sample\ttotal_contigs\tn_contigs>=500\tcontigs>=10000bp\tcontigs>=25000bp\tcontigs>=50000bp\tlargest_contig\tassembly_length\tN50\tL50\tGC\tn_total_reads\t%_mapped_reads\tquast_depth\tsamtools_horizontal_cov\tsamtools_depth\n")
		file2.write("sample\ttotal_contigs\tn_contigs>=500\tcontigs>=10000bp\tcontigs>=25000bp\tcontigs>=50000bp\tlargest_contig\tassembly_length\tN50\tL50\tGC\tn_total_reads\t%_mapped_reads\tquast_depth\n")
		file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\n")
		#file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\t" + str(samtools_coverage)+ "\t" + str(samtools_depth)+"\n")
		
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")


def prokka(dir, sample, consensus, threads=4):
	#do /home/vserver1/software/prokka/bin/prokka --outdir ./$c/annotation --force --prefix $c"_ann" ./$c/annotation/$c"_assembling.result.fasta"
	print("Start Prokka")
	print("------------")
	
	annotation_dir="/home/admingenomica/analisis/temp/"+sample
	
	cmd= "/home/admingenomica/software/prokka/bin/prokka --kingdom Viruses --proteins /home/admingenomica/software/resfinder/db_resfinder/resfinder_all_prot_db.fasta --outdir " + annotation_dir + " --centre X --compliant --gcode 11 --cpus " + threads + " --prefix " + sample+ " " + consensus

	print(cmd)
	os.system(cmd)	
	cmd2="mv /home/admingenomica/analisis/temp/"+sample+"/* "+dir+"/annotation/"
	
	print(cmd2)
	os.system(cmd2)


def definir_regiones_referencia(referencia):
	dict_regions={}
	for line in open(referencia, "r"):
		line=line.rstrip()
		if line.startswith(">"):
			cabecero=line.split(" ")[0][1:]
			dict_regions[cabecero]=[]
		else:
			dict_regions[cabecero]+=line

	for seg in dict_regions:
		length=len(dict_regions[seg])
		coordenadas=[1,length]
		dict_regions[seg]=coordenadas
		#print(seg, " ", dict_regions[seg][0], " ", dict_regions[seg][1], "\n")



	return dict_regions

'''
def cov_depth_by_segment(dir, sample,mergeresdir, region):
	
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario_best_ref_2")
	
	cmd="samtools index " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam"
	print(cmd)
	os.system(cmd)
	#cov horizontal (>20X) por segmento
	cmd1= "echo " + sample + " >> "+ mergeresdir+"horizontal_coverage_by_segment_samtools_depth_mayor20"
	os.system(cmd1)
	cmd2="while read -r chrom start end; do echo $chrom; samtools depth -m 0 -a "+ secondary_dir+"/"+sample+"_bwa_sort.bam -r $chrom | awk '{ if ($3 >= 20) {sum+=1}} END { printf '%.2f\n', (sum * 100)/NR}'; done <" + region +">>" + mergeresdir + "/horizontal_coverage_by_segment_samtools_depth_mayor20"
	os.system(cmd2)
	#cov horizontal (>10X) por segmento
	cmd1= "echo " + sample + " >> "+ mergeresdir+"horizontal_coverage_by_segment_samtools_depth_mayor10"
	os.system(cmd1)
	cmd2="while read -r chrom start end; do echo $chrom; samtools depth -m 0 -a "+ secondary_dir+"/"+sample+"_bwa_sort.bam -r $chrom | awk '{ if ($3 >= 10) {sum+=1}} END { printf '%.2f\n', (sum * 100)/NR}'; done <" + region +">>" + mergeresdir + "/horizontal_coverage_by_segment_samtools_depth_mayor10"
	os.system(cmd2)

'''




##################################################################################################################################################################
##################################################################################################################################################################

def main():
	#pasar el fichero conel directorio y el nombre de las muestras y recorrer cada muestras antes de hacer el analisis filogenetico.

	#Initiate the parser

	help_text="Pipeline for NGS analysis of illumina short reads from microrganism samples. This includes the following steps: raw reads quality analysis and trimming, assembling, genome annotation, mapping with reference genome, variant calling, identification of arm genes..."

	parser = argparse.ArgumentParser(description=help_text)

	#Add arguments
	#argumento sin valor
	#parser.add_argument("-V", "--version", help="show program version", action="store_true")

	required=parser.add_argument_group('required arguments')

	required.add_argument("--sample", "-s", help="sample name or id", required=True)
	required.add_argument("--output", "-o", help="output directory", required=True)
	parser.add_argument("--threads", "-t", help="number of threads")
	required.add_argument("--reads1", "-r1", help="fastq file with reads1 (paired-end)", required=True)
	required.add_argument("--reads2", "-r2", help="fastq file with reads2 (paired-end)", required=True)
	required.add_argument("--reference", "-Ref", help="reference genome (fasta file)") #En caso de no tener referencia, comentar todos los pasos del alineamiento con referencia y variant calling y ejecutar quast sin referencia (comentar linea correspondiente)
	required.add_argument("--mergeresdir", "-md", help="output directory for merged results", required=True)
	required.add_argument("--alignment_tool", "-alt", help="tool for aligning sequencing reads to long reference sequences: bwa or bowtie2")

	args=parser.parse_args()

	if args.sample:
		print("sample = %s" % args.sample)
		sample=args.sample

	if args.output:
		print("output dir = %s" % args.output)
		output=args.output

	if args.threads:
		print("threads = %s" % args.threads)
		threads=args.threads

	if args.reads1:
		print("reads1 = %s" % args.reads1)
		reads1=args.reads1

	if args.reads2:
		print("reads2 = %s" % args.reads2)
		reads2=args.reads2

	if args.mergeresdir:
		print("mergeresdir = %s" % args.mergeresdir)
		mergeresdir=args.mergeresdir

	if args.alignment_tool:
		print("alignment_tool=%s" % args.alignment_tool)
		alignment_tool=args.alignment_tool

	if args.reference:
		print("reference = %s" % args.reference)
		reference=args.reference
		regions=definir_regiones_referencia(reference)

	
	
	#create subdirectories
	createSubdirectories(output)
	'''
	#run fastqc
	fastqc(output, sample, threads, reads1, reads2)
	
	#run trimmomatic
	trimmed_paired_reads1, trimmed_paired_reads2, trimmed_unpaired_reads1, trimmed_unpaired_reads2=trimmomatic(output, sample, reads1, reads2, threads)
	
	#run fastqc
	fastqc(output, sample, threads, trimmed_paired_reads1, trimmed_paired_reads2)
	#fastqc(output, sample, threads, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz")
	
	#run kraken
	kraken(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", threads)
	
	#mash_screen (output, "/home/admingenomica/software/mash_screen_db/RefSeq88n.msh", reads1, reads2)
	mash_screen (output, "/home/admingenomica/software/mash_screen_db/refseq.genomes.k21s1000.msh", output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	
	if args.reference:
		if alignment_tool == "bwa":
			
			#MAPPING
			#run bwa
			#bwa(output, sample, trimmed_paired_reads1, trimmed_paired_reads2, reference, threads)
			bwa(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", reference, threads)
			#bwa(output, sample, reads1, reads2, reference, threads)
			#bwa(output, sample, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz", reference, threads)
			#bwa(output, sample, output+sample+"_ref_R1.fastq.gz", output+sample+"_ref_R2.fastq.gz", reference, threads)
			
			#run qualimap
			qualimap(output, sample, alignment_tool)
			
			#samtools coverage (for references with diferent segments or chromosomes)
			samtools_coverage(output, sample, mergeresdir, alignment_tool, regions)
			
			
			#run coverage_alignment
			coverage_alignment(output, sample, mergeresdir, alignment_tool, regions)		
			
			
			#run variantcalling
			bcftools_varcalling(output, sample, output + "/secundario/" + sample + "_bwa_sort_markedDup.bam", reference, threads)
			#bcftools_varcalling(output, sample, output + "/secundario/" + sample + "_bwa_sort.bam", reference, threads)
			#bcftools_varcalling(output, sample, output + "/secundario_ref_pOX_21_1088-10/" + sample + "_bwa_sort_markedDup.bam", reference, threads)
			
			generate_consensus_sequence(output, sample, output + "/terciario/" + sample + "_filter.vcf", reference, output + "/secundario/less_10_rango.bed")
			
			#cov_depth_by_segment(output, sample, mergeresdir, output + "/"+sample+ "/best_ref_2/regions.txt")
		if alignment_tool == "bowtie2":
			#bowtie(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", reference, threads)

			#run qualimap
			#qualimap(output, sample, alignment_tool)
			
			#samtools coverage (for references with diferent segments or chromosomes)
			samtools_coverage(output, sample, mergeresdir, alignment_tool)
			
			#run coverage_alignment
			coverage_alignment(output, sample, mergeresdir, alignment_tool)		
	
	
	#Anotacion secuencia consenso
	#prokka(output, sample, output+"/terciario_best_ref/"+sample+"_with_N_consensus.fasta", threads)
	prokka(output, sample, output+"/terciario/"+sample+"_with_N_consensus.fasta", threads)
	'''
	extract_reads(output, sample, reference, alignment_tool)
	
	#run spades
	spades(output, sample, output+"/"+sample+"_ref_R1.fastq.gz", output+"/"+sample+"_ref_R2.fastq.gz", threads)
	#spades(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", threads)
	#spades(output, sample, trimmed_paired_reads1, trimmed_paired_reads2,threads)
	#spades(output, sample, reads1, reads2, threads)

	#run 
	#con referencia
	
	if args.reference:
		quast(output, sample, output+"/assembling/scaffolds.fasta",  output+"/"+sample+"_ref_R1.fastq.gz", output+"/"+sample+"_ref_R2.fastq.gz", reference, threads)
		#quast(output, sample, output+"/assembling/scaffolds.fasta",  output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", reference, threads)
		#quast(output, sample, output+"/assembling/scaffolds.fasta",  reads1,  reads2, reference, threads)
	else:
		#sin referencia
		quast(output, sample, output+"/assembling/scaffolds.fasta", output+"/"+sample+"_ref_R1.fastq.gz", output+"/"+sample+"_ref_R2.fastq.gz", threads=threads)
		#quast(output, sample, output+"/assembling_metaspades/scaffolds.fasta", output +"/"+sample+"_1P.fastq.gz",  output +"/"+sample+"_2P.fastq.gz", threads=threads)
		#quast(output, sample, output+"/assembling/scaffolds.fasta", reads1, reads2, threads=threads)
	
	
	#run coverage_assembly
	coverage_assembly(output, sample, mergeresdir)
	
	
	
if __name__ == "__main__":

    main()