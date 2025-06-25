import sys
import subprocess
import argparse
import os
import gzip
import pandas as pd


#while read line; do python ../../pipeline_illumina.py -s $line -o ./$line -t 12 -r1 ./$line/*R1*fastq.gz -r2 ./$line/*R2*fastq.gz 
#-R ../campylobacter_reference/Campylobacter_jejuni_NCTC11168.fasta; done < ../id_cepas_analizar
###############################################FUNCTIONS######################################################

def createSubdirectories(output):
	parent_dir=output
	
	primary_dir=os.path.join(parent_dir, "primario")
	secundary_dir=os.path.join(parent_dir, "secundario")
	terciary_dir=os.path.join(parent_dir, "terciario")
	resfinder_dir=os.path.join(parent_dir, "resfinder")
	#card_dir=os.path.join(parent_dir, "card")
	plasmidfinder_dir=os.path.join(parent_dir, "plasmidfinder")
	#integronfinder_dir=os.path.join(parent_dir, "integronfinder")
	recycler_dir=os.path.join(parent_dir, "recycler")
	abricate_dir=os.path.join(parent_dir, "virulence_factors")
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")
	sistr_dir=os.path.join(parent_dir, "sistr")
	serotypefinder_dir=os.path.join(parent_dir, "serotypefinder")
	fimtype_dir=os.path.join(parent_dir, "fimtype")
	seqsero2_dir=os.path.join(parent_dir, "seqsero2")
	prokka_dir=os.path.join(parent_dir, "annotation")
	virulencefinder_dir=os.path.join(parent_dir, "virulencefinder")
	mlst_dir=os.path.join(parent_dir, "mlst")
	lissero_dir=os.path.join(parent_dir, "lissero")
	checkm_dir=os.path.join(parent_dir, "checkm")
	if not os.path.exists(primary_dir):
		os.mkdir(primary_dir)
	if not os.path.exists(secundary_dir):
		os.mkdir(secundary_dir)
	if not os.path.exists(terciary_dir):
		os.mkdir(terciary_dir)
	if not os.path.exists(resfinder_dir):
		os.mkdir(resfinder_dir)
	if not os.path.exists(plasmidfinder_dir):
		os.mkdir(plasmidfinder_dir)
	#if not os.path.exists(integronfinder_dir):
	#	os.mkdir(integronfinder_dir)
	if not os.path.exists(recycler_dir):
		os.mkdir(recycler_dir)
	if not os.path.exists(abricate_dir):
		os.mkdir(abricate_dir)
	#if not os.path.exists(card_dir):
	#	os.mkdir(card_dir)
	if not os.path.exists(taxonomy_dir):
		os.mkdir(taxonomy_dir)
	if not os.path.exists(sistr_dir):
		os.mkdir(sistr_dir)
	if not os.path.exists(serotypefinder_dir):
		os.mkdir(serotypefinder_dir)
	if not os.path.exists(fimtype_dir):
		os.mkdir(fimtype_dir)
	if not os.path.exists(seqsero2_dir):
		os.mkdir(seqsero2_dir)
	if not os.path.exists(prokka_dir):
		os.mkdir(prokka_dir)
	if not os.path.exists(virulencefinder_dir):
		os.mkdir(virulencefinder_dir)
	if not os.path.exists(mlst_dir):
		os.mkdir(mlst_dir)
	if not os.path.exists(lissero_dir):
		os.mkdir(lissero_dir)
	if not os.path.exists(checkm_dir):
		os.mkdir(checkm_dir)


def fastqc(dir, sample, threads, reads1, reads2):
	parent_dir=dir
	out_dir=os.path.join(parent_dir, "primario")
	#out_dir=os.path.join(parent_dir, "trimming")
	cmd="/home/admingenomica/software/FastQC-0.12.1/fastqc -o " + out_dir + " -t " + threads + " -f fastq " + reads1 + " " + reads2
	print("Start FastQC")
	print("------------")
	print(cmd)
	os.system(cmd)


#def trimmomatic(dir, sample, reads1, reads2, threads="4", phred="phred33", adapters="TruSeq3-PE.fa:2:30:10:2:keepBothReads", leading="30", trailing="30", slidingwindow="10:20", minlen="40"):
	#(dir, sample, reads1, reads2, threads="4", phred="phred33", adapters="TruSeq3-PE.fa:2:30:10:2:keepBothReads", leading="25", trailing="25", slidingwindow="10:30", minlen="40"):
def trimmomatic(dir, sample, reads1, reads2, threads="4", phred="phred33", adapters="NexteraPE-PE.fa:2:30:10:2:keepBothReads", leading="30", trailing="30", slidingwindow="5:30", minlen="40"):
	print("Start trimmomatic")
	print("-----------------")

	#java -jar trimmomatic SE -threads 4 -phred33 file_input file_output ILLUMINACLIP:adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	parent_dir=dir
	primary_dir=os.path.join(parent_dir, "primario")
	out_trimmed_reads1_paired=primary_dir+"/"+sample+"_trimmed_1P.fastq.gz"
	out_trimmed_reads2_paired=primary_dir+"/"+sample+"_trimmed_2P.fastq.gz"
	out_trimmed_reads1_unpaired=primary_dir+"/"+sample+"_trimmed_1U.fastq.gz"
	out_trimmed_reads2_unpaired=primary_dir+"/"+sample+"_trimmed_2U.fastq.gz"
	#out_trimmed_read=primary_dir+"/"+sample+"_trimmed.fastq.gz"
	s=" "
	cmd1=" /usr/lib/jvm/java-11-openjdk-amd64/bin/java -jar /usr/share/java/trimmomatic.jar PE -threads " + threads + s + "-phred33 " + reads1 + s + reads2 + s 
	#cmd1 = " java -jar /usr/share/java/trimmomatic.jar SE -threads " + threads + s + "-phred33 " + reads1 + s 
	cmd2=out_trimmed_reads1_paired + s + out_trimmed_reads1_unpaired + s + out_trimmed_reads2_paired + s + out_trimmed_reads2_unpaired + s
	#cmd2= out_trimmed_read + s
	#cmd3="ILLUMINACLIP:" + adapters + s + "LEADING:" + leading + " TRAILING:" + trailing + " SLIDINGWINDOW:" + slidingwindow + " MINLEN:" + minlen
	cmd3="ILLUMINACLIP:" + adapters + s + "LEADING:" + leading + " TRAILING:" + trailing + " MINLEN:" + minlen
	#cmd3="ILLUMINACLIP:TruSeq3-SE.fa:2:30:10" + s + "LEADING:" + leading + " TRAILING:" + trailing + " MINLEN:" + minlen
	cmd=cmd1 + cmd2 + cmd3
	print(cmd)
	os.system(cmd)
	return(out_trimmed_reads1_paired, out_trimmed_reads2_paired, out_trimmed_reads1_unpaired, out_trimmed_reads2_unpaired)


def q30(dir, sample, reads1, reads2, out_trimmed_reads1_paired, out_trimmed_reads2_paired):
	print("Calculating Q30 bases")
	print("---------------------")
	parent_dir=dir
	out_dir=os.path.join(parent_dir, "primario")

	#awk 'NR%4==0' your_file.fastq | awk '{for (i=1; i<=length($0); i++) {total++; if (substr($0,i,1) >= "!") q30_count++;}} END {printf "%.2f%%\n", (q30_count/total)*100}'
																																				
	cmd= "zcat " + reads1 + " " + reads2 + " | awk 'NR%4==0' | awk '{for (i=1; i<=length($0); i++) {total++; if (substr($0,i,1) >= \"!\") q30_count++;}} END {printf \"%.2f\\n\", (q30_count/total)*100}'"
	
	print(cmd)
	percentage_q30_raw_reads = subprocess.check_output(cmd, shell=True)
	percentage_q30_raw_reads_decode = percentage_q30_raw_reads.decode("utf-8").rstrip()

	cmd= "zcat " + out_trimmed_reads1_paired + " " + out_trimmed_reads2_paired + " | awk 'NR%4==0' | awk '{for (i=1; i<=length($0); i++) {total++; if (substr($0,i,1) >= \"!\") q30_count++;}} END {printf \"%.2f\\n\", (q30_count/total)*100}'"
	print(cmd)
	percentage_q30_trimmed_reads = subprocess.check_output(cmd, shell=True)
	percentage_q30_trimmed_reads_decode = percentage_q30_trimmed_reads.decode("utf-8").rstrip()

	q30_file=open(out_dir+"/"+ sample + "_q30_results.tab", "a")
	q30_file.write("%Q30_bases_raw_reads = " + percentage_q30_raw_reads_decode + "%\n%Q30_bases_trimmed_reads = " + percentage_q30_trimmed_reads_decode + "%")




def calculate_q30_percentage_from_compressed_fastq(dir, sample, reads1, reads2, out_trimmed_reads1_paired, out_trimmed_reads2_paired):
	print("Calculating Q30 bases")
	print("---------------------")
	parent_dir=dir
	out_dir=os.path.join(parent_dir, "primario")
	total_bases_raw_reads = 0
	q30_or_higher_bases_raw_reads = 0

	with gzip.open(reads1, 'rt') as f:
	    	for line_number, line in enumerate(f, start=1):
    			# Quality scores are in the 4th line of every 4-line FASTQ record
    			if line_number % 4 == 0:
    				for char in line.strip():
    					total_bases_raw_reads += 1
    					if ord(char) - 33 >= 30:
    						q30_or_higher_bases_raw_reads += 1

	with gzip.open(reads2, 'rt') as f:
		for line_number, line in enumerate(f, start=1):
			# Quality scores are in the 4th line of every 4-line FASTQ record
			if line_number % 4 == 0:
				for char in line.strip():
					total_bases_raw_reads += 1
					if ord(char) - 33 >= 30:
						q30_or_higher_bases_raw_reads += 1

	percentage_q30_raw_reads_decode_=(q30_or_higher_bases_raw_reads/total_bases_raw_reads)*100
	percentage_q30_raw_reads_decode="%.2f" % percentage_q30_raw_reads_decode_
	q30_file=open(out_dir+"/"+ sample + "_q30_results.tab", "a")
	q30_file.write("%Q30_bases_raw_reads = " + str(percentage_q30_raw_reads_decode) + "%")

def reads_PF_post_trimmeed(dir, sample, read1, read2, out_trimmed_reads1, out_trimmed_reads2):
	print("Calculating number PF reads post Trimmomatic")
	print("--------------------------------------------")
	parent_dir=dir
	out_dir=os.path.join(parent_dir, "primario")

	count_pre_read1 = 0
	with gzip.open(read1, 'r') as r1:
		for line in r1:
			line=line.rstrip()
			if line.startswith(b'@'):
				count_pre_read1+=1
			else:
				continue
	
	count_pre_read2 = 0
	with gzip.open(read2, 'r') as r2:
		for line in r2:
			line=line.rstrip()
			if line.startswith(b'@'):
				count_pre_read2+=1
			else:
				continue
	
	count_post_read1 = 0
	with gzip.open(out_trimmed_reads1, 'r') as t1:
		for line in t1:
			line=line.rstrip()
			if line.startswith(b'@'):
				count_post_read1+=1
			else:
				continue
	
	count_post_read2 = 0
	with gzip.open(out_trimmed_reads2, 'r') as t2:
		for line in t2:
			line=line.rstrip()
			if line.startswith(b'@'):
				count_post_read2+=1
			else:
				continue
	
	noPF_read1=count_pre_read1-count_post_read1
	noPF_read2=count_pre_read2-count_post_read2

	PF_file=open(out_dir+"/"+ sample + "_PF_trimmed.tab", "a")
	PF_file.write("Number PF R1\tNumber no PF R1\tNumber PF R2\tNumber no PF R2\n" + str(count_post_read1)+"\t"+str(noPF_read1)+"\t"+str(count_post_read2)+"\t"+str(noPF_read2))
				
	#PF_file.write("Number PF R1\t" + str(count_post_read1) +"\tNumber no PF R1\t"+ str(PF_read1) + "\tNumber PF R2\t" + str(count_post_read2) + "\tNumber no PF R2\t" + str(PF_read2))
	

def kraken(dir, sample, paired_reads1, paired_reads2, threads=4):
	print("Start kraken")
	print("------------")
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")
	
	#kraken standard db
	#/home/admingenomica/software/kraken2/kraken2 --db /labgenomica/kraken_viral_v01_12_2024/ --threads 10 --output ./taxonomy_new_viral_kraken_db/AHSV-PRUEBA-V_kraken.txt --gzip-compressed --paired ./primario/AHSV-PRUEBA-V_trimmed_1P.fastq.gz ./primario/AHSV-PRUEBA-V_trimmed_2P.fastq.gz --report ./taxonomy_new_viral_kraken_db/AHSV-PRUEBA-V_kraken_report.txt
	cmd="/home/admingenomica/software/kraken2-2.14/kraken2 --db /labgenomica/kraken_standard_db_v28_12_2024/ --threads " + threads + " --output " + taxonomy_dir + "/" + sample + "_kraken_standard_db.txt --gzip-compressed --paired " + paired_reads1 + " " + paired_reads2 + " --report " + taxonomy_dir + "/" + sample + "_kraken_standard_db_report.txt"
	print(cmd)
	os.system(cmd)
	
	#kraken viral db
	cmd2="/home/admingenomica/software/kraken2-2.14/kraken2 --db /labgenomica/kraken_viral_db_v28_12_2024/ --threads " + threads + " --output " + taxonomy_dir + "/" + sample + "_kraken_viral_db.txt --gzip-compressed --paired " + paired_reads1 + " " + paired_reads2 + " --report " + taxonomy_dir + "/" + sample + "_kraken_viral_db_report.txt"
	print(cmd2)
	os.system(cmd2)
	
	#mininusB db
	cmd="/home/admingenomica/software/kraken2-2.14/kraken2 --db /labgenomica/kraken_minusB_db_8G_v28_12_2024/ --threads " + threads + " --output " + taxonomy_dir + "/" + sample + "_minusB_db.txt --gzip-compressed --paired " + paired_reads1 + " " + paired_reads2 + " --report " + taxonomy_dir + "/" + sample + "_minusB_db_report.txt"
	#print(cmd)
	os.system(cmd)

	#minikraken
	cmd="/home/admingenomica/software/kraken/kraken --db /labgenomica/minikraken_db/minikraken_20171019_8GB --threads " + threads + " --output " + taxonomy_dir + "/" + sample + "_minikraken.txt --fastq-input --gzip-compressed --paired " + paired_reads1 + " " + paired_reads2
	print(cmd)
	os.system(cmd)
	cmd2="/home/admingenomica/software/kraken/kraken-report --db /labgenomica/minikraken_db/minikraken_20171019_8GB " + taxonomy_dir + "/" + sample + "_minikraken.txt > " + taxonomy_dir + "/" + sample+"_minikraken_report.txt"
	print(cmd2)
	os.system(cmd2)
	


def confindr (dir, path_db):
	print("Start Confindr")
	print("------------")
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")	
	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'confindr'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_confindr.py ' + parent_dir + " " + path_db + " " + taxonomy_dir
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)


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



def spades(dir, sample, paired_reads1, paired_reads2, unpaired_reads1, unpaired_reads2, threads=4):#  
	print("Start SPADES")
	print("------------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	#para correr plasmidspades
	#assembling_dir=os.path.join(parent_dir, "plasmidspades")
	cmd1="spades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2 + " -s " + unpaired_reads1 + " -s " + unpaired_reads2
	#cmd1="/home/admingenomica/software/SPAdes-4.1.0-Linux/bin/spades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2
	#para correr plasmidspades:
	#cmd1="/home/vserver1/software/SPAdes-4.1.0-Linux/bin/plasmidspades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2 + " --s1 " + unpaired_reads1
	#cmd2=" -o " + assembling_dir
	cmd2=" --cov-cutoff auto" + " -o " + assembling_dir + " --careful"
	#cmd2=" --cov-cutoff auto" + " -o " + assembling_dir + " --careful --plasmid"
	#cmd2=" --cov-cutoff auto" + " -o " + assembling_dir + " --isolate"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)


def pilon(dir, sample, scaffold, paired_reads1, paired_reads2, unpaired_reads1, unpaired_reads2, threads=4):
	print("Start Pilon")
	print("-----------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	
	s=" "
	
	cmd="bwa index " + scaffold
	print(cmd)
	os.system(cmd)

	cmd1="bwa mem -R \"@RG\\tID:"+ sample +"\\tSM:"+ sample +"\\tPL:ILLUMINA\" -t " + threads + s + scaffold + s + paired_reads1 + s + paired_reads2 
	cmd2= " | /home/admingenomica/software/samtools-1.21/bin/samtools view -Sbh - > "+ assembling_dir + "/" + sample + "_assembly_bwa.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools sort -@" + threads + s + assembling_dir + "/" + sample + "_assembly_bwa.bam -o " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam"
	os.system(cmd)
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools index " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam"
	os.system(cmd)
	
	cmd="pilon -Xmx200G --genome " + scaffold + " --frags " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam --output " + sample + " --outdir " + assembling_dir
	print(cmd)
	os.system(cmd)

	cmd="python /labgenomica/scripts/cambiar_contigs_names.py " + assembling_dir + "/" + sample + ".fasta"
	print(cmd)
	os.system(cmd)


def remove_short_contigs(dir, sample, scaffold, length="500"):
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'bbmap'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_reformat.py ' + sample + " " + scaffold + " " + assembling_dir + " " + length
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)


def checkm(dir, sample, scaffold, threads=4):
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	checkm_dir=os.path.join(parent_dir, "checkm")

	cmd="cp " + scaffold + " " + checkm_dir + "/"
	print(cmd)
	os.system(cmd)

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'checkm'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_checkm.py ' + checkm_dir + " " + threads
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)

	cmd="rm " + checkm_dir + "/" + sample + "_gt500.fasta"
	print(cmd)
	os.system(cmd)

	

def quast(dir, scaffold, paired_reads1, paired_reads2, unpaired_reads1, unpaired_reads2, reference="nan", threads=4):  #
	print("Start QUAST")
	print("-----------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	
	s=" "
	if not reference == "nan":
		parent_dir=dir
		# Replace 'my_environment' with the name of your Conda environment
		environment_name = 'quast'
		# Command to activate the Conda environment
		activate_command = f'conda run -n {environment_name}'
		# Command to execute your Python script
		python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_quast_bacterias.py ' + scaffold + " "+ paired_reads1 + " "+ paired_reads2 + " " + unpaired_reads1 + " " + unpaired_reads2 + " " + assembling_dir+ " " + threads + " " + reference
		# Construct the final command to run
		final_command = f'{activate_command} python {python_script}'
		# Execute the command
		subprocess.call(final_command, shell=True)
	else:
		parent_dir=dir
		# Replace 'my_environment' with the name of your Conda environment
		environment_name = 'quast'
		# Command to activate the Conda environment
		activate_command = f'conda run -n {environment_name}'
		# Command to execute your Python script
		python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_quast_bacterias.py ' + scaffold + " "+ paired_reads1 + " "+ paired_reads2 + " " + unpaired_reads1 + " " + unpaired_reads2 + " " + assembling_dir+ " " + threads
		# Construct the final command to run
		final_command = f'{activate_command} python {python_script}'
		# Execute the command
		subprocess.call(final_command, shell=True)




def coverage_assembly(dir, sample, coveragedir, out_trimmed_reads1_paired, out_trimmed_reads2_paired):
	print("Start merging coverage assemblies results")
	print("-----------------------------------------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	primary_dir=os.path.join(parent_dir, "primario")

	#zcat ./primario/2024-000676-1A_trimmed_1P.fastq.gz ./primario/2024-000676-1A_trimmed_2P.fastq.gz | awk 'NR%4 == 2 {sum+=length($0); n+=1} END { print "Average = "sum/n}'
	cmd="zcat " + out_trimmed_reads1_paired + " " + out_trimmed_reads2_paired + " | awk 'NR%4 == 2 {sum+=length($0); n+=1} END { print sum/n}'"
	average_length_trimmed_reads = subprocess.check_output(cmd, shell=True)
	average_length_trimmed_reads_decode = average_length_trimmed_reads.decode("utf-8").rstrip()

	for line in open(primary_dir+"/"+sample+"_q30_results.tab", "r"):
		field_list=line.split(" = ")
		q30=field_list[1]

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

			
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools depth -m 0 -a " + assembling_dir + "/" + sample + "_assembly_bwa_sort.bam" + " | awk \'{sum+=$3} END { printf \"%.2f\\n\", sum/NR}\' > " + assembling_dir + "/" + sample + "samtools_depth_result.tsv"
	print(cmd)
	os.system(cmd)
	for line in open(assembling_dir + "/" + sample + "samtools_depth_result.tsv", "r"):
		samtools_depth=line.rstrip()

	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools depth -m 0 -a " + assembling_dir + "/" + sample +"_assembly_bwa_sort.bam" + " | awk \'{if ($3 >= 20) {sum+=1}} END { printf \"%.2f\\n\", (sum * 100)/NR}\' > " + assembling_dir + "/" + sample + "samtools_cobertura_horizontal_result.tsv"
	print(cmd)
	os.system(cmd)
	for line in open(assembling_dir + "/" + sample + "samtools_cobertura_horizontal_result.tsv", "r"):
		samtools_coverage=line.rstrip()
	
	if os.path.isfile(coveragedir+ "/assembling_coverage_analysis_results.txt"):	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "a")
		file2.write(sample + "\t" + str(q30) + "\t" + str(n_total_reads) + "\t" + str(average_length_trimmed_reads_decode) + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\t" + str(samtools_coverage)+ "\t" + str(samtools_depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(average_length_trimmed_reads_decode) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\t" + str(samtools_coverage)+ "\t" + str(samtools_depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) +"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "w")
		file2.write("sample\t%q30bases\tn_total_reads\taverage_length_reads\ttotal_contigs\tn_contigs>=500\tcontigs>=10000bp\tcontigs>=25000bp\tcontigs>=50000bp\tlargest_contig\tassembly_length\tN50\tL50\tGC\t%_mapped_reads\tquast_depth\tsamtools_horizontal_cov\tsamtools_depth\n")
		#file2.write("sample\ttotal_contigs\tn_contigs>=500\tcontigs>=10000bp\tcontigs>=25000bp\tcontigs>=50000bp\tlargest_contig\tassembly_length\tN50\tL50\tGC\tn_total_reads\t%_mapped_reads\tquast_depth\n")
		#file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\n")
		file2.write(sample + "\t" + str(q30) + "\t" +  str(n_total_reads) + "\t" + str(average_length_trimmed_reads_decode) + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\t" + str(samtools_coverage)+ "\t" + str(samtools_depth)+"\n")
		
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")



def taxonomy_class_16S(dir, sample, scaffold, threads=4):
	print("Start barnnap")
	print("-----------")
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")
	#barrnap --outseq ./PT33-1_16S.fasta --threads 20 ../assembling/scaffolds.fasta
	cmd="barrnap --outseq " + taxonomy_dir + "/" + sample + "_16S.fasta --threads " + threads + " " + scaffold
	print(cmd)
	os.system(cmd)
	#rdp_classifier classify -o ./PT33_1_rdp_class ./PT33-1_16S.fasta
	cmd2="rdp_classifier classify -o " + taxonomy_dir + "/" + sample + "_rdp_class.tab -q " + taxonomy_dir + "/" + sample + "_16S.fasta"
	print(cmd2)
	os.system(cmd2)



def kmerfinder(dir, scaffold):
	print("Start kmerfinder")
	print("----------------")
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")

	cmd="mkdir " + taxonomy_dir +"/KmerFinder"
	print(cmd)
	os.system(cmd)
	kmerfinder_dir=taxonomy_dir + "/KmerFinder"

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'kmerfinder'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_kmerfinder.py ' + scaffold + " " + kmerfinder_dir
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)



def shigapass(dir, scaffold):
	print("Start ShigaPass")
	print("---------------")
	#/home/admingenomica/software/shigapass/ShigaPass/SCRIPT/ShigaPass.sh -l ./shiga_list_1.txt -o ./externas/SRR2994193/shigapass/ -p /home/admingenomica/software/shigapass/ShigaPass/SCRIPT/ShigaPass_DataBases
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")
	
	cmd="mkdir " + taxonomy_dir +"/ShigaPass"
	print(cmd)
	os.system(cmd)

	cmd="echo " + scaffold + "> " + taxonomy_dir +"/ShigaPass/path_assembly.txt"
	print(cmd)
	os.system(cmd)

	cmd="/home/admingenomica/software/shigapass/ShigaPass/SCRIPT/ShigaPass.sh -l " +  taxonomy_dir +"/ShigaPass/path_assembly.txt -o " + taxonomy_dir +"/ShigaPass/ -p /home/admingenomica/software/shigapass/ShigaPass/SCRIPT/ShigaPass_DataBases"
	print(cmd)
	os.system(cmd)



def shigatyper(dir, sample, paired_reads1, paired_reads2):
	print("Start ShigaTyper")
	print("---------------")
	parent_dir=dir
	taxonomy_dir=os.path.join(parent_dir, "taxonomy")
	
	#conda activate shigatyper
	#shigatyper --R1 ./SRR2994193/SRR2994193_1.fastq.gz --R2 ./SRR2994193/SRR2994193_2.fastq.gz -n SRR2994193 -o ./SRR2994193/shigatyper/
	ShigaTyper_dir=taxonomy_dir+"/ShigaTyper"
	cmd="mkdir " + ShigaTyper_dir
	print(cmd)
	os.system(cmd)

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'shigatyper'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_shigatyper.py ' + paired_reads1 + " " + paired_reads2 + " -n " + sample + " -o " + ShigaTyper_dir
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)



def mlst(mlstdir, scaffold):
	print("Start mlst")
	print("----------")
	
	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'mlst'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_mlst.py ' + scaffold + " " + mlstdir
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)
	


def mlst_DTU(dir, scaffold, specie):
	print("Start mlst_DTU")
	print("--------------")
	parent_dir=dir
	mlst_dir=os.path.join(parent_dir, "mlst")
	if specie=="salmonella":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s senterica -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="escherichia_coli":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s ecoli -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="listeria":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s lmonocytogenes -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="klebsiella":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s kpneumoniae -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="campylobacter":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s cjejuni -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="taylorella":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s taylorella -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="pmultocida":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s pmultocida -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="enterococcus_faecium":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s efaecium -p /home/admingenomica/software/mlst_db/ -x"
	if specie=="enterococcus_faecalis":
		cmd="python3 /home/admingenomica/software/mlst/mlst.py -i " + scaffold + " -o " + mlst_dir + " -t " + mlst_dir + " -s efaecalis -p /home/admingenomica/software/mlst_db/ -x"
	print(cmd)
	os.system(cmd)


def seqsero2(dir, scaffold, sample): #Salmonella serotyping
	print("Start SeqSero2")
	print("-------------")

	parent_dir=dir
	seqsero2_dir="/home/admingenomica/analisis/SeqSero2/" + sample

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'seqsero2'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = "/home/admingenomica/software/miniconda3/envs/seqsero2/bin/SeqSero2_package.py -m k -t 4 -i " + scaffold + " -d " + seqsero2_dir 
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)
	
	cmd="mv /home/admingenomica/analisis/SeqSero2/" + sample + "/* " + parent_dir + "/seqsero2"
	print(cmd)
	os.system(cmd)
	


def sistr(dir,sample,merge, scaffold): #Salmonella serotyping
	print("Start sistr")
	print("----------")
	#sistr --qc -vv --alleles-output ./sistr/1215_22-01_allele-results.json --novel-alleles ./sistr/1215_22-01_novel-alleles.fasta --cgmlst-profiles ./sistr/1215_22-01_cgmlst-profiles.csv -f tab -o ./sistr/1215_22-01_sistr-output.tab ./assembling/1215_22-01.fasta
	parent_dir=dir
	sistr_dir=os.path.join(parent_dir, "sistr")

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'sistr'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_sistr.py ' + scaffold + " " + sistr_dir + " " + sample + " " + merge
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)


def serotyfinder(dir, scaffold): #Ecoli serotyping
	print("Start serotyfinder")
	print("------------------")

	parent_dir=dir
	serotypefinder_dir=os.path.join(parent_dir, "serotypefinder")
	
	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'serotypefinder'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_serotypefinder.py ' + scaffold + " " + serotypefinder_dir
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)


def fimtype(dir, scaffold): #Characterize FimH type of E. coli.
	print("Start fimtype")
	print("-------------")
	parent_dir=dir
	fimtype_dir=os.path.join(parent_dir, "fimtype")

	cmd="perl /home/admingenomica/software/fimtyper/fimtyper.pl -i " + scaffold + " -d /home/admingenomica/software/fimtyper/fimtyper_db -k 95.00 -l 0.60 -o " + fimtype_dir
	print(cmd)
	os.system(cmd)


def patho_typing(dir, paired_reads1, paired_reads2, threads=4): #Pathotyping E. coli
	print("Start patho_typing")
	print("------------------")
	#patho_typing.py -f ./$line/$line*"_R1_"*".fastq.gz" ./$line/$line*"_R2_"*".fastq.gz" -j 10 -o ./$line/pathotyping/ -s Escherichia coli

	parent_dir=dir
	pathotyping_dir=os.path.join(parent_dir, "pathotyping")

	cmd="patho_typing.py -f " + paired_reads1 + " " + paired_reads2 + " -j " + threads + " -o " + pathotyping_dir + " -s Escherichia coli"
	print(cmd)
	os.system(cmd)


def lissero(dir, sample, scaffold): #Listeria serotyping
	print("Start lissero")
	print("-------------")
	parent_dir=dir
	lissero_dir=os.path.join(parent_dir, "lissero")
	cmd="lissero " + scaffold + " > " + lissero_dir + "/" + sample+"_serotyping.tab"

	print(cmd)
	os.system(cmd)


def resfinder(dir, scaffold, specie="nan"): #Identification of AMR genes.
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Resfinder")
	print("---------------")
	parent_dir=dir
	resfinder_dir=os.path.join(parent_dir, "resfinder")
	
	if (specie=="escherichia_coli"):
		specie="ecoli"

	if (specie == "nan"):
		cmd="python /home/admingenomica/software/resfinder/run_resfinder.py -ifa " + scaffold + " -o " + resfinder_dir +" -acq -db_res /home/admingenomica/software/resfinder_db/resfinder_db/"
		#cmd = "python /home/admingenomica/.local/lib/python3.9/site-packages/resfinder/run_resfinder.py -ifa " + scaffold + " -o " + resfinder_dir + " -acq -db_res /home/admingenomica/software/resfinder_db/resfinder_db/ -d -db_disinf /home/admingenomica/software/resfinder_db/disinfinder_db/"
	else:
		cmd = "python /home/admingenomica/software/resfinder/run_resfinder.py -ifa " + scaffold + " -o " + resfinder_dir + " -s " + specie + " -acq -t_p 0.90 -l_p 0.60 -db_res /home/admingenomica/software/resfinder_db/resfinder_db/ -c -db_point /home/admingenomica/software/resfinder_db/pointfinder_db/"
		#cmd = "python /home/admingenomica/.local/lib/python3.9/site-packages/resfinder/run_resfinder.py -ifa " + scaffold + " -o " + resfinder_dir + " -s " + specie + " -acq -db_res /home/admingenomica/software/resfinder_db/resfinder_db/ -d -db_disinf /home/admingenomica/software/resfinder_db/disinfinder_db/ -c -db_point /home/admingenomica/software/resfinder_db/pointfinder_db/"
	print(cmd)
	os.system(cmd)


def rgi_card(dir, sample, scaffold): #ESTA EN AMBIENTE CONDA. EDITAR FUNCION PARA QUE SE EJECUTE CON CONDA.
	print("Start rgi-card")
	print("--------------")
	parent_dir=dir
	card_dir=os.path.join(parent_dir, "card")
    #rgi main --input_sequence /path/to/nucleotide_input.fasta --output_file /path/to/output_file --input_type contig --local --clean
	cmd1="rgi load --card_json /home/vserver1/software/card/card.json"
	os.system(cmd1)

	cmd="rgi main --input_sequence " + scaffold + " --output_file " + card_dir + "/" + sample + "_card --input_type contig --clean" 
	print(cmd)
	os.system(cmd)


def resistencias(dir, sample, resistancedir, specie="nan"):
	print("Start Merging resistance info")
	print("-----------------------------")
	parent_dir=dir
	resfinder_file=parent_dir + "/resfinder/ResFinder_results_tab.txt"
	#resfinder_file=parent_dir + "/amr_genes/resfinder/results_tab.tsv"


	resistencias_dict={}

	listadegenes=[]
	for line in open (resfinder_file, "r"):
		if not line.startswith("Resistance gene"):
			fields=line.split("\t")
			ant=fields[7]
			gene=fields[0]
			listadegenes.append(gene)
			if ant in resistencias_dict: resistencias_dict[ant].append(gene)
			else: resistencias_dict[ant]=[gene]
	file_lista=open(resistancedir+"/resistances_results_list.txt", "a")
	file_lista.write(sample+"\t"+", ".join(listadegenes)+"\n")

	#Anadir datos de resistencias en el fichero de resistencias del directorio merge_res.
	if os.path.isfile(resistancedir + "/resistances_results.xlsx"):	
		df=pd.read_excel(resistancedir + "/resistances_results.xlsx", index_col=0)
		data=df.to_dict()

		gene_list=[]
		for ant in resistencias_dict:
			for gene in resistencias_dict[ant]: gene_list.append(gene)
		for gene in gene_list:
			if gene in data:
				data[gene]["_"+sample]="1"
			else:
				data[gene]={"_"+sample:"1"}



		file=pd.DataFrame(data)
		file.to_excel(resistancedir+ "/resistances_results.xlsx")

	else:
		data={"sample":["_"+sample]}
		for ant in resistencias_dict:
			for gene in resistencias_dict[ant]:
				data[gene]=["1"]

		file=pd.DataFrame(data)
		#print(file)
		file.to_excel(resistancedir+ "/resistances_results.xlsx", index=False)


	##############################################
	#Añadir los datos fenotipicos al fichero de fenotipos del directorio merge_res
	
	#Obtener los datos fenotipicos de la muestra
	pheno_muestra_list=[]
	genes_muestra_list=[]
	resfinder_pheno_file=parent_dir + "/resfinder/pheno_table.txt"
	for line in open (resfinder_pheno_file, "r"):
		if (not line.startswith("#")) and (line != "\n") and (line.split("\t")==5):
			#print(line)
			line=line.rstrip()
			fields=line.split("\t")
			#print(fields[4])
			is_resistant=fields[2]
			if (is_resistant=="Resistant"):
				res_class=fields[0]
				pheno=fields[1]
				gene=fields[4]
				if (len(gene.split(", ")) > 1):
					for i in gene.split(", "):
						list_to_append=[res_class, pheno, i]
						pheno_muestra_list.append(list_to_append)
						genes_muestra_list.append(i)
				else:
					list_to_append=[res_class, pheno, gene]
					pheno_muestra_list.append(list_to_append)
					genes_muestra_list.append(gene)
	

	#Si ya existe el dichero con los fenotipos de todas las muestras en el directorio del merge_res...
	if os.path.isfile(resistancedir + "/pheno_resistances_results.xlsx"):
		df=pd.read_excel(resistancedir + "/pheno_resistances_results.xlsx")
		data=df.to_dict(orient="list")
		
		#Creamos una lista con las muestras y los genes que ya estan en el excel del merge_res
		muestras_ya_en_excel_list= list(data.keys())[3:]
		genes_de_muestras_ya_en_excel= []
		for i in data["Genes"]:
			genes_de_muestras_ya_en_excel.append(i)
		
		data[sample]=[]
		#Iteración por cada gen que ya esta en el excel del merge_res para ver si la muestra que estamos chequeando lo tiene
		for gene_excel in genes_de_muestras_ya_en_excel:
			if gene_excel in genes_muestra_list:
				data[sample].append("1")
			else:
				data[sample].append("0")

		#Añadimos al excel los nuevos genes de la muestra que estamos chequeando y que no estaban incluidos en el excel del merge_res
		#already_checked_genes=[]
		for pheno in pheno_muestra_list:
			res_class=pheno[0]
			phenotype=pheno[1]
			gene=pheno[2]
			
			if (gene not in genes_de_muestras_ya_en_excel):
				for muestras_ya_en_excel in muestras_ya_en_excel_list: #añadimos un 0 en las muestras antiguas que no tienen el nuevo gen que estamos añadiendo
					data[muestras_ya_en_excel].append("0")
				#añadimos los datos del nuevo gen y un 1 en la muestra que estamos chequeando
				data["Resistance_class"].append(res_class)
				data["Phenotype"].append(phenotype)
				data["Genes"].append(gene)
				data[sample].append("1")

		file=pd.DataFrame(data)
		file.to_excel(resistancedir+ "/pheno_resistances_results.xlsx", index=False)
		
	else: #si el fichero con los datos fenotípicos de todas las muestras en el directorio merge_res aún no existe, lo creamos con los genes de la primera muestra
		data={"Resistance_class":[],
		"Phenotype":[],
		"Genes":[],
		}
		
		for pheno in pheno_muestra_list:
			res_class=pheno[0]
			phenotype=pheno[1]
			gene=pheno[2]
			data["Resistance_class"].append(res_class)
			data["Phenotype"].append(phenotype)
			data["Genes"].append(gene)
			if sample not in data:
				data[sample]=["1"]
			else:
				data[sample].append("1")
		
		file=pd.DataFrame(data)
		file.to_excel(resistancedir + "/pheno_resistances_results.xlsx", index=False)
		

def merge_pointfinder_results(dir, sample, pointfinderdir):
	print("Start Merging pointfinder info")
	print("--------------------------------")
	parent_dir=dir
	sample_ed=sample.split("_")[0]
	sample="_"+sample
	pointfinder_file=parent_dir + "/resfinder/PointFinder_results.txt"
	#pointfinder_file=parent_dir + "/pointfinder/assembly_blastn_results.tsv"
	#pointfinder_file=parent_dir + "/pointfinder/" + sample+ "_blastn_results.tsv"
	point_list=[]
	for line in open (pointfinder_file, "r"):
		if not line.startswith("Mutation"):
			fields=line.split("\t")
			point=fields[0]
			#point=fields[0]+"-"+fields[2] 
			point_list.append(point)



	#Anadir datos de mutaciones puntuales en el fichero de mutaciones puntuales de todas las cepas.
	if os.path.isfile(pointfinderdir + "/pointfinder_results.xlsx"):	
		df=pd.read_excel(pointfinderdir + "/pointfinder_results.xlsx", index_col=0)
		data=df.to_dict()

		for point in point_list:
			if point in data:
				data[point][str(sample)]="1"
			else:
				data[point]={str(sample):"1"}

		if len(point_list)==0:
			for point in data:
				data[point][str(sample)]=""

		file=pd.DataFrame(data)
		file.to_excel(pointfinderdir+ "/pointfinder_results.xlsx")

	else:
		data={"sample":[str(sample)]}
		for point in point_list:
			data[point]=["1"]

		file=pd.DataFrame(data)
		print(file)
		file.to_excel(pointfinderdir+ "/pointfinder_results.xlsx", index=False)
		


def resistencias2(dir, sample, resistancedir):
	print("Start Merging resistance info2")
	print("------------------------------")
	parent_dir=dir
	resfinder_file=parent_dir + "/resfinder/results_tab.tsv"
	#resfinder_file=parent_dir + "/amr_genes/resfinder/results_tab.tsv"

	resistencias_dict={}

	for line in open (resfinder_file, "r"):
		if not line.startswith("Database"):
			fields=line.split("\t")
			ant=fields[0]
			gene=fields[1]

			if ant in resistencias_dict: resistencias_dict[ant].append(gene)
			else: resistencias_dict[ant]=[gene]


	#Incluimos los genes de card
	card_file=parent_dir + "/card/" + sample + "_card.txt"
	#card_file=parent_dir + "/amr_genes/rgi_card/" + sample + "_card.txt"
	for line in open(card_file, "r"):
		if not line.startswith("ORF_ID"):
			fields=line.split("\t")
			ant=fields[14]
			gene=fields[8]
			if "OXA" in gene:
				gene="bla"+gene
			isalready="NO"
			for drugclass in resistencias_dict:
				if drugclass in ant:
					isalready="YES"
					if not gene in resistencias_dict[drugclass]:
						resistencias_dict[drugclass].append(gene)
			if isalready == "NO":	
				resistencias_dict[ant]=[gene]



	#Anadir datos de resistencias en el fichero de resistencias de todas las cepas.
	if os.path.isfile(resistancedir + "/resistances_results2.xlsx"):	
		df=pd.read_excel(resistancedir + "/resistances_results2.xlsx", index_col=0)
		data=df.to_dict()

		gene_list=[]
		for ant in resistencias_dict:
			for gene in resistencias_dict[ant]: gene_list.append(gene)
		for gene in gene_list:
			if gene in data:
				data[gene]["_"+sample]="1"
			else:
				data[gene]={"_"+sample:"1"}


		file=pd.DataFrame(data)
		file.to_excel(resistancedir+ "/resistances_results2.xlsx")

	else:
		data={"sample":["_"+sample]}
		for ant in resistencias_dict:
			for gene in resistencias_dict[ant]:
				data[gene]=["1"]

		file=pd.DataFrame(data)
		#print(file)
		file.to_excel(resistancedir + "/resistances_results2.xlsx", index=False)


def plasmidfinder(dir, scaffold):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Plasmidfinder")
	print("-------------------")
	parent_dir=dir
	plasmidfinder_dir=os.path.join(parent_dir, "plasmidfinder")
	#/home/vserver1/software/plasmidfinder/plasmidfinder.py -i ./MS3828/hybrid_assembling_2/assembly.fasta -o ./MS3828/plasmids/ -d /home/vserver1/software/plasmidfinder_db/plasmidfinder_db/
	cmd="/home/admingenomica/software/plasmidfinder/plasmidfinder.py -i " + scaffold + " -o " + plasmidfinder_dir + " -x -p /home/admingenomica/software/plasmidfinder_db/plasmidfinder_db/" + " -t 0.95 -l 0.60"
	print(cmd)
	os.system(cmd)



def merge_plasmidfinder_results(dir, sample, plasmidfinderdir):
	print("Start Merging plasmidfinder info")
	print("--------------------------------")
	parent_dir=dir
	plasmidfinder_file=parent_dir + "/plasmidfinder/results_tab.tsv"
	plasmid_list=[]
	for line in open (plasmidfinder_file, "r"):
		if not line.startswith("Database"):
			fields=line.split("\t")
			plasmid=fields[1]
			plasmid_list.append(plasmid)

	file_lista=open(plasmidfinderdir+"/plasmidfinder_results_list.txt", "a")
	file_lista.write(sample+"\t"+", ".join(plasmid_list)+"\n")

	#Anadir datos de plasmidos en el fichero de plasmidos de todas las cepas.
	if os.path.isfile(plasmidfinderdir + "/plasmidfinder_results.xlsx"):	
		df=pd.read_excel(plasmidfinderdir + "/plasmidfinder_results.xlsx", index_col=0)
		data=df.to_dict()

		for plasmid in plasmid_list:
			if plasmid in data:
				data[plasmid][sample]="1"
			else:
				data[plasmid]={sample:"1"}

		if len(plasmid_list)==0:
			for plasmid in data:
				data[plasmid][sample]=""

		file=pd.DataFrame(data)
		file.to_excel(plasmidfinderdir+ "/plasmidfinder_results.xlsx")

	else:
		data={"sample":[sample]}
		for plasmid in plasmid_list:
			data[plasmid]=["1"]

		file=pd.DataFrame(data)
		print(file)
		file.to_excel(plasmidfinderdir+ "/plasmidfinder_results.xlsx", index=False)



def abricate(dir, sample, scaffold):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Abricate - virulence factors")
	print("----------------------------------")
	parent_dir=dir
	abricate_dir=os.path.join(parent_dir, "virulence_factors")

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'abricate'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_abricate.py ' + scaffold + " " + abricate_dir + " " + sample
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)


def merge_abricate_results(dir, sample, abricatedir):
	print("Start Merging abricate virulence factors info")
	print("--------------------------------")
	parent_dir=dir
	abricate_file=parent_dir + "/virulence_factors/" + sample + "_vFactors.tsv"
	vfactors_list=[]
	for line in open (abricate_file, "r"):
		if not line.startswith("#"):
			fields=line.split("\t")
			vfactor=fields[5]
			vfactors_list.append(vfactor)

	file_lista=open(abricatedir+"/vFactors_results_list.txt", "a")


	if os.path.isfile(abricatedir + "/vFactors_results.xlsx"):	
		df=pd.read_excel(abricatedir + "/vFactors_results.xlsx", index_col=0)
		data=df.to_dict()
		#print("SEGUNDA")
		for vfactor in vfactors_list:
			if vfactor in data:
				data[vfactor][sample]="1"
			else:
				data[vfactor]={sample:"1"}

		if len(vfactors_list)==0:
			for vfactor in data:
				data[vfactor][sample]=""

		file=pd.DataFrame(data)
		file.to_excel(abricatedir+ "/vFactors_results.xlsx")

	else:
		#print("PRIMERA")
		data={"sample":[sample]}
		for vfactor in vfactors_list:
			data[vfactor]=["1"]

		file=pd.DataFrame(data)
		#print(file)
		file.to_excel(abricatedir + "/vFactors_results.xlsx", index=False)


	vfactors_ordered=[]
	for i in data:
		if i in vfactors_list:
			vfactors_ordered.append(i)
	for i in vfactors_list:
		if i not in vfactors_ordered:
			vfactors_ordered.append(i)

	file_lista.write(sample+"\t"+", ".join(vfactors_ordered)+"\n")

def virulencefinder(dir, scaffold, specie="nan"):
	print("Start virulencefinder")
	print("---------------------")
	parent_dir=dir
	virulencefinder_dir=os.path.join(parent_dir, "virulencefinder")
	if specie=="nan" or specie != "listeria" or specie != "escherichia_coli":
		cmd="python /home/admingenomica/software/virulencefinder/virulencefinder.py -i " + scaffold + " -o " + virulencefinder_dir + " -p /home/admingenomica/software/virulencefinder_db/virulencefinder_db/ -l 0.6 -t 0.9 -x"
		print(cmd)
		os.system(cmd)
	else:
		if specie == "listeria":
			cmd="python /home/admingenomica/software/virulencefinder/virulencefinder.py -i " + scaffold + " -o " + virulencefinder_dir + " -p /home/admingenomica/software/virulencefinder_db/virulencefinder_db/ -d listeria -l 0.6 -t 0.9 -x"
			
		if specie == "escherichia_coli":
			cmd="python /home/admingenomica/software/virulencefinder/virulencefinder.py -i " + scaffold + " -o " + virulencefinder_dir + " -p /home/admingenomica/software/virulencefinder_db/virulencefinder_db/ -d virulence_ecoli -l 0.6 -t 0.9 -x"

	print(cmd)
	os.system(cmd)


def merge_virulencefinder_results(dir, sample, virulencefinderdir):
	print("Start Merging virulencefinder virulence factors info")
	print("--------------------------------")
	parent_dir=dir
	virulencefinder_file=parent_dir + "/virulencefinder/results_tab.tsv"
	vfactors_list=[]
	for line in open (virulencefinder_file, "r"):
		if not line.startswith("Database"):
			fields=line.split("\t")
			vfactor=fields[1]
			vfactors_list.append(vfactor)

	file_lista=open(virulencefinderdir+"/virulencefinder_results_list.txt", "a")
	file_lista.write(sample+"\t"+", ".join(vfactors_list)+"\n")


	if os.path.isfile(virulencefinderdir + "/virulencefinder_results.xlsx"):	
		df=pd.read_excel(virulencefinderdir + "/virulencefinder_results.xlsx", index_col=0)
		data=df.to_dict()

		for vfactor in vfactors_list:
			if vfactor in data:
				data[vfactor][sample]="1"
			else:
				data[vfactor]={sample:"1"}

		if len(vfactors_list)==0:
			for vfactor in data:
				data[vfactor][sample]=""

		file=pd.DataFrame(data)
		file.to_excel(virulencefinderdir + "/virulencefinder_results.xlsx")

	else:
		data={"sample":[sample]}
		for vfactor in vfactors_list:
			data[vfactor]=["1"]

		file=pd.DataFrame(data)
		#print(file)
		file.to_excel(virulencefinderdir + "/virulencefinder_results.xlsx", index=False)


def recycler(dir, sample, assembly_graph):
	print("Start recycler")
	print("--------------")
	
	parent_dir=dir
	recycler_dir=os.path.join(parent_dir, "recycler")
	
	cmd="make_fasta_from_fastg.py -g " + assembly_graph + " -o " + recycler_dir + "/" + sample + "_assGraphNodes.fasta"
	print(cmd)
	os.system(cmd)
	
	cmd2="bwa index " + recycler_dir + "/" + sample + "_assGraphNodes.fasta"
	print(cmd2)
	os.system(cmd2)
	
	cmd3="bwa mem -t 20 " + recycler_dir + "/" + sample + "_assGraphNodes.fasta " + parent_dir + "/primario/" + sample + "_trimmed_1P.fastq.gz " + parent_dir + "/primario/" + sample + "_trimmed_2P.fastq.gz | /home/admingenomica/software/samtools-1.21/bin/samtools view -buS - > " + recycler_dir + "/" + sample + "_reads_pe.bam"
	print(cmd3)
	os.system(cmd3)
	cmd4="/home/admingenomica/software/samtools-1.21/bin/samtools view -bF 0x0800 " + recycler_dir + "/" + sample + "_reads_pe.bam > " + recycler_dir + "/" + sample + "_reads_pe_primary.bam"
	print(cmd4)
	os.system(cmd4)
	cmd5="/home/admingenomica/software/samtools-1.21/bin/samtools sort " + recycler_dir + "/" + sample + "_reads_pe_primary.bam -o " + recycler_dir + "/" + sample + "_reads_pe_primary.sort.bam"
	print(cmd5)
	os.system(cmd5)
	cmd6="/home/admingenomica/software/samtools-1.21/bin/samtools index " + recycler_dir + "/" + sample + "_reads_pe_primary.sort.bam"
	print(cmd6)
	os.system(cmd6)
	
	cmd7="rm -rf " + recycler_dir + "/" + sample + "_reads_pe.bam"
	os.system(cmd7)
	cmd8="rm -rf " + recycler_dir + "/" + sample + "_reads_pe_primary.bam"
	os.system(cmd8)
	
	cmd9="/home/admingenomica/software/Recycler/build/scripts-2.7/recycle.py -g " + assembly_graph + " -k 55 -b " + recycler_dir + "/" + sample + "_reads_pe_primary.sort.bam -i True -o " + recycler_dir
	print(cmd9)
	os.system(cmd9)
	
	plasmidfinder_recycler_dir=os.path.join(recycler_dir, "plasmidfinder")
	
	if not os.path.exists(plasmidfinder_recycler_dir):
		os.mkdir(plasmidfinder_recycler_dir)

	cmd="/home/admingenomica/software/plasmidfinder/plasmidfinder.py -i " + recycler_dir + "/assembly_graph.cycs.fasta -o " + plasmidfinder_recycler_dir + " -x -p /home/admingenomica/software/plasmidfinder_db/plasmidfinder_db/" + " -t 0.80"
	print(cmd)
	os.system(cmd)
	
	resfinder_recycler_dir=os.path.join(recycler_dir, "resfinder")
	cmd="python /home/admingenomica/software/resfinder/run_resfinder.py -ifa " + recycler_dir + "/assembly_graph.cycs.fasta -o " + resfinder_recycler_dir + " -acq -db_res /home/admingenomica/software/resfinder/db_resfinder_10_06_22" 
	print(cmd)
	os.system(cmd)



def mob_suite (dir, scaffold):
	print("Start mob_suite")
	print("------------")
	
	parent_dir=dir
	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'mob_suite'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_mob_suite.py ' + scaffold + " " + parent_dir
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)

	

def abacas(dir, sample, scaffold, reference):
	print("Start abacas")
	print("------------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	cmd="abacas -r " + reference + " -q " + scaffold + " -p nucmer -o " + assembling_dir + "/" + sample + " -m"
	print(cmd)
	os.system(cmd)



def integronfinder(dir, scaffold, threads):
	print("Start integronfinder")
	print("-----------------")
	parent_dir=dir
	integronfinder_dir=os.path.join(parent_dir, "integronfinder")
	#integron_finder --local-max --func-annot --cpu 20 --pdf ./assembly_graph.cycs.fasta
	cmd="integron_finder --local-max --func-annot --cpu " + threads + " --pdf  --outdir " + integronfinder_dir + " " + scaffold
	print(cmd)
	os.system(cmd)



def prokka(dir, sample, scaffold, threads=4):
	#do /home/vserver1/software/prokka/bin/prokka --outdir ./$c/annotation --force --prefix $c"_ann" ./$c/annotation/$c"_assembling.result.fasta"
	print("Start Prokka")
	print("------------")
	parent_dir=dir
	annotation_dir="/home/admingenomica/analisis/temp/"+sample
	
	cmd= "/home/admingenomica/software/prokka/bin/prokka --kingdom Bacteria --proteins /home/admingenomica/software/resfinder/db_resfinder/resfinder_all_prot_db.fasta --outdir " + annotation_dir + " --gcode 11 --cpus " + threads + " --prefix " + sample + " " + scaffold
	#cmd= "/home/admingenomica/software/prokka/bin/prokka --kingdom Bacteria --proteins /2023/EHDV_NGS_2023_001/referencia_culicoides/annotation/ena_PRJEB19938_protein.fasta --outdir " + annotation_dir + " --gcode 11 --cpus " + threads + " --prefix " + sample + " " + scaffold
	#cmd= "/home/admingenomica/software/prokka/bin/prokka --kingdom Viruses --proteins /home/admingenomica/software/resfinder/db_resfinder/resfinder_all_prot_db.fasta --outdir " + annotation_dir + " --gcode 11 --cpus " + threads + " --prefix " + sample + " " + scaffold
	
	print(cmd)
	os.system(cmd)	
	cmd2="mv /home/admingenomica/analisis/temp/"+sample+"/* "+dir+"/annotation/"
	#cmd2="mv /home/admingenomica/analisis/temp/"+sample+" "+dir+"/"+sample+"_ann"
	print(cmd2)
	os.system(cmd2)	


def bwa(dir, sample, reads1, reads2, reference, threads=4):
	print("Start bwa")
	print("---------")
	#bwa mem -t 14 ./data/Reference.fna ./data/illumina/Illumina_R1.fastq.gz ./data/illumina/Illumina_R2.fastq.gz > ./mapping/illumina.bwa.sam
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	#secondary_dir=os.path.join(parent_dir, "secundario_ref_pOX_21_1088-10")
	s=" "
	
	cmd1="/home/admingenomica/software/bwa/bwa mem -R \"@RG\\tID:"+ sample +"\\tSM:"+ sample +"\\tPL:ILLUMINA\" -t " + threads + s + reference + s + reads1 + s + reads2 
	cmd2= " | /home/admingenomica/software/samtools-1.21/bin/samtools view -Sbh -> "+ secondary_dir + "/" + sample + "_bwa.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools sort -@" + threads + s + secondary_dir + "/" + sample + "_bwa.bam -o " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools index " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	cmd="java -jar /home/admingenomica/software/picard-2-27-1/picard.jar MarkDuplicates -I " + secondary_dir + "/" + sample + "_bwa_sort.bam -O " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam -M picard_metrics.txt"
	os.system(cmd)
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools index " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam"

	
#Alineamiento con bowtie (bwa da mejores resultados)
def bowtie(dir, sample, reads1, reads2, reference, threads=4):
	print("Start bowtie")
	print("------------")
	#/home/vserver1/software/bowtie2-2.4.2-linux-x86_64/bowtie2
	#crear indice de referencia: bowtie2-build ./referencia_mpox/ON563414_monkey_pox_virus.fasta ON563414_monkey_pox_virus
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	s=" "
	cmd1="bowtie2 -p " + threads  + s + "--met-file" + s + secondary_dir + "/" + sample + "_bwt_metrics.txt" +  s + "--rg-id " + sample + " --rg SM:" + sample + " --rg PL:ILLUMINA  -x" + s + reference[:-4] + s + "-1" + s + reads1 + s + "-2" + s + reads2
	cmd2= " | /home/admingenomica/software/samtools-1.21/bin/samtools view -Sbh - > "+ secondary_dir + "/" + sample + "_bwa.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	#samtools sort MS3825.bam -o MS3825.sorted.bam
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools sort -@" + threads + s + secondary_dir + "/" + sample + "_bwa.bam -o " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	#samtools index MS3825.sorted.bam
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools index " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	cmd="java -jar /home/admingenomica/software/picard-2-27-1/picard.jar MarkDuplicates -I " + secondary_dir + "/" + sample + "_bwa_sort.bam -O " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam -M picard_metrics.txt"
	os.system(cmd)
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools index " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam"




def qualimap(dir, sample):
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	#secondary_dir=os.path.join(parent_dir, "secundario_ref_pOX_21_1088-10")
	cmd="qualimap bamqc -sd -bam " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam"
	#cmd="qualimap bamqc -sd -bam " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	print(cmd)
	os.system(cmd)


def picard_metrics(dir, sample, reference):
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	#java -jar /home/admingenomica/software/picard-2-27-1/picard.jar CollectWgsMetrics -I ./secundario/2022-003592-1_bwa_sort_markedDup.bam -O ./secundario/picard_metrics.tab -R ../SRR16920742/assembling/SRR16920742_less_500.fasta
	cmd="java -jar /home/admingenomica/software/picard-2-27-1/picard.jar CollectWgsMetrics -I " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam -O " + secondary_dir + "/" + sample + "_picard_metrics.tab -R " + reference
	print(cmd)
	os.system(cmd)


def samtools_coverage(dir, sample, coveragedir):
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools coverage " + secondary_dir + "/" + sample + "_bwa_sort_markedDup.bam > " + secondary_dir + "/" + sample + "_samtools_coverage.tab"
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

	if os.path.isfile(coveragedir+ "/coverage_by_segments.txt"):	
		file2=open(coveragedir + "/coverage_by_segments.txt", "a")
		file2.write(sample + "\t" + "\t".join(list_depth_by_seg)+ "\t" + "\t".join(list_cov_by_seg)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:	
		file2=open(coveragedir + "/coverage_by_segments.txt", "a")
		file2.write("sample\tdepth_1\tdepth_2\tdepth_3\tdepth_4\tdepth_5\tdepth_6\tdepth_7\tdepth_8\tdepth_9\tdepth_10\tcov_1\tcov_2\tcov_3\tcov_4\tcov_5\tcov_6\tcov_7\tcov_8\tcov_9\tcov_10\n")
		file2.write(sample + "\t" + "\t".join(list_depth_by_seg)+ "\t" + "\t".join(list_cov_by_seg)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")


def coverage_alignment(dir, sample, coveragedir):
	print("Start merging coverage alignment results")
	print("-----------------------------------------")
	parent_dir=dir
	secundario_dir=os.path.join(parent_dir, "secundario")
	ref_length=0
	ass_length=0
	n50=0
	hay_dup="no"
	for line in open(secundario_dir+"/"+sample+"_bwa_sort_markedDup_stats/genome_results.txt", "r"):
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
			depth=fields_list[-1]
		elif "of reference with a coverageData >= 20X" in line:
			fields_list=line.split(" ")
			coverage=fields_list[-8]

		elif "NC_004297.1" in line:
			fields_list=line.split("\t")
			depth_1=fields_list[4][0:6]
		elif "NC_004296.1" in line:
			fields_list=line.split("\t")
			depth_2=fields_list[4][0:6]
		

	if hay_dup=="no":
		dupreads="0"	
	if os.path.isfile(coveragedir+ "/alignment_coverage_stats.txt"):	
		file2=open(coveragedir + "/alignment_coverage_stats.txt", "a")
		file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(dupreads) + "\t" + str(coverage) + "\t" + str(depth)+"\n")#"\t" +str(depth_1) +"\t"+str(depth_2)+"\t"+str(depth_3)+"\t"+str(depth_4)+"\t"+str(depth_5)+"\t"+str(depth_6)+"\t"+str(depth_7)+"\t"+str(depth_8)+"\t"+str(depth_9)+"\t"+str(depth_10)+"\n")
		#file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(dupreads) + "\t" + str(coverage) + "\t" + str(depth)+"\t" +str(depth_1) +"\t"+str(depth_2)+"\t"+str(depth_3)+"\t"+str(depth_4)+"\t"+str(depth_5)+"\t"+str(depth_6)+"\t"+str(depth_7)+"\t"+str(depth_8)+"\t"+str(depth_9)+"\t"+str(depth_10)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:
		file2=open(coveragedir + "/alignment_coverage_stats.txt", "w")
		file2.write("sample\ttotal_reads\tmapped_reads\tpercent_mapped_reads\tdup_reads\thorizontal_coverage\tdepth\n")#\tdepth_1\tdepth_2\tdepth_3\tdepth_4\tdepth_5\tdepth_6\tdepth_7\tdepth_8\tdepth_9\tdepth_10\n")
		#file2.write("sample\ttotal_reads\tmapped_reads\tpercent_mapped_reads\tdup_reads\thorizontal_coverage\tdepth\tdepth_1\tdepth_2\tdepth_3\tdepth_4\tdepth_5\tdepth_6\tdepth_7\tdepth_8\tdepth_9\tdepth_10\n")
		file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(dupreads) + "\t" + str(coverage) + "\t" + str(depth)+"\n")#\t"+str(depth_1) +"\t"+str(depth_2)+"\t"+str(depth_3)+"\t"+str(depth_4)+"\t"+str(depth_5)+"\t"+str(depth_6)+"\t"+str(depth_7)+"\t"+str(depth_8)+"\t"+str(depth_9)+"\t"+str(depth_10)+"\n")
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
	cmd1="/home/admingenomica/software/bcftools_1.21/bcftools mpileup --min-MQ 40 --min-BQ 30 -I --threads " + threads + " -s " + sample + s + " --annotate INFO/AD -f " + reference + s + bam + s + "| /home/admingenomica/software/bcftools_1.21/bcftools call -mv --ploidy 1 --threads " + threads + " -Ov -o " + terciary_dir + "/" + sample + "_raw.vcf"

	print(cmd1)
	os.system(cmd1)
	
	#cmd2="bcftools norm -Ov -f " + reference + " -d all -o " + terciary_dir + "/" + sample + "_norm_raw.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#print(cmd2)
	#os.system(cmd2)
	
	cmd3="/home/admingenomica/software/bcftools_1.21/bcftools filter -Ov -e 'QUAL<100||PL[0:0]<150||DP<5' --threads " + threads + " -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#cmd3="bcftools filter -Ov -e 'QUAL<100||PL[0:0]<130||DP<5' --threads " + threads + " -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#cmd3="bcftools filter -Ov -e 'QUAL<30 || DP<10' -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_norm_raw.vcf"
	print(cmd3)
	os.system(cmd3)
	
	#cmd4="bcftools view -O v -V indels -o " + terciary_dir + "/" + sample + "_filter.vcf " + parent_dir + "/" + sample + "_AE017334_q30.vcf"
	#print(cmd4)
	#os.system(cmd4)


#OPCIONAL: Si quieres tener un vcf con los snps de las regiones fágicas filtradas, ejecutar esta función (aunque para hacer el análisis filogenético, las regiones fágicas se enmascaran en la función generate_consensus_sequence pasándole el fichero bed con las coordenadas fágicas)
def bedtools_filterPhageRegion(dir, sample, vcf, phages_bed):
	print("Start bedtools - filter snps in phage regions")
	print("---------------------------------------------")
	#bedtools intersect -header -v -a ./variant_calling/brucella_samples_filtered.vcf -b ./phage_regions.bed -wa > ./variant_calling/brucella_samples_filtered_Nophages.vcf
	parent_dir=dir
	terciary_dir=os.path.join(parent_dir, "terciario")
	cmd="bedtools intersect -header -v -a " + vcf + " -b " + phages_bed + " -wa > " + terciary_dir + "/" + sample + "_filter_Nophages.vcf"
	print(cmd)
	os.system(cmd)


def generate_consensus_sequence(dir, sample, vcf, reference, phages_bed="nan"):
	print("Start bcftools - generate consensus sequence")
	print("--------------------------------------------")
	#bcftools consensus -i -s NA001 -f in.fa in.vcf.gz > out.fa
	parent_dir=dir
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
		cmd2="/home/admingenomica/software/bcftools_1.21/bcftools consensus -p " + sample + " -m " + phages_bed + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"
		print(cmd2)
		os.system(cmd2)
	else:
		cmd2="/home/admingenomica/software/bcftools_1.21/bcftools consensus -p " + sample + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"
		print(cmd2)
		os.system(cmd2)	
	#cmd2="bcftools consensus -p " + sample + " -m " + phages_bed + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter_Nophages.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"


##################################MAIN##############################################


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
	required.add_argument("--phages_bed", "-pb", help="Bed file with phages coordinates")
	required.add_argument("--specie", "-sp", help="salmonella, escherichia_coli, campylobacter, staphylococcus_aureus,enterococcus_faecium, enterococcus_faecalis, klebsiella, listeria or pmultocida")
	

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

	if args.reference:
		print("reference = %s" % args.reference)
		reference=args.reference

	if args.mergeresdir:
		print("mergeresdir = %s" % args.mergeresdir)
		mergeresdir=args.mergeresdir

	if args.phages_bed:
		print("phages_bed = %s" % args.phages_bed)
		phages_bed=args.phages_bed

	if args.specie:
		print("specie = %s" % args.specie)
		specie=args.specie
		specie=specie.rstrip()
	
	
	#create subdirectories
	createSubdirectories(output)
	
	#run fastqc
	fastqc(output, sample, threads, reads1, reads2)
	
	#run trimmomatic
	trimmed_paired_reads1, trimmed_paired_reads2, trimmed_unpaired_reads1, trimmed_unpaired_reads2=trimmomatic(output, sample, reads1, reads2, threads)
	
	#run fastqc
	#fastqc(output, sample, threads, trimmed_paired_reads1, trimmed_paired_reads2)
	fastqc(output, sample, threads, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	#fastqc(output, sample, threads, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz")
	
	#run calculate q30 results
	#q30(output, sample, reads1, reads2, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz") #con bash - awk
	calculate_q30_percentage_from_compressed_fastq(output, sample, reads1, reads2, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	
	#run calculate number PF reads post trimmed
	reads_PF_post_trimmeed(output, sample, reads1, reads2, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	
	#run kraken
	kraken(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", threads)
	#kraken(output, sample, output+"/Analisis_contaminacion/all_sin_innocua/"+sample+"_all_sin_innocua_1P.fastq.gz", output+"/Analisis_contaminacion/all_sin_innocua/"+sample+"_all_sin_innocua_2P.fastq.gz", threads)
	
	#run confindr
	confindr(output, "/home/admingenomica/software/confindr_db/")
	
	#run mash_screen
	mash_screen (output, "/home/admingenomica/software/mash_screen_db/RefSeq88n.msh", output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	#mash_screen (output, "/home/admingenomica/software/mash_screen_db/RefSeq88n.msh", output+"/primario/"+sample+"_trimmed_1P_id_modificados.fastq.gz", output+"/primario/"+sample+"_trimmed_2P_id_modificados.fastq.gz")
	#mash_screen (output, "/home/admingenomica/software/mash_screen_db/refseq.genomes.k21s1000.msh", output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	
	
	#run spades
	#spades(output, sample, trimmed_paired_reads1, trimmed_paired_reads2, trimmed_unpaired_reads1, trimmed_unpaired_reads2, threads)
	#spades(output, sample, reads1, reads2, threads)
	spades(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", output+"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads)
	#spades(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", threads)
	#spades(output, sample, output+"/"+sample+"_chaparevirus_with_bwa_ref_6_R1.fastq.gz", output+"/"+sample+"_chaparevirus_with_bwa_ref_6_R2.fastq.gz", threads)
	#spades(output, sample, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz", output + "/trimming/" + sample + "_trimmed_L1_SE.fastq.gz", threads)
	
	#run pilon
	pilon(output, sample, output+"/assembling/scaffolds.fasta", output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", output +"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads)
	#pilon(output, sample, output+"/assembling/scaffolds.fasta", paired_reads1, paired_reads2, unpaired_reads1, unpaired_reads2, threads):
	
	
	#run quast
	#con referencia
	if args.reference:
		quast(output, output+"/assembling/" + sample + ".fasta",  output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", output +"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads, reference)
		#quast(output, sample, output+"/assembling/scaffolds.fasta",  output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", reference, threads)
	else:
		#sin referencia
		quast(output, output+"/assembling/" + sample + ".fasta", output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", output+"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads=threads)
		#quast(output, sample, output+"/assembling/scaffolds.fasta", output +"/"+sample+"_1P.fastq.gz",  output +"/"+sample+"_2P.fastq.gz", threads=threads)
		#quast(output, sample, output+"/assembling/scaffolds.fasta", output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", threads=threads)
		#quast(output, sample, output+"/assembling/scaffolds.fasta", output+"/"+sample+"_chaparevirus_with_bwa_ref_6_R1.fastq.gz", output+"/"+sample+"_chaparevirus_with_bwa_ref_6_R2.fastq.gz", threads=threads)
	
	
	#run coverage_assembly
	coverage_assembly(output, sample, mergeresdir, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	#coverage_assembly(output, sample, mergeresdir, output+"/"+ sample+"_chaparevirus_with_bwa_ref_6_R1.fastq.gz", output+"/"+sample+"_chaparevirus_with_bwa_ref_6_R2.fastq.gz")	

	
	#remove short contigs (<500 pb)
	remove_short_contigs(output, sample, output+"/assembling/" + sample + ".fasta")

	
	#Check contamination and completedness of assembly
	checkm(output, sample, output+"/assembling/" + sample + "_gt500.fasta", threads)
	
	
	#run taxonomy_class_16S
	taxonomy_class_16S(output, sample, output+"/assembling/" + sample + "_gt500.fasta", threads)
	
	
	#run kmerfinder for specie classification
	kmerfinder(output, output+"/assembling/" + sample + "_gt500.fasta")
	
	
	#run ShigaPass and ShigaTyper to distinguish ecoli from shigella
	if (args.specie) and (args.specie == "escherichia_coli"):
		shigapass(output, output+"/assembling/" + sample + "_gt500.fasta")
		shigatyper(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz")
	
	
	#run mlst
	mlst(mergeresdir, output+"/assembling/" + sample + "_gt500.fasta")
	#mlst(mergeresdir, output+"/hybrid_assembling/assembly.fasta")
	#mlst(mlstdir, output+"/assembling/"+sample+"_assembling.result.fasta")
	
	if (args.specie) and (args.specie != "nan"):
		#run mlst_DTU
		mlst_DTU(output, output+"/assembling/" + sample + "_gt500.fasta", specie)	
	
	
	#Salmonella ssp. - Salmonella serotyping
	if (args.specie) and (specie == "salmonella"):
		seqsero2(output, output+"/assembling/" + sample + "_gt500.fasta", sample)
		#sistr(output, sample, mergeresdir, output+"/assembling/" + sample + "_gt500.fasta")
	
	
	#Escherichia coli
	if (args.specie) and (specie == "escherichia_coli"):
		serotyfinder(output, output+"/assembling/" + sample + "_gt500.fasta") #run serotypefinder - E. coli serotyping
		fimtype(output, output+"/assembling/" + sample + "_gt500.fasta") #run firmtyper - Characterize FimH type in E. coli
		patho_typing(output, output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", threads) #run pathotyper - characterize VTEC pathotypes.
	
	
	#run lissero - Listeria serotyping
	if (args.specie) and (specie == "listeria"):
		lissero(output, output+"/assembling/" + sample + "_gt500.fasta")
	
	
	#run resfinder
	species_resfinder=["escherichia_coli", "salmonella", "campylobacter", "klebsiella", "enterococcus_faecalis", "enterococcus_faecium", "staphylococcus_aureus"]
	
	if (args.specie) and (specie in species_resfinder):
		resfinder(output, output+"/assembling/" + sample + "_gt500.fasta", specie)
		#resfinder(output, output+"/assembling/assembly.fasta", specie)
	else:
		resfinder(output, output+"/assembling/" + sample + "_gt500.fasta")
	#resfinder(output, output+"/hybrid_assembling/assembly.fasta")
	#resfinder(output, output+"/"+sample+".fasta")
	#resfinder(output, output+"/assembly_graph.cycs.fasta")
	

	#run rgi-card
	#rgi_card(output, sample, output+"/assembling/scaffolds.fasta")
	#rgi_card(output, sample, output+"/assembling/" + sample + ".contigs_modifiedId.fasta")
	
	
	#run merge resistences info
	if (args.specie) and (specie in species_resfinder):
		resistencias(output, sample, mergeresdir, specie)
		merge_pointfinder_results(output, sample, mergeresdir)
	else: 
		resistencias(output, sample, mergeresdir)
	#Usar resistencias2 cuando activamos rgi-card
	#resistencias2(output, sample, mergeresdir)
	
	
	#run plasmidfinder - Identifying plasmid replicons
	plasmidfinder(output, output+"/assembling/" + sample + "_gt500.fasta")
	#plasmidfinder(output, output+"/assembling/assembly.fasta")
	#plasmidfinder(output, output+"/scaffolds.fasta")
	#plasmidfinder(output, output+"/assembling/" + sample + "_assembling.result.fasta")
	#plasmidfinder(output, output+"/assembly_graph.cycs.fasta")
	
	
	#run merge plasmidfinder results
	merge_plasmidfinder_results(output, sample, mergeresdir)
	
	#run abricate
	abricate(output, sample, output+"/assembling/" + sample + "_gt500.fasta")
	#abricate(output, sample, output+"/hybrid_assembling/assembly.fasta")
	#abricate(output, sample, output+"/assembling/"+sample+"_assembling.result.fasta")
	
	
	#run_merge_abricate_results
	merge_abricate_results(output, sample, mergeresdir)
	
	
	#run virulencefinder
	if (args.specie) and (args.specie != "nan"):
		virulencefinder(output, output+"/assembling/" + sample + "_gt500.fasta", specie)
	else:
		virulencefinder(output, output+"/assembling/" + sample + "_gt500.fasta")
	
	merge_virulencefinder_results(output, sample, mergeresdir)
	

	#Characterization of complete plasmod sequences.
	#run recycler
	recycler(output, sample, output+"/assembling/assembly_graph.fastg")
	
	
	#run mob_suite
	#mob_suite(output, output+"/assembling/" + sample + "_gt500.fasta")
	
	
	#run abacas
	#abacas(output, sample,  output+"/assembling/" + sample + "_gt500.fasta", reference)
	#abacas(output, sample, output+"/assembling/" + sample + ".contigs.fasta", reference)
	
	
	#run integronfinder
	#integronfinder(output, output+"/assembling/" + sample + "_gt500.fasta", threads)
	
	
	
	#run prokka
	#prokka(output, sample, output+"/assembling/"+sample+"_gt500.fasta", threads)
	#prokka(output, sample, output+"/assembling/scaffolds.fasta", threads)
	prokka(output, sample, output+"/assembling/" + sample + "_gt500.fasta", threads)
	#prokka(output, sample, output+"/hybrid_assembling/assembly.fasta", threads)
	#prokka(output, sample, output+"/assembling/"+sample+".MULTIFASTA.fa")
		
	
	
	if args.reference:
		
		#MAPPING
		#run bwa
		#bwa(output, sample, trimmed_paired_reads1, trimmed_paired_reads2, reference, threads)
		bwa(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", reference, threads)
		#bwa(output, sample, output+"/primario/"+sample+"_trimmed_1P_id_modificados.fastq.gz", output+"/primario/"+sample+"_trimmed_2P_id_modificados.fastq.gz", reference, threads)
		#bwa(output, sample, reads1, reads2, reference, threads)
		#bwa(output, sample, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz", reference, threads)
		
		#run bowtie
		#bowtie(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", reference, threads)
		
		#run qualimap
		qualimap(output, sample)
		
		#run picard
		picard_metrics(output, sample, reference)
		
		#samtools coverage (for references with diferent segments or chromosomes)
		#samtools_coverage(output, sample, mergeresdir)
		
		
		#run variantcalling
		#bcftools_varcalling(output, sample, output + "/secundario/" + sample + "_nodup_sort_mapped_rg.bam", reference, threads)
		bcftools_varcalling(output, sample, output + "/secundario/" + sample + "_bwa_sort_markedDup.bam", reference, threads) #bwa
		
		#run bedtools to remove snps in phage coordinates from vcf files
		if args.phages_bed:
			bedtools_filterPhageRegion(output, sample, output + "/terciario/" + sample + "_filter.vcf", phages_bed)
		
		#run bcftools consensus to generate consensus sequence
		if args.phages_bed:
			generate_consensus_sequence(output, sample, output + "/terciario/" + sample + "_filter.vcf", reference, phages_bed)
		else:
			generate_consensus_sequence(output, sample, output + "/terciario/" + sample + "_filter.vcf", reference)
			#generate_consensus_sequence(output, sample, output + "/terciario_ref_pOX_21_1088-10/" + sample + "_filter.vcf", reference)
		
		#run coverage_alignment
		coverage_alignment(output, sample, mergeresdir)		
	
	
if __name__ == "__main__":

    main()
