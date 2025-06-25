import sys
import subprocess
import argparse
import os
import pandas as pd

## Nanopore desde fastq

def createSubdirectories(output_nanopore):
	parent_dir=output_nanopore
	
	primary_dir=os.path.join(parent_dir,"primario")
	assembling_dir=os.path.join(parent_dir, "assembling")
	plasmidfinder_dir=os.path.join(parent_dir, "plasmidfinder")
	resfinder_dir=os.path.join(parent_dir, "resfinder")
	sistr_dir=os.path.join(parent_dir, "sistr")
	
	
	if not os.path.exists(primary_dir):
		os.mkdir(primary_dir)
	if not os.path.exists(assembling_dir):
		os.mkdir(assembling_dir)
	if not os.path.exists(plasmidfinder_dir):
		os.mkdir(plasmidfinder_dir)
	if not os.path.exists(resfinder_dir):
		os.mkdir(resfinder_dir)
	if not os.path.exists(sistr_dir):
		os.mkdir(sistr_dir)

	

def longqc(output_nanopore, minion_read, sample, threads):
	parent_dir=output_nanopore
	out_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="mkdir" + s + out_dir + "/longqc"
	print(cmd1)
	os.system(cmd1)
	cmd2= "python /home/admingenomica/software/LongQC/longQC.py sampleqc -o " + out_dir + "/longqc" + " -x ont-ligation -s " + sample + " -p " + threads + s + parent_dir + "/" + sample + ".fastq.gz"
	print(cmd2)
	os.system(cmd2)

def fastqc(output_nanopore, minion_read, sample, threads):
	parent_dir=output_nanopore
	fastqc_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="mkdir" + s + fastqc_dir + "/fastqc"
	print(cmd1)
	os.system(cmd1)
	cmd2="fastqc -t " + threads + " -o " + fastqc_dir + "/fastqc/" + s + parent_dir + "/" + sample + ".fastq.gz"
	print(cmd2)
	os.system(cmd2)

def porechop(output_nanopore, sample, threads):
	parent_dir=output_nanopore
	porechop_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="porechop -i" + s + parent_dir + "/" + sample + ".fastq.gz" + s +  "-t " + threads + s + "-v 1 -o" + s + porechop_dir + "/fastqc/" + sample + "_trimmed.fastq.gz"
	print(cmd1)
	os.system(cmd1)

def nanofilt(output_nanopore, trimming, sample):
	parent_dir=output_nanopore
	nanofilt_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="gunzip -c " + nanofilt_dir + "/fastqc/" + sample + "_trimmed.fastq.gz "
	cmd2="| NanoFilt -q 10 -l 500 "
	cmd3="| gzip > " + nanofilt_dir + "/fastqc/" + sample + "_trimmed_filtered.fastq.gz"
	cmd=cmd1+cmd2+cmd3
	
	print(cmd)
	os.system(cmd)

def fasqc_trimming(output_nanopore, filtered_trimming, sample, threads):
	parent_dir=output_nanopore
	fastqc_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="fastqc -t " + threads + " -o " + fastqc_dir + "/fastqc/" + s + fastqc_dir + "/fastqc/" + sample + "_trimmed_filtered.fastq.gz"

	print(cmd1)
	os.system(cmd1)


def longqc_timming(output_nanopore, filtered_trimming, sample, threads):
	parent_dir=output_nanopore
	longqc_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd="python /home/admingenomica/software/LongQC/longQC.py sampleqc -o " + longqc_dir + s + "-x ont-ligation -s " + sample + s + "-p " + threads + s + longqc_dir + "/" + sample + "_trimmed_filtered.fastq.gz"

	print(cmd)
	os.system(cmd)


def assembl_hibrid(output_illumina, output_nanopore, unpair1, unpair2, sample, pair1, pair2, filtered_trimming, threads):
	
	print("----------------------")
	print("Start hibrid assembly ")
	print("----------------------")
	
	parent_dir_Illumina=output_illumina
	primario_dir_illumina=os.path.join(parent_dir_Illumina, "primario")
	parent_dir_Nanopore=output_nanopore
	assembling_dir=os.path.join(parent_dir_Nanopore, "assembling")
	s=" "
	#Concatener unpairs reads illumina
	cmd="cat " + unpair1 + s + unpair2 + ">" + s + primario_dir_illumina + "/" + sample + "_trimmed_allU.fastq.gz"
	print(cmd)
	os.system(cmd)
	cmd1="gzip -d" + s + filtered_trimming
	print(cmd1)
	os.system(cmd1)
	#con unpairs sin concatenados
	#cmd2="unicycler -1 " + pair1 + s + "-2 " + pair2 + s + "-s " + unpair1 + " -s " + unpair2 + s + "-l " + filtered_trimming + s + "-o" + s + assembl_hibrid + s + "-t " + threads
	#con unpairs concatenados (mejor opcion)
	cmd2="unicycler -1 " + pair1 + s + "-2 " + pair2 + s + "-s " + output_illumina+"/primario/"+sample+"_trimmed_allU.fastq.gz" + s + "-l " + output_nanopore+"/primario/fastqc/"+ sample+"_trimmed_filtered.fastq" + s + "-o" + s + assembling_dir + s + "-t " + threads
	#sin unpairs
	#cmd2="unicycler -1 " + pair1 + s + "-2 " + pair2 + s + "-l " + filtered_trimming + s + "-o" + s + assembl_hibrid + s + "-t " + threads
	print(cmd2)
	os.system(cmd2)
	
	cmd3="gzip -q" + s + output_nanopore+"/primario/fastqc/"+ sample+"_trimmed_filtered.fastq"
	print(cmd3)
	os.system(cmd3)
	
def quast(output_illumina, output_nanopore, sample, pair1, pair2, unpair_concatened, filtered_trimming, threads ):
	print("----------------------")
	print("Start Quast ")
	print("----------------------")
	
	parent_dir_Illumina=output_illumina
	primario_dir_illumina=os.path.join(parent_dir_Illumina, "primario")
	parent_dir_Nanopore=output_nanopore
	assembling_dir=os.path.join(parent_dir_Nanopore, "assembling")
	s=" "

	cmd="quast.py -o " + assembling_dir + s + "-t "+ threads + " --pe1 " + pair1 + s + "--pe2 " + pair2 + s + "--single " + output_illumina+"/primario/"+sample+"_trimmed_allU.fastq.gz --nanopore " + parent_dir_Nanopore + "/primario/fastqc/" + sample + "_trimmed_filtered.fastq.gz" + s + assembling_dir + "/assembly.fasta"
	print(cmd)
	os.system(cmd)



def coverage_assembly(output_nanopore, sample, coveragedir):
	print("Start merging coverage assemblies results")
	print("-----------------------------------------")
	parent_dir=output_nanopore
	assembling_dir=os.path.join(parent_dir, "assembling")
	
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
		#estaaafile2.write("sample\ttotal_contigs\tn_contigs>=500\tcontigs>=10000bp\tcontigs>=25000bp\tcontigs>=50000bp\tlargest_contig\tassembly_length\tN50\tL50\tGC\tn_total_reads\t%_mapped_reads\tquast_depth\tsamtools_horizontal_cov\tsamtools_depth\n")
		file2.write("sample\ttotal_contigs\tn_contigs>=500\tcontigs>=10000bp\tcontigs>=25000bp\tcontigs>=50000bp\tlargest_contig\tassembly_length\tN50\tL50\tGC\tn_total_reads\t%_mapped_reads\tquast_depth\n")
		file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\n")
		#estaafile2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\t" + str(samtools_coverage)+ "\t" + str(samtools_depth)+"\n")
		
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")


def plasmidfinder(output_nanopore, assembly):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Plasmidfinder")
	print("-------------------")
	parent_dir=output_nanopore
	plasmidfinder_dir=os.path.join(parent_dir, "plasmidfinder")
	#/home/vserver1/software/plasmidfinder/plasmidfinder.py -i ./MS3828/hybrid_assembling_2/assembly.fasta -o ./MS3828/plasmids/ -d /home/vserver1/software/plasmidfinder_db/
	cmd="/home/admingenomica/software/plasmidfinder/plasmidfinder.py -i " + assembly + " -o " + plasmidfinder_dir + " -x -p /home/admingenomica/software/plasmidfinder_db/" + " -t 0.80"
	print(cmd)
	os.system(cmd)


def merge_plasmidfinder_results(output_nanopore, sample, plasmidfinderdir, mergedir):
	print("Start Merging plasmidfinder info")
	print("--------------------------------")
	parent_dir=output_nanopore
	plasmidfinder_file=parent_dir + "/plasmidfinder/results_tab.tsv"
	plasmid_list=[]
	for line in open (plasmidfinder_file, "r"):
		if not line.startswith("Database"):
			fields=line.split("\t")
			plasmid=fields[1]
			plasmid_list.append(plasmid)

	file_lista=open(mergedir+"/plasmidfinder_results_list.txt", "a")
	file_lista.write(sample+"\t"+", ".join(plasmid_list)+"\n")

	#Anadir datos de plasmidos en el fichero de plasmidos de todas las cepas.
	if os.path.isfile(mergedir + "/plasmidfinder_results.xlsx"):	
		df=pd.read_excel(mergedir + "/plasmidfinder_results.xlsx", index_col=0)
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
		file.to_excel(mergedir+ "/plasmidfinder_results.xlsx")

	else:
		data={"sample":[sample]}
		for plasmid in plasmid_list:
			data[plasmid]=["1"]

		file=pd.DataFrame(data)
		print(file)
		file.to_excel(mergedir+ "/plasmidfinder_results.xlsx", index=False)


def mlst(mlstdir, assembly):
	print("Start mlst")
	print("----------")
	if os.path.isfile(mlstdir + "/mlst_results.tab"):
		cmd="/home/admingenomica/software/miniconda3/pkgs/mlst-2.11-0/bin/mlst " + assembly + " >>" + mlstdir + "/mlst_results.tab"
	else:
		cmd="/home/admingenomica/software/miniconda3/pkgs/mlst-2.11-0/bin/mlst " + assembly + " >" + mlstdir + "/mlst_results.tab"
	print(cmd)
	os.system(cmd)


def sistr(output_nanopore,sample,mergedir, assembly):
	print("Start sistr")
	print("----------")
	#sistr --qc -vv --alleles-output ./sistr/1215_22-01_allele-results.json --novel-alleles ./sistr/1215_22-01_novel-alleles.fasta --cgmlst-profiles ./sistr/1215_22-01_cgmlst-profiles.csv -f tab -o ./sistr/1215_22-01_sistr-output.tab ./assembling/1215_22-01.fasta
	parent_dir=output_nanopore
	sistr_dir=os.path.join(parent_dir, "sistr")
	
	cmd="sistr --qc -vv --alleles-output " + sistr_dir + "/" + sample + "_allele-results.json --novel-alleles " + sistr_dir + "/" + sample + "_novel-alleles.fasta --cgmlst-profiles " + sistr_dir+ "/" + sample + "_cgmlst-profiles.csv -f tab -o " + sistr_dir + "/" + sample + "_sistr-output.tab " + assembly
	print(cmd)
	os.system(cmd)
	
	file_sistr_merge=open(mergedir + "/sistr_results.tab", "a")

	#if os.path.isfile(sistrdir + "/sistr_results.tab"):
	for line in open(sistr_dir + "/" + sample + "_sistr-output.tab", "r"):
		if not line.startswith("cgmlst_ST"):
			file_sistr_merge.write(line)
	#else:
	#	cmd="cat " + sistr_dir + "/" + sample + "_sistr-output.tab >>" + sistrdir + "/sistr_results.tab"
	#print(cmd)
	#os.system(cmd)


def resfinder(output_nanopore, assembly, specie="nan"):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Resfinder")
	print("---------------")
	parent_dir=output_nanopore
	resfinder_dir=os.path.join(parent_dir, "resfinder")
	if (specie=="nan"):
		cmd="python /home/admingenomica/software/resfinder/run_resfinder.py -ifa " + assembly + " -o " + resfinder_dir + " -acq -db_res /home/admingenomica/software/resfinder/db_resfinder_10_06_22" 
	else:
		if specie == "escherichia_coli":
			specie = "ecoli"
		cmd="python /home/admingenomica/software/resfinder/run_resfinder.py -ifa " + assembly + " -o " + resfinder_dir + " -acq -db_res /home/admingenomica/software/resfinder/db_resfinder_10_06_22 -s " + specie 
	print(cmd)
	os.system(cmd)




##################################################################################################################################################################
##################################################################################################################################################################

def main():
	#pasar el fichero conel directorio y el nombre de las muestras y recorrer cada muestras antes de hacer el analisis filogenetico.

	#Initiate the parser

	help_text="Pipeline for NGS analysis of Nanopore long reads from microrganism samples. This includes the following steps: raw reads quality analysis and trimming, assembling, genome annotation, mapping with reference genome, variant calling, identification of arm genes..."

	parser = argparse.ArgumentParser(description=help_text)

	#Add arguments
	#argumento sin valor
	#parser.add_argument("-V", "--version", help="show program version", action="store_true")

	required=parser.add_argument_group('required arguments')

	required.add_argument("--sample", "-s", help="sample name or id", required=True)
	required.add_argument("--output_illumina", "-oI", help="output directory for Illumina reads", required=True)
	required.add_argument("--output_nanopore", "-oN", help="output directory for MinIon reads", required=True)
	parser.add_argument("--threads", "-t", help="number of threads")
	required.add_argument("--reads1", "-r1", help="fastq file with reads1 (paired-end)", required=True)
	required.add_argument("--reads2", "-r2", help="fastq file with reads2 (paired-end)", required=True)
	required.add_argument("--reference", "-Ref", help="reference genome (fasta file)") #En caso de no tener referencia, comentar todos los pasos del alineamiento con referencia y variant calling y ejecutar quast sin referencia (comentar linea correspondiente)
	required.add_argument("--mergeresdir", "-md", help="output directory for merged results", required=True)
	required.add_argument("--minion_read", "-mr", help="fastq file nanopore reads", required=True)
	required.add_argument("--specie", "-sp", help="salmonella, escherichia_coli, campylobacter or klebsiella")
	

	args=parser.parse_args()

	if args.sample:
		print("sample = %s" % args.sample)
		sample=args.sample

	if args.output_illumina:
		print("output_illumnina = %s" % args.output_illumina)
		output_illumina=args.output_illumina

	if args.output_nanopore:
		print("output_nanopore = %s" % args.output_nanopore)
		output_nanopore=args.output_nanopore

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
	
	if args.minion_read:
		print("minion_read = %s" % args.minion_read)
		minion_read=args.minion_read

	if args.specie:
		print("specie = %s" % args.specie)
		specie=args.specie
		specie=specie.rstrip()
	
	#create subdirectories
	createSubdirectories(output_nanopore)
	
    #longqc
	longqc(output_nanopore, minion_read, sample, threads)  #Da problemas con la ruta del output
	
    #fastqc
	fastqc(output_nanopore, minion_read, sample, threads)
	
	#Trimming with porechop
	porechop(output_nanopore, sample, threads)
	
	#Filtering with nanofilt
	nanofilt(output_nanopore, output_nanopore+"/primario/fastqc/"+sample+"_trimmed.fastq.gz", sample)
	
	#fasqc trimming
	fasqc_trimming(output_nanopore, output_nanopore+"/primario/fastqc/"+sample+"_trimmed_filtered.fastq.gz", sample, threads)
	
	#longqc trimming
	#longqc_timming(output_nanopore, output_nanopore+"/primario/fastqc/"+sample+"_trimmed_filtered.fastq.gz", sample, threads)
	
	#Ensamblado hibrido
	assembl_hibrid(output_illumina, output_nanopore, output_illumina+"/primario/"+sample +"_trimmed_1U.fastq.gz", output_illumina+"/primario/"+sample+"_trimmed_2U.fastq.gz", sample, output_illumina+"/primario/"+sample+"_trimmed_1P.fastq.gz", output_illumina+"/primario/"+sample+"_trimmed_2P.fastq.gz", output_nanopore+"/primario/fastqc/"+ sample+"_trimmed_filtered.fastq.gz", threads)
	
	#Calidad ensamblado
	quast(output_illumina, output_nanopore, sample, output_illumina+"/primario/"+sample+"_trimmed_1P.fastq.gz", output_illumina+"/primario/"+sample+"_trimmed_2P.fastq.gz", output_illumina + "/primario/" + sample + "trimmed_allU.fastq.gz", output_nanopore+"/primario/fastqc/"+ sample+"_trimmed_filtered.fastq.gz", threads)
	
	#merge assembly quality
	coverage_assembly(output_nanopore, sample, mergeresdir)
	
	
	plasmidfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta")
	
	merge_plasmidfinder_results(output_nanopore, sample, output_nanopore+"/plasmidfinder", mergeresdir)
	

	mlst(mergeresdir, output_nanopore+"/assembling/assembly.fasta")
	
	#run sistr
	if (args.specie) and (specie == "salmonella"):
		sistr(output_nanopore, sample, mergeresdir, output_nanopore+"/assembling/assembly.fasta")
	
	#run resfinder
	if (args.specie) and (specie != "nan"):
		#estaaaaresfinder(output_nanopore, output_nanopore+"/assembling/"+sample+".fasta", specie)
		resfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta", specie)
	else:
		#resfinder(output_nanopore, output_nanopore+"/assembling/"+sample+".fasta")
		resfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta")
	
	

if __name__ == "__main__":

    main()