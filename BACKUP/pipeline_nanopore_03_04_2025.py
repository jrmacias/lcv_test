import sys
import subprocess
import argparse
import os
import pandas as pd

## Nanopore desde fastq
def createSubdirectories(output_nanopore):
	parent_dir=output_nanopore
	
	primary_dir=os.path.join(parent_dir,"primario")
	secundary_dir=os.path.join(parent_dir, "secundario")
	terciary_dir=os.path.join(parent_dir, "terciario")
	assembling_dir=os.path.join(parent_dir, "assembling")
	plasmidfinder_dir=os.path.join(parent_dir, "plasmidfinder")
	resfinder_dir=os.path.join(parent_dir, "resfinder")
	sistr_dir=os.path.join(parent_dir, "sistr")
	pointfinder_dir=os.path.join(parent_dir, "pointfinder")
	virfactor_dir=os.path.join(parent_dir, "virulence_factors")
	prokka_dir=os.path.join(parent_dir, "annotation")
	
	if not os.path.exists(primary_dir):
		os.mkdir(primary_dir)
	if not os.path.exists(secundary_dir):
		os.mkdir(secundary_dir)
	if not os.path.exists(terciary_dir):
		os.mkdir(terciary_dir)
	if not os.path.exists(assembling_dir):
		os.mkdir(assembling_dir)
	if not os.path.exists(plasmidfinder_dir):
		os.mkdir(plasmidfinder_dir)
	if not os.path.exists(resfinder_dir):
		os.mkdir(resfinder_dir)
	if not os.path.exists(sistr_dir):
		os.mkdir(sistr_dir)
	if not os.path.exists(pointfinder_dir):
		os.mkdir(pointfinder_dir)
	if not os.path.exists(virfactor_dir):
		os.mkdir(virfactor_dir)
	if not os.path.exists(prokka_dir):
		os.mkdir(prokka_dir)

def longqc(output_nanopore, minion_read, sample, threads):
	parent_dir=output_nanopore
	out_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd2= "python /home/admingenomica/software/LongQC/longQC.py sampleqc -o " + out_dir + "/longqc_raw_reads" + " -x ont-ligation -s " + sample + " -p " + threads + s + minion_read
	print(cmd2)
	os.system(cmd2)

def nanoplot(output_nanopore, minion_read, sample, threads):
	parent_dir=output_nanopore
	out_dir=os.path.join(parent_dir, "primario")
	s=" "
	nanoplot_dir=out_dir + "/nanoplot"
	cmd1="mkdir" + s + nanoplot_dir
	print(cmd1)
	os.system(cmd1)

	# Replace 'my_environment' with the name of your Conda environment
	environment_name = 'nanoplot'
	# Command to activate the Conda environment
	activate_command = f'conda run -n {environment_name}'
	# Command to execute your Python script
	python_script = '/labgenomica/scripts/pruebas_conda/ejecutar_nanoplot.py ' + minion_read + " " + nanoplot_dir + " " + sample + " " + threads
	# Construct the final command to run
	final_command = f'{activate_command} python {python_script}'
	# Execute the command
	subprocess.call(final_command, shell=True)
	cmd2= "python /home/admingenomica/software/LongQC/longQC.py sampleqc -o " + out_dir + "/longqc_raw_reads" + " -x ont-ligation -s " + sample + " -p " + threads + s + minion_read
	print(cmd2)
	os.system(cmd2)

def fastqc(output_nanopore, minion_read, sample, threads):
	parent_dir=output_nanopore
	fastqc_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="mkdir" + s + fastqc_dir + "/fastqc"
	print(cmd1)
	os.system(cmd1)
	cmd2="/home/admingenomica/software/FastQC-0.12.1/fastqc -t " + threads + " -o " + fastqc_dir + "/fastqc/" + s + minion_read
	print(cmd2)
	os.system(cmd2)

def porechop(output_nanopore, sample, minion_read, threads):
	parent_dir=output_nanopore
	porechop_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="porechop -i" + s + minion_read + s +  "-t " + threads + s + "-v 1 -o" + s + porechop_dir + "/" + sample + "_trimmed.fastq.gz"
	print(cmd1)
	os.system(cmd1)

def nanofilt(output_nanopore, trimming, sample):
	parent_dir=output_nanopore
	nanofilt_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="gunzip -c " + nanofilt_dir + "/" + sample + "_trimmed.fastq.gz "
	cmd2="| NanoFilt -q 7 -l 300 "
	cmd3="| gzip > " + nanofilt_dir + "/" + sample + "_trimmed_filtered.fastq.gz"
	cmd=cmd1+cmd2+cmd3
	print(cmd)
	os.system(cmd)

def fasqc_trimming(output_nanopore, filtered_trimming, sample, threads):
	parent_dir=output_nanopore
	fastqc_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd1="/home/admingenomica/software/FastQC-0.12.1/fastqc -t " + threads + " -o " + fastqc_dir + "/fastqc/" + s + fastqc_dir + "/" + sample + "_trimmed_filtered.fastq.gz"
	print(cmd1)
	os.system(cmd1)


def longqc_trimming(output_nanopore, filtered_trimming, sample, threads):
	parent_dir=output_nanopore
	longqc_dir=os.path.join(parent_dir, "primario")
	s=" "
	cmd="python /home/admingenomica/software/LongQC/longQC.py sampleqc -o " + longqc_dir + "/longQC_trimmed_reads" + s + "-x ont-ligation -s " + sample + s + "-p " + threads + s + filtered_trimming

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
	#cmd1="gzip -d" + s + filtered_trimming
	#print(cmd1)
	#os.system(cmd1)
	#con unpairs sin concatenados
	#cmd2="unicycler -1 " + pair1 + s + "-2 " + pair2 + s + "-s " + unpair1 + " -s " + unpair2 + s + "-l " + filtered_trimming + s + "-o" + s + assembl_hibrid + s + "-t " + threads
	#con unpairs concatenados (mejor opcion)
	cmd2="unicycler -1 " + pair1 + s + "-2 " + pair2 + s + "-s " + output_illumina+"/primario/"+sample+"_trimmed_allU.fastq.gz" + s + "-l " + filtered_trimming + s + "-o" + s + assembling_dir + s + "-t " + threads
	#sin unpairs
	#cmd2="unicycler -1 " + pair1 + s + "-2 " + pair2 + s + "-l " + filtered_trimming + s + "-o" + s + assembl_hibrid + s + "-t " + threads
	print(cmd2)
	os.system(cmd2)
	cmd3="gzip -q" + s + output_nanopore+"/primario/fastqc/"+ sample+"_trimmed_filtered.fastq"
	print(cmd3)
	os.system(cmd3)
	
def assembl_nanopore(output_nanopore, sample, filtered_trimming, threads):
	print("-------------------------")
	print("Start assembling Nanopore")
	print("-------------------------")
	parent_dir_Nanopore=output_nanopore
	assembling_dir=os.path.join(parent_dir_Nanopore, "assembling")
	s=" "
	#flye -g 200k --nano-raw ./Sample-06-X-2022_monkeypox.fastq.gz --out-dir ./assembling_flye --threads 24
	cmd2="flye" + s + "--nano-raw " + filtered_trimming + s + "-o" + s + assembling_dir + s + "-t " + threads
	print(cmd2)
	os.system(cmd2)

def quast(output_illumina, output_nanopore, sample, pair1, pair2, unpair_concatened, filtered_trimming, assembly, threads ):
	print("----------------------")
	print("Start Quast ")
	print("----------------------")
	
	parent_dir_Illumina=output_illumina
	primario_dir_illumina=os.path.join(parent_dir_Illumina, "primario")
	parent_dir_Nanopore=output_nanopore
	assembling_dir=os.path.join(parent_dir_Nanopore, "assembling")
	s=" "

	cmd="quast.py -o " + assembling_dir + s + "-t "+ threads + " --pe1 " + pair1 + s + "--pe2 " + pair2 + s + "--single " + output_illumina+"/primario/"+sample+"_trimmed_allU.fastq.gz --nanopore " + filtered_trimming + s + assembly
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

	if os.path.isfile(coveragedir+ "/assembling_coverage_analysis_results.txt"):	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "a")
		file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) +"\n")
	
	else:	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "w")
		file2.write("sample\ttotal_contigs\tn_contigs>=500\tcontigs>=10000bp\tcontigs>=25000bp\tcontigs>=50000bp\tlargest_contig\tassembly_length\tN50\tL50\tGC\tn_total_reads\t%_mapped_reads\tquast_depth\n")
		file2.write(sample + "\t" + str(n_contigs_total) + "\t" + str(n_contigs) + "\t" + str(n_contigs_10000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(largest_contig) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(gc_content) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth) + "\n")
		

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

#Es obligatorio la especie
def pointfinder(dir, scaffold, specie):
	print("Start pointfinder")
	print("-----------------")
	parent_dir=dir
	#print(scaffold,"\n",specie,"\n")
	pointfinder_dir=os.path.join(parent_dir, "pointfinder")
	cmd="python /home/admingenomica/software/pointfinder/PointFinder.py -i " + scaffold + " -o " + pointfinder_dir + " -s " + specie + " -p /home/admingenomica/software/pointfinder_db -m blastn -t 0.90 -l 0.60 -m_p /home/admingenomica/software/ncbi-blast-2.9.0+/bin/blastn"
	
	print(cmd)
	os.system(cmd)

def merge_pointfinder_results(dir, sample, pointfinderdir):
	print("Start Merging pointfinder info")
	print("--------------------------------")
	parent_dir=dir
	sample_ed=sample.split("_")[0]
	sample="_"+sample
	pointfinder_file=parent_dir + "/pointfinder/assembly_blastn_results.tsv"
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

def resistencias(output_nanopore, sample, resistancedir, specie="nan"):
	print("Start Merging resistance info")
	print("-----------------------------")
	parent_dir=output_nanopore
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

	#Anadir datos de resistencias en el fichero de resistencias de todas las cepas.
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


def abricate(output_nanopore, sample, scaffold):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Abricate - virulence factors")
	print("----------------------------------")
	parent_dir=output_nanopore
	abricate_dir=os.path.join(parent_dir, "virulence_factors")
	'''
	#abricate --db vfdb --quiet ./assembling/contigs.fast
	cmd="/home/admingenomica/software/miniconda3/pkgs/abricate-0.8-pl526ha92aebf_1/bin/abricate --db vfdb --quiet --minid 95 --mincov 95 " + scaffold + " > " + abricate_dir + "/" + sample + "_vFactors.tsv"
	print(cmd)
	os.system(cmd)
	'''
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


def merge_abricate_results(output_nanopore, sample, abricatedir):
	print("Start Merging abricate virulence factors info")
	print("--------------------------------")
	parent_dir=output_nanopore
	abricate_file=parent_dir + "/virulence_factors/" + sample + "_vFactors.tsv"
	vfactors_list=[]
	for line in open (abricate_file, "r"):
		if not line.startswith("#"):
			fields=line.split("\t")
			vfactor=fields[4]
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


def minimap(output_nanopore, sample, filtered_trimming, reference, threads):
	print("Start mapping with minimap")
	print("--------------------------")
	parent_dir=output_nanopore
	secundario_dir=os.path.join(parent_dir, "secundario")
	#/home/admingenomica/software/minimap2-2.22_x64-linux/minimap2 -ax map-ont reference/lambda_phage.mmi ./run1_lambda/nanopore/run1_lambda.fastq.gz -t 10 -o ./run1_lambda/nanopore/secundario/run1_Lambda.sam
	cmd="/home/admingenomica/software/minimap2-2.28_x64-linux/minimap2 -ax map-ont -R \"@RG\\tID:"+ sample +"\\tSM:"+ sample +"\\tPL:NANOPORE\" " + reference + " " + filtered_trimming + " -t " + threads + " -o " + secundario_dir + "/" + sample + ".sam"
	print(cmd)
	os.system(cmd)

	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools view -Sbh " + secundario_dir + "/" + sample + ".sam -o " + secundario_dir + "/" + sample + ".bam"
	print(cmd)
	os.system(cmd)

	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools sort " + secundario_dir + "/" + sample + ".bam -o " + secundario_dir + "/" + sample + "_sort.bam"
	print(cmd)
	os.system(cmd)

	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools index " + secundario_dir + "/" + sample + "_sort.bam"
	print(cmd)
	os.system(cmd)


def qualimap(output_nanopore, sample):
	parent_dir=output_nanopore
	secondary_dir=os.path.join(parent_dir, "secundario")
	#secondary_dir=os.path.join(parent_dir, "secundario_ref_pOX_21_1088-10")
	cmd="qualimap bamqc -sd -bam " + secondary_dir + "/" + sample + "_sort.bam"
	#cmd="qualimap bamqc -sd -bam " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	print(cmd)
	os.system(cmd)


def samtools_coverage(output_nanopore, sample, coveragedir):
	parent_dir=output_nanopore
	secondary_dir=os.path.join(parent_dir, "secundario")
	cmd="/home/admingenomica/software/samtools-1.21/bin/samtools coverage " + secondary_dir + "/" + sample + "_sort.bam > " + secondary_dir + "/" + sample + "_samtools_coverage.tab"
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


def coverage_alignment(output_nanopore, sample, coveragedir):
	print("Start merging coverage alignment results")
	print("-----------------------------------------")
	parent_dir=output_nanopore
	secundario_dir=os.path.join(parent_dir, "secundario")
	ref_length=0
	ass_length=0
	n50=0
	hay_dup="no"
	for line in open(secundario_dir+"/"+sample+"_sort_stats/genome_results.txt", "r"):
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


def bcftools_varcalling(output_nanopore, sample, bam, reference, threads="4"):
	print("Start bcftools_variantcalling")
	print("-----------------------------")
	#bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
	parent_dir=output_nanopore
	terciary_dir=os.path.join(parent_dir, "terciario")
	reference=reference[:-4] + ".fasta"
	print("la referencia es ", reference)

	s=" "
	
	#cmd1="/home/admingenomica/software/bcftools_1.21/bcftools mpileup --min-MQ 50 --min-BQ 30 -I --threads " + threads + " -s " + sample + s + " --annotate INFO/AD -f " + reference + s + bam + s + "| /home/admingenomica/software/bcftools_1.21/bcftools call -mv --ploidy 1 --threads " + threads + " -Ov -o " + terciary_dir + "/" + sample + "_raw.vcf"
	cmd1="/home/admingenomica/software/bcftools_1.21/bcftools mpileup --min-MQ 40 --min-BQ 30 -I --threads " + threads + " -s " + sample + s + " --annotate INFO/AD -f " + reference + s + bam + s + "| /home/admingenomica/software/bcftools_1.21/bcftools call -mv --ploidy 1 --threads " + threads + " -Ov -o " + terciary_dir + "/" + sample + "_raw.vcf"

	print(cmd1)
	os.system(cmd1)
	
	#cmd2="/home/admingenomica/software/bcftools_1.21/bcftools norm -Ov -f " + reference + " -d all -o " + terciary_dir + "/" + sample + "_norm_raw.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#print(cmd2)
	#os.system(cmd2)
	
	cmd3="/home/admingenomica/software/bcftools_1.21/bcftools filter -Ov -e 'QUAL<100||PL[0:0]<150||DP<5' --threads " + threads + " -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#cmd3="/home/admingenomica/software/bcftools_1.21/bcftools filter -Ov -e 'QUAL<100||PL[0:0]<130||DP<5' --threads " + threads + " -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#cmd3="/home/admingenomica/software/bcftools_1.21/bcftools filter -Ov -e 'QUAL<30 || DP<10' -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_norm_raw.vcf"
	print(cmd3)
	os.system(cmd3)
	
	#cmd4="/home/admingenomica/software/bcftools_1.21/bcftools view -O v -V indels -o " + terciary_dir + "/" + sample + "_filter.vcf " + parent_dir + "/" + sample + "_AE017334_q30.vcf"
	#print(cmd4)
	#os.system(cmd4)


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
	required.add_argument("--mergeresdir", "-md", help="output directory for merged results", required=True)
	required.add_argument("--minion_read", "-mr", help="fastq file nanopore reads", required=True)
	required.add_argument("--assembly", "-as", help="hybrid or only_nanopore", required=True)
	required.add_argument("--specie", "-sp", help="salmonella, escherichia_coli, campylobacter or klebsiella")
	required.add_argument("--reference", "-Ref", help="reference genome (fasta file)") #En caso de no tener referencia, comentar todos los pasos del alineamiento con referencia y variant calling y ejecutar quast sin referencia (comentar linea correspondiente)

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

	if args.reference:
		print("reference = %s" % args.reference)
		reference=args.reference

	if args.assembly:
		print("assembly = %s" % args.assembly)
		assembly=args.assembly

	
	#create subdirectories
	createSubdirectories(output_nanopore)
	
    #longqc
	longqc(output_nanopore, minion_read, sample, threads)  #Da problemas con la ruta del output
	
	#nanoplot
	nanoplot(output_nanopore, minion_read, sample, threads)
	
    #fastqc
	fastqc(output_nanopore, minion_read, sample, threads)
	
	#Trimming with porechop
	porechop(output_nanopore, sample, minion_read, threads)
	
	#Filtering with nanofilt
	nanofilt(output_nanopore, output_nanopore+"/primario/"+sample+"_trimmed.fastq.gz", sample)
	
	#fasqc trimming
	fasqc_trimming(output_nanopore, output_nanopore+"/primario/"+sample+"_trimmed_filtered.fastq.gz", sample, threads)
	
	#longqc trimming
	longqc_trimming(output_nanopore, output_nanopore+"/primario/"+sample+"_trimmed_filtered.fastq.gz", sample, threads)
	

	if assembly == "hybrid":
		#Ensamblado hibrido
		assembl_hibrid(output_illumina, output_nanopore, output_illumina+"/primario/"+sample +"_trimmed_1U.fastq.gz", output_illumina+"/primario/"+sample+"_trimmed_2U.fastq.gz", sample, output_illumina+"/primario/"+sample+"_trimmed_1P.fastq.gz", output_illumina+"/primario/"+sample+"_trimmed_2P.fastq.gz", output_nanopore+"/primario/"+ sample+"_trimmed_filtered.fastq.gz", threads)
	elif assembly == "only_nanopore":
		#Ensamblado solo Nanopore
		assembl_nanopore(output_nanopore, sample, output_nanopore+"/primario/"+ sample+"_trimmed_filtered.fastq.gz", threads)
	
	#Calidad ensamblado
	quast(output_illumina, output_nanopore, sample, output_illumina+"/primario/"+sample+"_trimmed_1P.fastq.gz", output_illumina+"/primario/"+sample+"_trimmed_2P.fastq.gz", output_illumina + "/primario/" + sample + "trimmed_allU.fastq.gz", output_nanopore+"/primario/"+ sample+"_trimmed_filtered.fastq.gz", threads)
	
	#Prokka
	prokka(output_nanopore, sample, output_nanopore+"/assembling/assembly.fasta", threads)

	
	#merge assembly quality
	coverage_assembly(output_nanopore, sample, mergeresdir)
	
	
	plasmidfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta")
	
	merge_plasmidfinder_results(output_nanopore, sample, output_nanopore+"/plasmidfinder", mergeresdir)
	
	if (args.specie) and (specie != "nan"):
		pointfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta", specie)
		#pointfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta", specie)

		#run merge_pointfinder_results	
		merge_pointfinder_results(output_nanopore, sample, mergeresdir)
	
	
	#mlst(mergeresdir, output_nanopore+"/assembling/assembly.fasta")
	mlst(mergeresdir, output_nanopore+"/assembling/"+ sample+".fasta")
	
	#run sistr
	if (args.specie) and (specie == "salmonella"):
		sistr(output_nanopore, sample, mergeresdir, output_nanopore+"/assembling/assembly.fasta")
	
	#run resfinder
	if (args.specie) and (specie != "nan"):
		resfinder(output_nanopore, output_nanopore+"/assembling/"+sample+".fasta", specie)
		#resfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta", specie)
	else:
		#resfinder(output_nanopore, output_nanopore+"/assembling/"+sample+".fasta")
		resfinder(output_nanopore, output_nanopore+"/assembling/assembly.fasta")
	
	#run merge resistences info
	if (args.specie):
		resistencias(output_nanopore, sample, mergeresdir, specie)
	else: 
		resistencias(output_nanopore, sample, mergeresdir)
	
	#run abricate
	abricate(output_nanopore, sample, output_nanopore+"/assembling/assembly.fasta")
	#run_merge_abricate_results
	merge_abricate_results(output_nanopore, sample, mergeresdir)

	
	if args.reference:
		#mapping
		minimap(output_nanopore, sample, output_nanopore+"/primario/"+ sample+"_trimmed_filtered.fastq.gz", reference, threads)
		
		#run qualimap
		qualimap(output_nanopore, sample)

		#run variantcalling
		bcftools_varcalling(output_nanopore, sample, output_nanopore + "/secundario/" + sample + "_sort.bam", reference, threads)

		#run coverage_alignment
		coverage_alignment(output_nanopore, sample, mergeresdir)	
		
	
	
if __name__ == "__main__":

    main()