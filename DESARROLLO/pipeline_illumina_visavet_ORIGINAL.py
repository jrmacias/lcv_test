import sys
import subprocess
import argparse
import os
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
	
	pointfinder_dir=os.path.join(parent_dir, "pointfinder")
	abricate_dir=os.path.join(parent_dir, "virulence_factors")
	
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
	
	if not os.path.exists(pointfinder_dir):
		os.mkdir(pointfinder_dir)
	if not os.path.exists(abricate_dir):
		os.mkdir(abricate_dir)
	#if not os.path.exists(card_dir):
	#	os.mkdir(card_dir)
	

def fastqc(dir, sample, threads, reads1, reads2):
	parent_dir=dir
	out_dir=os.path.join(parent_dir, "primario")
	#out_dir=os.path.join(parent_dir, "trimming")
	cmd="fastqc -o " + out_dir + " -t " + threads + " -f fastq " + reads1 + " " + reads2
	print("Start FastQC")
	print("------------")
	print(cmd)
	os.system(cmd)


#def trimmomatic(dir, sample, reads1, reads2, threads="4", phred="phred33", adapters="TruSeq3-PE.fa:2:30:10:2:keepBothReads", leading="25", trailing="25", slidingwindow="4:20", minlen="40"):
def trimmomatic(dir, sample, reads1, reads2, threads="4", phred="phred33", adapters="TruSeq3-PE.fa:2:30:10:2:keepBothReads", leading="30", trailing="30", slidingwindow="4:20", minlen="40"):	
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
	cmd1=" java -jar /usr/share/java/trimmomatic.jar PE -threads " + threads + s + reads1 + s + reads2 + s 
	cmd2=out_trimmed_reads1_paired + s + out_trimmed_reads1_unpaired + s + out_trimmed_reads2_paired + s + out_trimmed_reads2_unpaired + s
	cmd3="ILLUMINACLIP:" + adapters + s + "LEADING:" + leading + " TRAILING:" + trailing + " MINLEN:" + minlen

	#cmd3="ILLUMINACLIP:" + adapters + s + "LEADING:" + leading + " TRAILING:" + trailing + " SLIDINGWINDOW:" + slidingwindow + " MINLEN:" + minlen
	cmd=cmd1 + cmd2 + cmd3
	print(cmd)
	os.system(cmd)
	return(out_trimmed_reads1_paired, out_trimmed_reads2_paired, out_trimmed_reads1_unpaired, out_trimmed_reads2_unpaired)


def kraken(dir, reads1, reads2, threads="4"):
	print("Start kraken")
	print("------------")
	parent_dir=dir
	primary_dir=os.path.join(parent_dir, "primario")
	#/home/vserver1/software/kraken2-2.0.8-beta/kraken2 --db /home/vserver1/software/kraken2-2.0.8-beta/taxoDB --threads 12 --output ./ --gzip-compressed --paired ./AV-01_S32_L001_R1_001.fastq.gz ./AV-01_S32_L001_R2_001.fastq.gz  
	cmd="/home/vserver1/software/kraken2-2.0.8-beta/kraken2 --db /home/vserver1/software/kraken2-2.0.8-beta/taxoDB --threads "+threads+ " --output " + primary_dir + " --gzip-compressed --paired " + reads1 + " " + reads2
	os.system(cmd)



def spades(dir, sample, paired_reads1, paired_reads2, unpaired_reads1, unpaired_reads2, threads="4"):
#def spades(dir, sample, paired_reads1, paired_reads2, unpaired_reads1, threads="4"):
	print("Start SPADES")
	print("------------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	#para correr plasmidspades
	#assembling_dir=os.path.join(parent_dir, "plasmidspades")
	cmd1="/home/vserver1/software/SPAdes-3.14.1/bin/spades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2 + " -s " + unpaired_reads1 + " -s " + unpaired_reads2
	#cmd1="/home/vserver1/software/SPAdes-3.14.1/bin/spades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2
	#para correr plasmidspades:
	#cmd1="/home/vserver1/software/SPAdes-3.14.1/bin/plasmidspades.py -t " + threads + " -1 " + paired_reads1 + " -2 " + paired_reads2 + " --s1 " + unpaired_reads1
	#cmd2=" -o " + assembling_dir
	cmd2=" --cov-cutoff auto" + " -o " + assembling_dir + " --careful"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)


def quast(dir, sample, scaffold, paired_reads1, paired_reads2, unpaired_reads1, unpaired_reads2, reference="nan", threads="4"):
	print("Start QUAST")
	print("-----------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	if not reference == "nan":
		cmd="/home/vserver1/software/quast-5.0.2/quast.py " + " -o " + assembling_dir + " -t " +threads + " -r " + reference + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " --single " + unpaired_reads1 + " --single " + unpaired_reads2 + " " + scaffold
		#cmd="/home/vserver1/software/quast-5.0.2/quast.py " + " -o " + assembling_dir + " -t " +threads + " -r " + reference + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " --single " + unpaired_reads1 + " " + scaffold
	else:
		cmd="/home/vserver1/software/quast-5.0.2/quast.py " + " -o " + assembling_dir + " -t " +threads + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " --single " + unpaired_reads1 + " --single " + unpaired_reads2 + " " + scaffold
		#cmd="/home/vserver1/software/quast-5.0.2/quast.py " + " -o " + assembling_dir + " -t " +threads + " --pe1 " + paired_reads1 + " --pe2 " + paired_reads2 + " " + scaffold

	print(cmd)
	os.system(cmd)



def abacas(dir, sample, scaffold, reference):
	print("Start abacas")
	print("------------")
	parent_dir=dir
	assembling_dir=os.path.join(parent_dir, "assembling")
	cmd="abacas -r " + reference + " -q " + scaffold + " -p nucmer -o " + assembling_dir + "/" + sample + " -m"
	print(cmd)
	os.system(cmd)


def prokka(dir, sample, scaffold, threads="4"):
	#do /home/vserver1/software/prokka/bin/prokka --outdir ./$c/annotation --force --prefix $c"_ann" ./$c/annotation/$c"_assembling.result.fasta"
	print("Start Prokka")
	print("------------")
	parent_dir=dir
	annotation_dir=os.path.join(parent_dir, "annotation")	
	cmd= "prokka --kingdom Bacteria --proteins /home/vserver1/software/resfinder_db/resfinder_all_prot_db.fasta --outdir " + annotation_dir + " --gcode 11 --cpus " + threads + " --prefix " + sample + " " + scaffold
	#cmd= "prokka --kingdom Viruses --proteins /home/vserver1/software/resfinder_db/resfinder_all_prot_db.fasta --outdir " + annotation_dir + " --gcode 11 --cpus " + threads + " --prefix " + sample + " " + scaffold 	
	print(cmd)
	os.system(cmd)	


def mlst(mlstdir, scaffold):
	print("Start mlst")
	print("----------")
	if os.path.isfile(mlstdir + "/mlst_results.tab"):
		cmd="mlst " + scaffold + " >>" + mlstdir + "/mlst_results.tab"
	else:
		cmd="mlst " + scaffold + " >" + mlstdir + "/mlst_results.tab"
	print(cmd)
	os.system(cmd)


def resfinder(dir, scaffold):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Resfinder")
	print("---------------")
	parent_dir=dir
	resfinder_dir=os.path.join(parent_dir, "resfinder")
	cmd="python /home/vserver1/software/resfinder/resfinder.py -i " + scaffold + " -o " + resfinder_dir + " -x -p /home/vserver1/software/resfinder_db -mp /home/vserver1/software/miniconda3/bin/blastn"
	print(cmd)
	os.system(cmd)


def pointfinder(dir, scaffold):
	print("Start pointfinder")
	print("-----------------")
	parent_dir=dir
	pointfinder_dir=os.path.join(parent_dir, "pointfinder")
	#cmd="python /home/vserver1/software/pointfinder/PointFinder.py -i " + scaffold + " -o " + pointfinder_dir + " -s salmonella -p /home/vserver1/software/pointfinder_db -m blastn -t 0.90 -l 0.60 -m_p /home/vserver1/software/miniconda3/bin/blastn"
	#cmd="python /home/vserver1/software/pointfinder/PointFinder.py -i " + scaffold + " -o " + pointfinder_dir + " -s campylobacter -p /home/vserver1/software/pointfinder_db -m blastn -t 0.90 -l 0.60 -m_p /home/vserver1/software/miniconda3/bin/blastn"
	cmd="python /home/vserver1/software/pointfinder/PointFinder.py -i " + scaffold + " -o " + pointfinder_dir + " -s escherichia_coli -p /home/vserver1/software/pointfinder_db -m blastn -t 0.90 -l 0.60 -m_p /home/vserver1/software/miniconda3/bin/blastn"
	print(cmd)
	os.system(cmd)


def abricate(dir, sample, scaffold):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Abricate - virulence factors")
	print("----------------------------------")
	parent_dir=dir
	abricate_dir=os.path.join(parent_dir, "virulence_factors")
	#abricate --db vfdb --quiet ./assembling/contigs.fast
	cmd="abricate --db vfdb --quiet --minid 95 --mincov 95 " + scaffold + " > " + abricate_dir + "/" + sample + "_vFactors.tsv"
	print(cmd)
	os.system(cmd)


def merge_abricate_results(dir, sample, abricatedir):
	print("Start Merging abricate virulence factors info")
	print("--------------------------------")
	parent_dir=dir
	abricate_file=parent_dir + "/virulence_factors/" + sample + "_vFactors.tsv"
	vfactors_list=[]
	for line in open (abricate_file, "r"):
		if not line.startswith("#"):
			fields=line.split("\t")
			vfactor=fields[4]
			vfactors_list.append(vfactor)

	file_lista=open(abricatedir+"/vFactors_results_list.txt", "a")
	file_lista.write(sample+"\t"+", ".join(vfactors_list)+"\n")


	if os.path.isfile(abricatedir + "/vFactors_results.xlsx"):	
		df=pd.read_excel(abricatedir + "/vFactors_results.xlsx", index_col=0)
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
		file.to_excel(abricatedir+ "/vFactors_results.xlsx")

	else:
		data={"sample":[sample]}
		for vfactor in vfactors_list:
			data[vfactor]=["1"]

		file=pd.DataFrame(data)
		#print(file)
		file.to_excel(abricatedir + "/vFactors_results.xlsx", index=False)



def integronfinder(dir, scaffold, threads):
	print("Start integronfinder")
	print("-----------------")
	parent_dir=dir
	integronfinder_dir=os.path.join(parent_dir, "integronfinder")
	#integron_finder --local-max --func-annot --cpu 20 --pdf ./assembly_graph.cycs.fasta
	cmd="integron_finder --local-max --func-annot --cpu " + threads + " --pdf  --outdir " + integronfinder_dir + " " + scaffold
	print(cmd)
	os.system(cmd)



def recycler(dir, sample, assembly_graph):
	print("Start recycler")
	print("-------------------")

	#Correr estos comandos para preparar el fichero bam que necesita recycler antes de lanzarlo.
	#make_fasta_from_fastg.py -g assembly_graph.fastg [-o assembly_graph.nodes.fasta]
	#bwa index assembly_graph.nodes.fasta
	#bwa mem assembly_graph.nodes.fasta R1.fastq.gz R2.fastq.gz | samtools view -buS - > reads_pe.bam
	#samtools view -bF 0x0800 reads_pe.bam > reads_pe_primary.bam
	#samtools sort reads_pe_primary.bam -o reads_pe_primary.sort.bam
	#samtools index reads_pe_primary.sort.bam
	#python2.7 /home/vserver1/software/Recycler/bin/recycle.py -g ./assembling/ZTA20_00785_assembling/spades/assembly_graph.fastg -k 55 -b ./recycler/ZTA20_00785_reads_pe_primary.sort.bam -i True -o ./recycler/

	parent_dir=dir
	recycler_dir=os.path.join(parent_dir, "recycler")
	
	cmd="python2.7 /home/vserver1/software/Recycler/bin/make_fasta_from_fastg.py -g " + assembly_graph + " -o " + recycler_dir + "/" + sample + "_assGraphNodes.fasta"
	print(cmd)
	os.system(cmd)
	cmd2="bwa index " + recycler_dir + "/" + sample + "_assGraphNodes.fasta"
	print(cmd2)
	os.system(cmd2)
	cmd3="bwa mem " + recycler_dir + "/" + sample + "_assGraphNodes.fasta " + parent_dir + "/primario/" + sample + "_trimmed_1P_ed.fastq " + parent_dir + "/primario/" + sample + "_trimmed_2P_ed.fastq | samtools view -buS - > " + recycler_dir + "/" + sample + "_reads_pe.bam"
	print(cmd3)
	os.system(cmd3)
	cmd4="samtools view -bF 0x0800 " + recycler_dir + "/" + sample + "_reads_pe.bam > " + recycler_dir + "/" + sample + "_reads_pe_primary.bam"
	print(cmd4)
	os.system(cmd4)
	cmd5="samtools sort " + recycler_dir + "/" + sample + "_reads_pe_primary.bam -o " + recycler_dir + "/" + sample + "_reads_pe_primary.sort.bam"
	print(cmd5)
	os.system(cmd5)
	cmd6="samtools index " + recycler_dir + "/" + sample + "_reads_pe_primary.sort.bam"
	print(cmd6)
	os.system(cmd6)
	cmd7="rm -rf " + recycler_dir + "/" + sample + "_reads_pe.bam"
	os.system(cmd7)
	cmd8="rm -rf " + recycler_dir + "/" + sample + "_reads_pe_primary.bam"
	os.system(cmd8)
	cmd9="python2.7 /home/vserver1/software/Recycler/bin/recycle.py -g " + assembly_graph + " -k 55 -b " + recycler_dir + "/" + sample + "_reads_pe_primary.sort.bam -i True -o " + recycler_dir
	print(cmd9)
	os.system(cmd9)

def plasmidfinder(dir, scaffold):
	#python3 resfinder.py -i test.fsa -o . -p /path/to/resfinder_db -mp /path/to/blastn -d aminoglycoside -t 0.90 -l 0.60
	print("Start Plasmidfinder")
	print("-------------------")
	parent_dir=dir
	plasmidfinder_dir=os.path.join(parent_dir, "plasmidfinder")
	#/home/vserver1/software/plasmidfinder/plasmidfinder.py -i ./MS3828/hybrid_assembling_2/assembly.fasta -o ./MS3828/plasmids/ -d /home/vserver1/software/plasmidfinder_db/
	cmd="plasmidfinder.py -i " + scaffold + " -o " + plasmidfinder_dir + " -x -p /home/vserver1/software/genomicepidemiology-plasmidfinder_db-1307168b1ce7/" + " -t 0.80"
	print(cmd)
	os.system(cmd)


def merge_plasmidfinder_results(dir, sample, plasmidfinderdir):
	print("Start Merging plasmidfinder info")
	print("--------------------------------")
	parent_dir=dir
	#plasmidfinder_file=parent_dir + "/plasmidfinder/results_tab.tsv"
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


'''Identificación de genes de resistencia con arg_ann y card (LO DESCARTAMOS DEL ANÁLISIS)
def arg_ann(dir, sample, scaffold, threads="4"):
	print("Start argann")
	print("------------")
	parent_dir=dir
	resistencias_dir=os.path.join(parent_dir, "amr_genes")
	arg_dir=os.path.join(resistencias_dir, "arg_ann_db")
	cmd="blastx -query " + scaffold + " -out " + arg_dir + "/" + sample + "_argAnn.txt -db /home/vserver1/software/arg_ann_db/arg_ann_BlastDB -num_threads " + threads + " -evalue 1e-3 -outfmt 6" 
	print(cmd)
	os.system(cmd)
'''

def rgi_card(dir, sample, scaffold):
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


def resistencias(dir, sample, resistancedir):
	print("Start Merging resistance info")
	print("-----------------------------")
	parent_dir=dir
	resfinder_file=parent_dir + "/resfinder/results_tab.tsv"
	#resfinder_file=parent_dir + "/amr_genes/resfinder/results_tab.tsv"


	resistencias_dict={}

	listadegenes=[]
	for line in open (resfinder_file, "r"):
		if not line.startswith("Database"):
			fields=line.split("\t")
			ant=fields[0]
			gene=fields[1]
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


	


def bwa(dir, sample, reads1, reads2, reference, threads="4"):
	print("Start bwa")
	print("---------")
	#bwa mem -t 14 ./data/Reference.fna ./data/illumina/Illumina_R1.fastq.gz ./data/illumina/Illumina_R2.fastq.gz > ./mapping/illumina.bwa.sam
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario")
	s=" "
	cmd1="bwa mem -t " + threads + s + reference + s + reads1 + s + reads2 
	cmd2= " | samtools view -Sbh - > "+ secondary_dir + "/" + sample + "_bwa.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	cmd="samtools sort -@ " + threads + s + secondary_dir + "/" + sample + "_bwa.bam -o " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	cmd="samtools index " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)



'''Alineamiento con bowtie (bwa da mejores resultados)
def bowtie(dir, sample, reads1, reads2, reference, threads="4"):
	print("Start bowtie")
	print("------------")
	#/home/vserver1/software/bowtie2-2.4.2-linux-x86_64/bowtie2
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "secundario_2")
	s=" "
	cmd1="/home/vserver1/software/bowtie2-2.4.2-linux-x86_64/bowtie2 -p " + threads  + s + "--met-file" + s + secondary_dir + "/" + sample + "_bwt_metrics.txt" +  s + "-x" + s + "/home/vserver1/Documentos/salmonellas_enteritidis/ref_enteritidis/s_ent_bwt" + s + "-1" + s + reads1 + s + "-2" + s + reads2
	cmd2= " | samtools view -Sbh - > "+ secondary_dir + "/" + sample + "_bwa.bam"
	cmd=cmd1+cmd2
	print(cmd)
	os.system(cmd)
	#samtools sort MS3825.bam -o MS3825.sorted.bam
	cmd="samtools sort -@" + threads + s + secondary_dir + "/" + sample + "_bwa.bam -o " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
	#samtools index MS3825.sorted.bam
	cmd="samtools index " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	os.system(cmd)
'''

'''Analisis de cobertura del bam (NO SE USA, SE HACE CON QUIALIMAP)
def coverage(dir, sample, reads1, reads2, bam, coveragedir):
	parent_dir=dir
	#assembling_dir=os.path.join(parent_dir, "assembling")
	secondary_dir=os.path.join(parent_dir, "secundario")
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
		elif field=="# contigs":
			n_contigs=int(value)
		elif field=="N50":
			n50=int(value)
		
	n_mapped_reads_with_reference= int(subprocess.check_output("samtools view -F 0x4 " + bam + " | cut -f 1 | sort | uniq | wc -l", shell=True))
	n_total_reads=int(subprocess.check_output("samtools view " + bam + " | cut -f 1 | sort | uniq | wc -l", shell=True))
	read1_length= float(subprocess.check_output("zcat " + reads1 + " | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'", shell=True))
	read2_length= float(subprocess.check_output("zcat " + reads2 + " | awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}'", shell=True))
	ave_read_length=(read1_length+read2_length)/2
	alignment_depth=(n_mapped_reads_with_reference*ave_read_length*2)/ref_length
	assembly_depth=(n_total_reads*ave_read_length*2)/ass_length
	print(sample,"\t",ave_read_length, "\t", str(ave_read_length), "\n")

	if os.path.isfile(coveragedir+ "/alignment_coverage_analysis_results.txt"):	
		file=open(coveragedir + "/alignment_coverage_analysis_results.txt", "a")
		file.write(sample + "\t" + str(ref_length) + "\t" + str(n_mapped_reads_with_reference) + "\t" + str(ave_read_length) + "\t" + str(alignment_depth)+"\n")
	else:	
		file=open(coveragedir + "/alignment_coverage_analysis_results.txt", "w")
		file.write("sample\tref_length\treads_mapped_on_ref\tave_read_length\tdepth\n")
		file.write(sample + "\t" + str(ref_length) + "\t" + str(n_mapped_reads_with_reference) + "\t" + str(ave_read_length) + "\t" + str(alignment_depth)+"\n")
	

	if os.path.isfile(coveragedir+ "/assembling_coverage_analysis_results.txt"):	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "a")
		file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(n_total_reads) + "\t" + str(ave_read_length) + "\t" + str(assembly_depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "w")
		file2.write("sample\tn_contigs\tassembly_length\tN50\tn_total_reads\tave_read_length\tdepth\n")
		file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(n_total_reads) + "\t" + str(ave_read_length) + "\t" + str(assembly_depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
'''

def coverage_assembly(dir, sample, coveragedir):
	print("Start merging coverage assemblies results")
	print("-----------------------------------------")
	parent_dir=dir
	#assembling_dir=os.path.join(parent_dir, "flye")
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
		elif field=="# contigs":
			n_contigs=int(value)
		elif field=="# contigs (>= 1000 bp)":
			n_contigs_1000=int(value)
		elif field=="# contigs (>= 5000 bp)":
			n_contigs_5000=int(value)
		elif field=="# contigs (>= 25000 bp)":
			n_contigs_25000=int(value)	
		elif field=="# contigs (>= 50000 bp)":
			n_contigs_50000=int(value)

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
		file2.write(sample + "\t" + str(n_contigs) + "\t" + str(n_contigs_1000) + "\t" + str(n_contigs_5000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:	
		file2=open(coveragedir + "/assembling_coverage_analysis_results.txt", "w")
		file2.write("sample\tn_contigs\tcontigs>=1000bp\tcontigs>=5000bp\tcontigs>=25000bp\tcontigs>=50000bp\tassembly_length\tN50\tL50\tn_total_reads\t%_mapped_reads\tdepth\n")
		file2.write(sample + "\t" + str(n_contigs) + "\t" + str(n_contigs_1000) + "\t" + str(n_contigs_5000) + "\t" + str(n_contigs_25000) + "\t" + str(n_contigs_50000) + "\t" + str(ass_length) + "\t" + str(n50) + "\t" + str(l50) + "\t" + str(n_total_reads) + "\t" + str(percentage_mapped_reads_with_assembly) + "\t" + str(assembly_depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")


def qualimap(dir, sample):
	parent_dir=dir
	secondary_dir=os.path.join(parent_dir, "alineamiento_plasmid_21_1088-10_pOX2")
	#secondary_dir=os.path.join(parent_dir, "secundario")

	cmd="/home/vserver1/software/qualimap_v2.2.1/qualimap bamqc -bam " + secondary_dir + "/" + sample + "_bwa_sort.bam"
	print(cmd)
	os.system(cmd)


def coverage_alignment(dir, sample, coveragedir):
	print("Start merging coverage alignment results")
	print("-----------------------------------------")
	parent_dir=dir
	secundario_dir=os.path.join(parent_dir, "alineamiento_plasmid_21_1088-10_pOX2")
	ref_length=0
	ass_length=0
	n50=0

	for line in open(secundario_dir+"/"+sample+"_bwa_sort_stats/genome_results.txt", "r"):
		line=line.rstrip()
		if "number of reads = " in line:
			fields_list=line.split(" ")
			nreads=fields_list[-1]
		elif "number of mapped reads = " in line:
			fields_list=line.split(" ")
			nmappedreads=fields_list[-2]
			percentage_mapped_reads=fields_list[-1][1:-1]
		elif "mean coverageData = " in line:
			fields_list=line.split(" ")
			depth=fields_list[-1]
		elif "of reference with a coverageData >= 20X" in line:
			fields_list=line.split(" ")
			coverage=fields_list[-8]


	if os.path.isfile(coveragedir+ "/alignment_coverage_stats.txt"):	
		file2=open(coveragedir + "/alignment_coverage_stats.txt", "a")
		file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(coverage) + "\t" + str(depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")
	
	else:	
		file2=open(coveragedir + "/alignment_coverage_stats.txt", "w")
		file2.write("sample\ttotal_reads\tmapped_reads\thorizontal_coverage\tdepth\n")
		file2.write(sample + "\t" + str(nreads) + "\t" + str(nmappedreads) + "\t" + str(percentage_mapped_reads) + "\t" + str(coverage) + "\t" + str(depth)+"\n")
		#file2.write(sample + "\t" + str(n_contigs) + "\t" + str(ass_length) + "\t" + str(n50) + "\n")




def bcftools_varcalling(dir, sample, bam, reference, threads="4"):
	print("Start bcftools_variantcalling")
	print("-----------------------------")
	#bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
	parent_dir=dir
	terciary_dir=os.path.join(parent_dir, "terciario")
	s=" "
	
	cmd1="bcftools mpileup --min-MQ 30 --min-BQ 30 -I --threads " + threads + " -f " + reference + s + bam + s + "| bcftools call -mv --ploidy 1 --threads " + threads + " -Ov -o " + terciary_dir + "/" + sample + "_raw.vcf"
	print(cmd1)
	os.system(cmd1)
	
	#cmd2="bcftools norm -Ov -f " + reference + " -d all -o " + terciary_dir + "/" + sample + "_norm_raw.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
	#print(cmd2)
	#os.system(cmd2)
	
	cmd3="bcftools filter -Ov -e 'QUAL<50||DP<5' --threads " + threads + " -o " + terciary_dir + "/" + sample + "_filter.vcf " + terciary_dir + "/" + sample + "_raw.vcf"
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
		cmd2="bcftools consensus -p " + sample + " -m " + phages_bed + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"
		print(cmd2)
		os.system(cmd2)
	else:
		cmd2="bcftools consensus -p " + sample + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"
		print(cmd2)
		os.system(cmd2)	
	#cmd2="bcftools consensus -p " + sample + " -m " + phages_bed + " -f " + reference + " " + terciary_dir + "/" + sample + "_filter_Nophages.vcf.gz >" + terciary_dir + "/" + sample + "_consensus.fasta"




def merge_pointfinder_results(dir, sample, pointfinderdir):
	print("Start Merging pointfinder info")
	print("--------------------------------")
	parent_dir=dir
	pointfinder_file=parent_dir + "/pointfinder/scaffolds_blastn_results.tsv"
	#pointfinder_file=parent_dir + "/pointfinder/" + sample+ "_blastn_results.tsv"
	point_list=[]
	for line in open (pointfinder_file, "r"):
		if not line.startswith("Mutation"):
			fields=line.split("\t")
			point=fields[0]+"-"+fields[2] 
			point_list.append(point)



	#Anadir datos de mutaciones puntuales en el fichero de mutaciones puntuales de todas las cepas.
	if os.path.isfile(pointfinderdir + "/pointfinder_results.xlsx"):	
		df=pd.read_excel(pointfinderdir + "/pointfinder_results.xlsx", index_col=0)
		data=df.to_dict()

		for point in point_list:
			if point in data:
				data[point][sample]="1"
			else:
				data[point]={sample:"1"}

		if len(point_list)==0:
			for point in data:
				data[point][sample]=""

		file=pd.DataFrame(data)
		file.to_excel(pointfinderdir+ "/pointfinder_results.xlsx")

	else:
		data={"sample":[sample]}
		for point in point_list:
			data[point]=["1"]

		file=pd.DataFrame(data)
		print(file)
		file.to_excel(pointfinderdir+ "/pointfinder_results.xlsx", index=False)
		
		
		







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
	required.add_argument("--reference", "-R", help="reference genome (fasta file)") #En caso de no tener referencia, comentar todos los pasos del alineamiento con referencia y variant calling y ejecutar quast sin referencia (comentar linea correspondiente)
	required.add_argument("--mergeresdir", "-md", help="output directory for merged results", required=True)
	required.add_argument("--phages_bed", "-pb", help="Bed file with phages coordinates")
	'''
	required.add_argument("--coveragedir", "-cd", help="output directory for coverage analysis results", required=True)
	required.add_argument("--mlstdir", "-md", help="output directory for mlst results", required=True)
	required.add_argument("--resistancesdir", "-rd", help="output directory for gene resistance analysis results", required=True)
	required.add_argument("--plasmidfinderdir", "-pd", help="output directory for plasmidfinder results", required=True)
	required.add_argument("--abricatedir", "-ad", help="output directory for abricate virulence factors results", required=True)
	'''

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

	'''
	#create subdirectories
	createSubdirectories(output)
	
	#run fastqc
	fastqc(output, sample, threads, reads1, reads2)
	
	#run trimmomatic
	trimmed_paired_reads1, trimmed_paired_reads2, trimmed_unpaired_reads1, trimmed_unpaired_reads2=trimmomatic(output, sample, reads1, reads2, threads)

	#run fastqc
	fastqc(output, sample, threads, trimmed_paired_reads1, trimmed_paired_reads2)
	#fastqc(output, sample, threads, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz")
	
	#run spades
	spades(output, sample, trimmed_paired_reads1, trimmed_paired_reads2, trimmed_unpaired_reads1, trimmed_unpaired_reads2, threads)
	#spades(output, sample, output+"/primario/"+sample+"_trimmed_1P_ed.fastq", output+"/primario/"+sample+"_trimmed_2P_ed.fastq", output+"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads)
	#spades(output, sample, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz", output + "/trimming/" + sample + "_trimmed_L1_SE.fastq.gz", threads)
	
	#run quast
	if args.reference:
		quast(output, sample, output+"/assembling/scaffolds.fasta",  output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", output +"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", reference, threads)
	else:
		#sin referencia
		quast(output, sample, output+"/assembling/scaffolds.fasta", output +"/primario/"+sample+"_trimmed_1P.fastq",  output +"/primario/"+sample+"_trimmed_2P.fastq", output+"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads=threads)
		#quast(output, sample, output+"/hybrid_assembly_prueba3/assembly.fasta", output +"/primario/"+sample+"_trimmed_1P.fastq",  output +"/primario/"+sample+"_trimmed_2P.fastq", output+"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads=threads)		
		#quast(output, sample, output+"/assembling_spades_meta/scaffolds.fasta", output +"/primario/"+sample+"_trimmed_1P.fastq.gz",  output +"/primario/"+sample+"_trimmed_2P.fastq.gz", output+"/primario/"+sample+"_trimmed_1U.fastq.gz", output+"/primario/"+sample+"_trimmed_2U.fastq.gz", threads=threads)
	
	#run mlst
	mlst(mergeresdir, output+"/assembling/scaffolds.fasta")
	#mlst(mlstdir, output+"/assembling/"+sample+"_assembling.result.fasta")
	
	#run resfinder
	resfinder(output, output+"/assembling/scaffolds.fasta")
	#resfinder(output, output+"/minimap_miniasm/" + sample +"_assembly.fasta")
	#resfinder(output, output+"/assembling/"+sample+"_assembling.result.fasta")
	#resfinder(output, output+"/hybrid_assembling_prueba2/assembly.fasta")
	#resfinder(output, output+"/"+sample+".fasta")
	#resfinder(output, output+"/assembly_graph.cycs.fasta")
	
	#run abricate
	abricate(output, sample, output+"/assembling/scaffolds.fasta")
	#abricate(output, sample, output+"/minimap_miniasm/" + sample +"_assembly.fasta")
	#abricate(output, sample, output+"/assembling/"+sample+"_assembling.result.fasta")
	
	#run_merge_abricate_results
	merge_abricate_results(output, sample, mergeresdir)
	

	#run arg_ann_db
	#arg_ann(output, sample, output+"/assembling/" + sample + ".contigs_modifiedId.fasta", threads)

	#run rgi-card
	#rgi_card(output, sample, output+"/assembling/scaffolds.fasta")
	#rgi_card(output, sample, output+"/assembling/" + sample + ".contigs_modifiedId.fasta")
	
	#run merge resistences info
	resistencias(output, sample, mergeresdir)
	#Usar resistencias2 cuando activamos rgi-card
	#resistencias2(output, sample, mergeresdir)
	
	#run pointfinder
	pointfinder(output, output+"/assembling/scaffolds.fasta")
	#pointfinder(output, output+"/minimap_miniasm/" + sample +"_assembly.fasta")
	#pointfinder(output, output+"/assembling/" + sample + ".contigs.fasta")
	#pointfinder(output, output+"/assembling/"+sample+"_assembling.result.fasta")
	
	#run recycler
	recycler(output, sample, output+"/assembling/assembly_graph.fastg")
	
	
	#run plasmidfinder
	plasmidfinder(output, output+"/assembling/scaffolds.fasta")
	#plasmidfinder(output, output+"/minimap_miniasm/" + sample +"_assembly.fasta")
	#plasmidfinder(output, output+"/scaffolds.fasta")
	#plasmidfinder(output, output+"/hybrid_assembling_prueba2/assembly.fasta")
	#plasmidfinder(output, output+"/assembling/" + sample + "_assembling.result.fasta")
	#plasmidfinder(output, output+"/assembly_graph.cycs.fasta")
	
	#run merge plasmidfinder results
	merge_plasmidfinder_results(output, sample, mergeresdir)
	
	#run integronfinder
	#integronfinder(output, output+"/assembling/scaffolds.fasta", threads)
	


	#run abacas
	#abacas(output, sample, output+"/assembling/scaffolds.fasta", reference)
	#abacas(output, sample, output+"/assembling/" + sample + ".contigs.fasta", reference)
	
	#run prokka
	prokka(output, sample, output+"/assembling/scaffolds.fasta", threads)
	#prokka(output, sample, output+"/assembling/"+sample+".MULTIFASTA.fa")
	#prokka(output, sample, output+"/hybrid_assembling_prueba2/assembly.fasta", threads)
	'''

	#run kraken
	#kraken(output, trimmed_reads1, trimmed_reads2, threads)
	

	if args.reference:
		'''
		#MAPPING
		#run bwa
		#bwa(output, sample, trimmed_paired_reads1, trimmed_paired_reads2, reference, threads)
		bwa(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", reference, threads)
		#bwa(output, sample, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz", reference, threads)

		#run bowtie
		#bowtie(output, sample, output + "/trimming/" + sample + "_trimmed_L1_R1.fastq.gz",  output + "/trimming/" + sample + "_trimmed_L1_R2.fastq.gz", reference, threads)
		'''
		#run qualimap
		qualimap(output, sample)
		'''
		#run variantcalling
		bcftools_varcalling(output, sample, output + "/secundario/" + sample + "_bwa_sort.bam", reference, threads)
		
		#run bedtools to remove snps in phage coordinates from vcf files
		if args.phages_bed:
			bedtools_filterPhageRegion(output, sample, output + "/terciario/" + sample + "_filter.vcf", phages_bed)
		
		#run bcftools consensus to generate consensus sequence
		if args.phages_bed:
			generate_consensus_sequence(output, sample, output + "/terciario/" + sample + "_filter.vcf", reference, phages_bed)
		else:
			generate_consensus_sequence(output, sample, output + "/terciario/" + sample + "_filter.vcf", reference)
		


		#run coverage
		#coverage(output, sample, output+"/primario/"+sample+"_trimmed_1P.fastq.gz", output+"/primario/"+sample+"_trimmed_2P.fastq.gz", output + "/secundario/" + sample + "_bwa_sort.bam", mergeresdir)
		'''
		
		#run coverage_alignment
		coverage_alignment(output, sample, mergeresdir)
	
	'''
	#run coverage_assembly
	coverage_assembly(output, sample, mergeresdir)
	
	#run merge_pointfinder_results
	merge_pointfinder_results(output, sample, mergeresdir)
	'''




if __name__ == "__main__":

    main()
