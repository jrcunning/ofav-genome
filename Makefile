all: gag_out/genome.stats

# Generate annotation statistics using GAG
gag_out/genome.stats: ofav.maker.output/ofav.all.renamed.function.gff
	python ~/local/GAG/gag.py --fasta data/ofav_genome.fa --gff ofav.maker.output/ofav.all.gff --fix_start_stop --out gag_out

# Add putative gene function annotations
ofav.maker.output/ofav.all.renamed.function.gff: ofav.maker.output/blastp.output.renamed
	maker_functional_gff data/ref/uniprot_sprot.fasta ofav.maker.output/blastp.output.renamed ofav.maker.output/ofav.all.renamed.gff > ofav.maker.output/ofav.all.renamed.function.gff
	maker_functional_fasta data/ref/uniprot_sprot.fasta ofav.maker.output/blastp.output.renamed ofav.maker.output/ofav.all.maker.proteins.renamed.fasta > ofav.maker.output/ofav.all.maker.proteins.renamed.function.fasta
	maker_functional_fasta data/ref/uniprot_sprot.fasta ofav.maker.output/blastp.output.renamed ofav.maker.output/ofav.all.maker.transcripts.renamed.fasta > ofav.maker.output/ofav.all.maker.transcripts.renamed.function.fasta


# Rename genes in output
ofav.maker.output/blastp.output.renamed: ofav.maker.output/blastp.output
	maker_map_ids --prefix ofav_ --justify 8 ofav.maker.output/ofav.all.gff > ofav.maker.output/ofav.all.id.map
	cp ofav.maker.output/ofav.all.gff ofav.maker.output/ofav.all.renamed.gff
	cp ofav.maker.output/ofav.all.maker.proteins.fasta ofav.maker.output/ofav.all.maker.proteins.renamed.fasta
	cp ofav.maker.output/ofav.all.maker.transcripts.fasta ofav.maker.output/ofav.all.maker.transcripts.renamed.fasta
	cp ofav.maker.output/blastp.output ofav.maker.output/blastp.output.renamed
	map_gff_ids ofav.maker.output/ofav.all.id.map ofav.maker.output/ofav.all.renamed.gff
	map_fasta_ids ofav.maker.output/ofav.all.id.map ofav.maker.output/ofav.all.maker.proteins.renamed.fasta
	map_fasta_ids ofav.maker.output/ofav.all.id.map ofav.maker.output/ofav.all.maker.transcripts.renamed.fasta
	map_data_ids ofav.maker.output/ofav.all.id.map ofav.maker.output/blastp.output.renamed

# BLAST maker proteins against uniprot-sprot
ofav.maker.output/blastp.output: ofav.maker.output/ofav.all.gff
	cd data/ref && makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot.db
	blastp -db data/ref/uniprot_sprot.db -query ofav.maker.output/ofav.all.maker.proteins.fasta -outfmt 6 -num_threads 96 -out ofav.maker.output/blastp.output

# Re-run MAKER using SNAP HMM file, collect results
ofav.maker.output/ofav.all.gff: snap/ofav.hmm
	mpiexec -mca btl ^openib -n 96 maker -base ofav maker_ctl/1/maker_opts.ctl maker_ctl/1/maker_bopts.ctl maker_ctl/1/maker_exe.ctl -fix_nucleotides
	cd ofav.maker.output && gff3_merge -d ofav_master_datastore_index.log && fasta_merge -d ofav_master_datastore_index.log

# Train SNAP gene predictor based on initial MAKER run output
snap/ofav.hmm: ofav.maker.output/ofav0.all.gff
	rm -rf snap
	mkdir snap
	cd snap && maker2zff ../ofav.maker.output/ofav0.all.gff && \
	fathom -categorize 1000 genome.ann genome.dna && \
	fathom -export 1000 -plus uni.ann uni.dna && \
	forge export.ann export.dna && \
	hmm-assembler.pl ofav . > ofav.hmm

# Run initial MAKER annotation, collect results into ofav0.all.gff
ofav.maker.output/ofav0.all.gff: busco/run_ofav/short_summary_ofav.txt maker_ctl/0/maker_opts.ctl
	mpiexec -mca btl ^openib -n 96 maker -base ofav maker_ctl/0/maker_opts.ctl maker_ctl/0/maker_bopts.ctl maker_ctl/0/maker_exe.ctl -fix_nucleotides
	cd ofav.maker.output && gff3_merge -d ofav_master_datastore_index.log && \
	mv ofav.all.gff ofav0.all.gff


# Run BUSCO on filtered ofav assembly and train Augustus
busco/run_ofav/short_summary_ofav.txt: data/ofav_genome.fa
	cd /scratch/projects/crf/ofav-genome/busco && \
	python ~/local/busco/BUSCO.py -f -c 96 --long -i ../data/ofav_genome.fa -o ofav -l ~/local/busco/metazoa_odb9 -m geno

# Generate contigs and summarize fasta files
data/contigs.fasta.summary: data/ofav_genome.fa
	cat data/ofav_genome.fa | seqkit fx2tab | cut -f 2 | sed -r 's/n+/\n/gi' | cat -n | seqkit tab2fx | seqkit replace -p "(.+)" -r "Contig{nr}" > data/contigs.fasta
	fasta_tool --nt_count --summary data/ofav_genome.fa > data/ofav_genome.fa.summary
	fasta_tool --nt_count --summary data/contigs.fasta > data/contigs.fasta.summary

# Reformat names of scaffolds in Orbicella genome
data/ofav_genome.fa: data/48498_ref_ofav_dov_v1_unplaced.fa
	sed 's/\ .*$//' data/48498_ref_ofav_dov_v1_unplaced.fa > data/ofav_genome.fa
