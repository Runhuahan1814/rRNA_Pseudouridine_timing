## Convert muli_fast5 to single_fast5 files using ont_fast5_api v4.1.3 (https://github.com/nanoporetech/ont_fast5_api)
``` bash
multi_to_single_fast5 \
    --input_path \   # path folder containing multi_read_fast5 files  
    --save_path \    # path to folder where single_read fast5 files will be output  
    --recursive \    # recursively search sub-directories  
    --threads \      # number of CPU threads to use  
```

## Basecalling using Guppy (https://timkahlke.github.io/LongRead_tutorials/BS_G.html)
We used Version 6.5.7 for basecalling of all of our reads:
``` bash
guppy_basecaller \
--input_path ${input} \          # input path    
--save_path ${output} \          # output path      
--c rna_r9.4.1_70bps_hac.cfg  \   # high-accuracy DRS config
--flowcell FLO-MIN106 \            # name of flow cell
--kit SQK-RNA02 \                   # name of kit
--compress_fastq \
--fast5_out \
--recursive \
```

## Mapping to the 35S rRNA reference sequence using minimap2 (https://github.com/lh3/minimap2)
We used Version minimap2 (v2.28) with parameters -ax splice -uf -k14 
``` bash
minimap2 -ax splice -uf -k14 ref.fa direct-rna.fq > aln.sam  # Nanopore Direct RNA-seq
```

## Convert mapping results (in SAM files) into BAM files using SAMtools (https://www.htslib.org/)
``` bash
samtools view -bS aligned.sam -o aligned.bam
samtools sort aligned.bam -o aligned.sorted.bam
samtools index aligned.sorted.bam
```

## Evaluation of sequencing quality using Nanoplot (https://github.com/wdecoster/NanoPlot)
``` bash
NanoPlot --summary sequencing_summary.txt --loglength -o summary-plots-log-transformed
```

## Processing of classified reads from R script output
1. Convert R script output files of read names to usable versions to cluster actual reads
``` bash
tail -c +2 input.txt | tr -d '()"' | tr -d '\n ' | tr -d '\r' | tr ',' '\n' > output_list.txt
```
2. Sort and index resulting clustered .bam files
``` bash
samtools view -bN output_list.txt input.sorted.bam > classified.sorted.bam
samtools index classified.sorted.bam
```
3. Coverage analysis (if desired, can be done pre- and post-clustering)
``` bash
samtools depth -a \
	    -J {name of sorted.bam file} \
		-o {name of output file .tsv}
```
4. Use clustering results to separate read fast5 read files into folders for RNA modification analysis
``` bash
for i in $(output_list.txt); do find /input_directory_fast5_files -type f -name "$i".fast5 -exec cp {} ./output_directory_fast5_files \; ; done
```	

## Analysis of pseudouridines using NanoPsu
Details of the analysis pipeline and model implementation are available on the Nanopore_psU GitHub repository: https://github.com/sihaohuanguc/Nanopore_psU (Huang et al.).

## Analysis of rRNA modifications using SignalAlign
Details can be found [in https://github.com/UCSC-nanopore-cgl/signalAlign
](https://github.com/UCSC-nanopore-cgl/signalAlign.)
