## basecalling using Guppy (v6.5.7)



## mapping


## filtering to remove low-quality reads

## Evaluation of sequencing quality using Nanoplot

## Analysis of pseudouridines using NanoPsu
Details of the analysis pipeline and model implementation are available on the Nanopore_psU GitHub repository: https://github.com/sihaohuanguc/Nanopore_psU (Huang et al.).

## Analysis of rRNA modifications using SignalAlign
Details can be found [in https://github.com/UCSC-nanopore-cgl/signalAlign
](https://github.com/UCSC-nanopore-cgl/signalAlign.)

1. Basecalling
2. Alignment
Combine fastq.gz files
cat *.fastq.gz > bigfile.fastq.gz
Unzip fastq.gz files
Gunzip *. fastq.gz
Run alignment in NanoPSA
nanospa alignment -i ./fastq_pre1_primary.fastq -r 35SrDNA.fasta -o ./alignment
Remove Introns
nanospa remove_intron -i ./alignment/
Extract features
nanospa extract_features -i ./alignment/ -o ./alignment/feature_extraction.csv
Psi prediction
nanospa prediction_psU -i ./alignment/feature_extraction.csv -o ./alignment/Cbf5Gal_psi_prediction.csv
Re-process m6A
nanospa preprocess_m6A -i ./alignment/feature_extraction.csv -o ./alignment/m6A/
m6A prediction
nanospa prediction_m6A -i ./alignment/m6A
<img width="468" height="458" alt="image" src="https://github.com/user-attachments/assets/eba6d06f-a0a5-42d4-93fa-aabb80428b9f" />
