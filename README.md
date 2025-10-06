# Dip-Patel
hello world 
This repo is a demonstration of using Git with RStudio.
 
SLE777 – Applied Bioinformatics (Assessment 4)

Overview

This repository contains all R scripts, data files, and outputs created for Assessment 4 of the SLE777 – Applied Bioinformatics unit. The purpose of this project is to apply R programming skills to analyse and visualise biological data, and to explore sequence diversity between bacterial genomes. The assessment is divided into two major parts: Part 1, which focuses on data wrangling and visualisation of experimental datasets, and Part 2, which performs comparative genomic and proteomic analyses. The entire workflow is implemented in RStudio using an R Markdown report (Assessment4.Rmd), which sources all scripts and compiles outputs into a single reproducible document. This README provides a detailed description of the purpose, inputs, and outputs of every script used in the analysis, as required by the assessment instructions and marking rubric.

# Part 1 – Data Wrangling and Visualisation

Part 1 develops skills in importing, cleaning, and analysing biological data using R. It involves two datasets: gene_expression.tsv, containing RNA-seq count data, and growth_data.csv, containing tree growth measurements over time. Each script in this part is designed to solve specific questions (Q1–Q10) from the assignment brief and to generate visual and numerical summaries suitable for inclusion in the report.

Script: part1_gene_expression.R

Purpose:
This script focuses on exploring RNA-seq gene expression data to calculate descriptive statistics and visualise expression patterns across genes. It answers Questions 1 to 5 from Part 1 of the assessment by performing mathematical and visual operations on the dataset.

Inputs:
The input file for this script is data/gene_expression.tsv. This tab-separated text file contains RNA-seq read count data for three different samples, with each row representing a gene and each column showing the read counts for that gene across samples.

Outputs:
The script generates several outputs stored in the outputs/part1/ folder. It produces a CSV file containing the first six rows of the dataset (q1_head_gene_expression.csv), a modified dataset that includes a new column for the mean expression value (q2_mean_added.csv), a table listing the ten genes with the highest mean expression values (q3_top10_genes.csv), and a text file recording the number of genes with mean expression less than 10 (q4_count_low_mean.txt). In addition, it creates a histogram figure that visualises the distribution of mean expression values across all genes (q5_hist_mean_expression.png). Together, these outputs provide both numerical and graphical insights into gene expression variation.

Script: part1_growth_data.R

Purpose:
This script analyses long-term tree growth data to compare growth patterns between control and treatment sites. It answers Questions 6 to 10 from the assessment and applies descriptive statistics, visualisation, and hypothesis testing.

Inputs:
The input file for this script is data/growth_data.csv. It contains tree circumference measurements recorded at the start and end of a 20-year period for both control and treatment sites. Each record includes the site type, timepoint, and corresponding tree circumference data.

Outputs:
This script produces multiple output files stored in outputs/part1/. It first saves a list of column names to q6_column_names.txt to confirm the dataset structure. Next, it calculates the mean and standard deviation of tree circumference at both the start and end of the study for each site, saving the results in q7_summary_start_end.csv. It then generates a boxplot (q8_boxplot_circumference.png) that visually compares the distributions of tree circumference between the control and treatment sites at different timepoints. The script also computes mean growth over the last ten years at each site (q9_mean_growth_10y.csv) and performs a two-sample t-test to test for significant differences in ten-year growth between the two sites, saving the statistical results to q10_ttest_result.txt. These outputs collectively summarise growth patterns, variability, and site-specific differences.


# Part 2 – Examining Sequence Diversity

Part 2 examines sequence diversity and comparative genomics between Escherichia coli and the assigned organism, Alphaproteobacteria bacterium (GCA_002697805). It involves downloading genome data, calculating DNA and protein composition, and assessing codon usage and k-mer patterns. Each script in this part contributes to building a complete comparative genomic workflow that integrates bioinformatics tools with R-based analysis and visualisation.

Script: part2_download_sequences.R

Purpose:
This script automates the retrieval of coding DNA sequences (CDS) and protein FASTA files for both E. coli and the assigned organism Alphaproteobacteria bacterium. It ensures that the sequence data are downloaded, stored in the correct format, and available for all downstream analyses.

Inputs:
The script does not require a user-supplied dataset but instead uses genome accession numbers and NCBI or ENA links to fetch sequence data directly from public databases. It requires an active internet connection for data retrieval.

Outputs:
The downloaded files are saved in the data/ directory as ecoli_cds.fna.gz, ecoli_proteins.faa.gz, target_cds.fna.gz, and target_proteins.faa.gz. In addition, the script generates a download log file (download_log.txt) in the outputs/part2/ folder, recording the date, time, and completion status of each download.

Script: part2_sequence_stats.R

Purpose:
This script performs basic genome-level statistical analysis for both species. It calculates the number of coding sequences, total coding DNA length, and gene-length distributions to provide an overview of genome architecture.

Inputs:
The script requires the CDS FASTA files of both organisms (ecoli_cds.fna.gz and target_cds.fna.gz) and their corresponding protein FASTA files (ecoli_proteins.faa.gz and target_proteins.faa.gz).

Outputs:
The results are stored in outputs/part2/. The script produces tables summarising the number of coding sequences per organism (q1_cds_counts.csv), the total length of all coding sequences in base pairs (q2_total_coding_bp.csv), and descriptive statistics such as mean and median gene lengths (q3_cds_length_summary.csv). It also generates a boxplot (q3_cds_length_boxplot.png) that visually compares the distribution of coding sequence lengths between E. coli and Alphaproteobacteria bacterium.

Script: part2_nuc_aa_frequency.R

Purpose:
This script calculates and compares the nucleotide and amino-acid composition of the genomes and proteomes of both organisms, identifying compositional biases and potential evolutionary differences.

Inputs:
The required inputs are the CDS and protein FASTA files of both species that were downloaded by the previous script.

Outputs:
The script outputs CSV files containing nucleotide frequency results (q4_nucleotide_frequency.csv) and amino-acid frequency data (q4_amino_acid_frequency.csv). It also generates two barplots – one for nucleotide composition (q4_nucleotide_frequency_barplot.png) and one for amino-acid composition (q4_amino_acid_frequency_barplot.png) – that visually compare the base and amino-acid proportions between the two organisms.

Script: part2_codon_usage.R

Purpose:
This script analyses codon usage bias by counting codon occurrences and calculating Relative Synonymous Codon Usage (RSCU) values, which quantify preferences for specific codons over synonymous alternatives.

Inputs:
The inputs are the CDS FASTA files of both organisms (ecoli_cds.fna.gz and target_cds.fna.gz).

Outputs:
Outputs are stored in outputs/part2/. They include raw codon count data (q5_codon_counts.csv), RSCU values (q5_rscu_values.csv), and a summary of total codon usage across all genes (q5_total_codon_summary.csv). The script also produces a visual representation of codon usage bias through a barplot highlighting the top 20 most frequently used codons (q5_rscu_barplot_top20.png).

Script: part2_kmer_enrichment.R

Purpose:
This script explores short peptide motifs (k-mers) of lengths 3 to 5 in the protein sequences of both organisms. It identifies motifs that are statistically over-represented or under-represented and compares these patterns between E. coli and the target bacterium.

Inputs:
The required inputs are the protein FASTA files (ecoli_proteins.faa.gz and target_proteins.faa.gz) downloaded in the first Part 2 script.

Outputs:
The script produces separate CSV files for each k-mer size (q6_kmer_enrichment_k3.csv, q6_kmer_enrichment_k4.csv, q6_kmer_enrichment_k5.csv) summarising enrichment scores for all possible motifs. It also generates a comparative plot (q6_kmer_top10_over_under.png) that visually highlights the ten most over-represented and under-represented motifs in the allocated organism and compares them to E. coli.

Reproducibility and Organisation

All scripts are written for R version 4.3 or later and are compatible with standard RStudio environments. They use the here::here() function to ensure that all file paths work consistently across systems. Outputs are automatically saved into organised folders under outputs/part1/ and outputs/part2/. The R Markdown report (Assessment4.Rmd) sources these scripts sequentially, embedding their outputs and plots directly into the final report. Each script includes meaningful comments explaining every step to satisfy the rubric’s requirement for well-documented code.

Allocated Organism

According to the assignment instructions, the allocated organism for Student ID S225592308 is Alphaproteobacteria bacterium (GCA_002697805). This species is compared against the reference organism Escherichia coli in all analyses conducted in Part 2.

Conclusion

This repository provides a complete, well-structured R workflow for both data analysis and sequence diversity studies. Each script has a clearly defined purpose, with detailed explanations of inputs and outputs to ensure reproducibility and transparency. The outputs include all required plots, tables, and summary files, and the code is fully documented and functional. Together, these components meet the expectations of the assessment manual and rubric for a comprehensive and professional submission.