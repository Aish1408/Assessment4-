# Assessment4-
RNA-Seq Analysis, Growth Data Analysis, and Biological Sequence Diversity in Comparative Genomics

Contributors: Aishwarya, Shrestha, Swanitha
Date: 2024-10-11
Overview

This project involves RNA-Seq data analysis, focusing on gene expression and comparative genomics. The primary goal is to import gene expression data, calculate expression levels, and visualize the distribution of gene expression across multiple samples. Additionally, biological sequence diversity is explored for comparative genomics purposes.
Objectives

Gene Expression Data Import:
Imported gene expression data using read.table() to read TSV files.
The resulting data frame contains genes as rows and samples as columns.
Outcome: A structured data frame for downstream analysis.
Mean Expression Calculation:
Calculated mean expression values across all samples for each gene.
Added a new column, Mean_Expression, to store these values.
Outcome: A data frame with mean expression levels per gene.
Top Gene Identification:
Sorted genes by their mean expression values in descending order.
Extracted the top 10 genes with the highest average expression levels.
Outcome: A list of the top 10 genes with the highest mean expression.
Gene Expression Threshold Analysis:
Determined how many genes have a mean expression value below 10.
Outcome: A count of genes with low expression values, providing insights into the expression landscape.
Histogram of Mean Expression:
Created a histogram to visualize the distribution of mean expression values.
Outcome: A visual representation of the frequency of gene expression levels.
Code Snippets and Explanation

Gene Expression Data Import
r
Copy code
pip = read.table("gene_expression_data.tsv", header=TRUE, row.names=1)
Objective: Import gene expression data.
Outcome: Structured data frame ready for analysis.
Mean Expression Calculation
r
Copy code
pip$Mean_Expression <- rowMeans(pip)
Objective: Calculate the mean expression for each gene.
Outcome: Updated data frame with Mean_Expression.
Top 10 Genes
r
Copy code
ordered_pip <- pip[order(-pip$Mean_Expression), ]
top_10_genes <- head(ordered_pip, 10)
print(top_10_genes)
Objective: Identify top 10 genes with the highest mean expression.
Outcome: List of top 10 genes by expression level.
Gene Expression Threshold
r
Copy code
num_genes_low_mean <- sum(pip$Mean_Expression < 10)
print(num_genes_low_mean)
Objective: Count genes with mean expression below 10.
Outcome: Number of low-expression genes.
Histogram of Mean Expression
r
Copy code
hist(pip$Mean_Expression, main="Histogram of Mean Expression Values", xlab="Mean Expression", col="lightblue", border="black")
Objective: Visualize distribution of gene expression values.
Outcome: Histogram showing the distribution of mean expression.
Required Libraries

R.utils
seqinr
r
Copy code
install.packages("R.utils")
install.packages("seqinr")
library(R.utils)
library(seqinr)
Outcome

This analysis provides insights into the gene expression landscape, highlighting top expressed genes and the distribution of expression levels across the dataset. It offers a foundational approach to RNA-Seq data analysis for growth data and sequence diversity studies in comparative genomics.

Tree Growth Data Analysis README

Introduction

This repository contains R code to analyze tree growth data from two sites (northeast and southwest) over a period of 15 years (2005-2020). The data is read from a CSV file and analyzed for statistical insights on tree circumference and growth rates. Additionally, a box plot is created to visualize the circumference changes, and a t-test is performed to assess the significance of the growth differences between the two sites.
Prerequisites

Before running the analysis, ensure that you have the following packages installed in R:
readr
ggplot2
You can install the necessary packages by running the following commands in R:
r
Copy code
install.packages("readr")
install.packages("ggplot2")
Code Overview

1. Reading the Data
The growth data is read from a remote CSV file using the read_csv() function from the readr package. The URL of the dataset is defined, and the dataset is loaded into a data frame called growth_data.
2. Displaying Initial Data
The first few rows of the dataset are displayed using the head() function, and the column names are printed with colnames(). The dataset contains the following columns:
Site: The region where the tree was measured (northeast or southwest).
TreeID: A unique identifier for each tree.
Circumf_2005_cm, Circumf_2010_cm, Circumf_2015_cm, Circumf_2020_cm: Circumference measurements (in cm) recorded in 2005, 2010, 2015, and 2020.
3. Calculating Mean and Standard Deviation
The mean and standard deviation of tree circumferences for 2005 and 2020 are calculated using the mean() and sd() functions. The results are checked for NA values and stored in a data frame called circumference_stats.
4. Box Plot Visualization
To visualize the tree circumference at the start and end of the study, the data is reshaped into a long format using data.frame(). A box plot is then created using ggplot2, showing the distribution of tree circumferences in 2005 and 2020, with fill colors representing the sites (northeast and southwest).
5. Growth Calculation
The growth of each tree over the 15-year period (from 2005 to 2020) is calculated by subtracting the 2005 circumference from the 2020 circumference. This growth value is stored in a new column Growth. The mean growth for each site is then calculated using the aggregate() function and printed.
6. Statistical Analysis
A Welch Two Sample t-test is conducted to compare the growth rates between the two sites (northeast and southwest). The p-value is used to determine if the growth difference is statistically significant. The confidence interval and sample estimates are also reported.
7. Results
The analysis found that:
Trees in the northeast showed a higher average growth (~48.94 cm) compared to the southwest (~40.73 cm).
The p-value from the t-test was 0.059, suggesting that the difference in growth between the sites is not statistically significant at the 5% significance level.
Files

growth_data.R: The R script file that performs the analysis.
README.md: This file, which provides an overview of the analysis.
How to Run

Open R or RStudio.
Ensure the required packages are installed (see Prerequisites).
Load the R script by sourcing it or running each block of code sequentially.
Review the output statistics and visualizations.
Conclusion

This analysis provides insight into tree growth patterns across two regions over 15 years. The data shows substantial growth in both regions, but no statistically significant difference between them. The box plot and statistical outputs help visualize and quantify these results.
Examining Biological Sequence Diversity in Escherichia coli and Anaerococcus tetradius

Title: Examining Biological Sequence Diversity in E. coli and Anaerococcus tetradius
Authors: Swanitha, Aish, and Shreshta
Date: 2024-10-11
Output: HTML Document
Introduction

This repository contains R scripts that analyze the coding DNA sequences (CDS) of Escherichia coli and Anaerococcus tetradius, exploring their sequence diversity, coding content, and nucleotide compositions. The analysis compares these two microorganisms with distinct ecological niches to highlight evolutionary strategies and genomic adaptations.
The analysis includes:
Counting the total CDS for both organisms.
Summing total coding DNA content.
Plotting the distribution of coding sequence lengths.
Calculating and comparing nucleotide and amino acid frequencies.
Prerequisites

Ensure the following R packages are installed:
readr
seqinr
R.utils
You can install them by running:
r
Copy code
install.packages("readr")
install.packages("seqinr")
install.packages("R.utils")
Code Overview

1. Downloading and Decompressing Coding DNA Sequences
The scripts start by downloading the CDS files for both E. coli and Anaerococcus tetradius from the Ensembl database. The download.file() function is used to retrieve the files, and gunzip() is applied to decompress them.
Inputs:
E. coli CDS URL: <http://ftp.ensemblgenomes.org/pub/bacteria/release-53/fasta/bacteria_0_collection/esc>
Anaerococcus tetradius CDS URL: <https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_12_collection/anaerococcus_tetradius>
Outputs:
Decompressed CDS files ready for further analysis.
2. Reading the CDS Files
The CDS sequences are read using the seqinr::read.fasta() function.
Inputs:
Decompressed CDS files: "ecoli_cds.fa", "anaerococcus_cds.fa"
Outputs:
Variables containing the CDS sequences for each organism (ecoli_cds, anaerococcus_cds).
3. Counting the CDS
The total number of CDS in E. coli and Anaerococcus tetradius is calculated using the length() function.
Outputs:
A data frame summarizing CDS counts:
r
Copy code
Organism                 CDS_Count
Anaerococcus tetradius    1836
E. coli                   4239
4. Calculating Total Coding DNA Content
The lengths of all coding sequences are summed to calculate the total coding DNA content for each organism.
Outputs:
A summary table displaying the total coding DNA content:
r
Copy code
Organism                 Total_Coding_DNA
Anaerococcus tetradius    1,758,063
E. coli                   3,978,528
5. Boxplot of CDS Length Distribution
The distribution of coding sequence lengths for both organisms is visualized using a boxplot. This helps compare the CDS length variations between E. coli and Anaerococcus tetradius.
Outputs:
Boxplots illustrating the CDS length distributions for each organism.
Example summary statistics:
r
Copy code
Mean CDS Length (Anaerococcus): 957 bp
Median CDS Length (Anaerococcus): 804 bp
Mean CDS Length (E. coli): 938 bp
Median CDS Length (E. coli): 831 bp
6. Nucleotide and Amino Acid Composition Analysis
The CDS sequences are analyzed to compute nucleotide frequencies (A, T, G, C) and amino acid frequencies after translating the DNA into protein sequences.
Outputs:
A table displaying the nucleotide composition for both organisms.
A table showing amino acid frequencies for each organism, helping to explore their functional diversity.
Files

sequence_analysis.R: Contains the code for downloading CDS files, analyzing their content, and performing sequence diversity analysis.
README.md: This file, detailing the structure and purpose of each script.


