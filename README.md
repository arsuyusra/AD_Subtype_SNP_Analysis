# Alzheimer’s Disease Subtypes SNP Analysis

## A Comprehensive Study of Genetic Variations in Alzheimer's Disease Subtypes

This repository contains the code for analyzing and visualizing genetic variations in Alzheimer's disease (AD) subtypes using whole-genome sequencing (WGS) data. 
The scripts provided enable the identification of single-nucleotide polymorphisms (SNPs) associated with different AD subtypes (atypical, intermediate, and typical)
and healthy control groups. Using Python, bioinformatics tools like PLINK, and visualization libraries such as Seaborn, this code conducts statistical analysis 
(including Fisher’s exact test) and generates heatmaps to visualize allele frequencies and significance of SNPs across AD subtypes. The goal is to offer a reproducible
analysis pipeline for genetic research in AD.

## Setup and Developer Instructions
1.  Clone the repository:

   ```bash
   git clone https://github.com/yourusername/AD_Subtype_SNP_Analysis.git
   cd AD_Subtype_SNP_Analysis

2. Create a virtual environment:

python3 -m venv venv
source venv/bin/activate  # On Windows, use `venv\Scripts\activate`

3. Install the required tendencies:

pip3 install numpy pandas seaborn matplotlib scipy  

To contribute to this project, ensure that you have Python 3.9 or higher installed on your computer, along with the required dependencies
listed in requirements.txt. If you wish to add or modify analysis scripts, make sure that your changes maintain the code's structure,
particularly in data processing and visualization. Follow the proper coding standards, write clear comments, and document any new
functionality added.









