# SPLINTR Barcode Library Pre-processing. 

For more information see the paper [Fennell and Vassiliadis et al., 2021](https://www.nature.com/articles/s41586-021-04206-7)

## Software Dependencies
- [Python 3.7+](https://www.python.org/downloads/)
  - [Biopython](https://biopython.org/wiki/Download)
- [FLASh](https://ccb.jhu.edu/software/FLASH/)
- [Starcode](https://github.com/gui11aume/starcode) 
- [R 4.0+](https://cran.r-project.org/)
  - [Tidyverse](https://www.tidyverse.org/)
  
## Procedure to generate SPLINTR reference barcode libraries:
This process assumes familiarity with and access to a HPC cluster running the SLURM scheduler. 

1. Amplify SPLINTR barcodes from each library plasmid pool by PCR (see attached PCR protocol).
  - Perform PCRs in technical duplicate
  - Use a different index combination in PCR2 for each replicate

2. Sequence library pool to a depth of at least 100-200M reads (~100X coverage for a 1 million barcode library) per PCR reaction.
  - Lower diversity libraries will require less sequencing. 
  - Ideally perform single end 150bp or paired end 75bp sequencing, DO NOT perform single end 75bp as reads will not cover the entire barcode.	

3. Generate a text file containing full path of each library fastq file. See example `fastqs.txt`

4. Confirm sequencing quality by running FastQC
  -  **[OPTIONAL]** if paired end 75bp sequencing was performed. Run `flash.sh` to combine paried-end reads into a single overlapped read that covers the entire barcode region.

5. Extract barcode reads from fastq files using `filter_barcodes.sbatch` which submits `extractBarcodeReads.py` scripts to a SLURM HPC

6. Run Starcode on the extracted barcode reads to cluster barcodes according to Edit distance. See `starcode.sbatch`

7. Perform QC and analysis of Starcode output in R using `library-analysis-template.Rmd`

Output will include a fasta file of barcode sequences within each library which can be used to generate a reference library for subsequent barcode alignment.
See example output in `library-analysis-template.html` and `test`
