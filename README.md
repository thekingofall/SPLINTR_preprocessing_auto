
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


### 封装步骤


```
/home/YangZongmian/SPLINTR/SPLINTR_preprocessing/run_splintr_pipeline.sh -i /home/YangZongmian/SPLINTR/X101SC24127971-Z01-J027/01.RawData -o /home/YangZongmian/SPLINTR/AnalysisResults
```

```
(SPLINTR) [YangZongmian@master01 01.RawData]$ echo "alias run_splintr='bash /home/YangZongmian/SPLINTR/SPLINTR_preprocessing/run_splintr_pipeline.sh'" >> ~/.bashrc
(SPLINTR) [YangZongmian@master01 01.RawData]$ source ~/.bashrc
(base) [YangZongmian@master01 01.RawData]$ run_splintr -h
Usage: /home/YangZongmian/SPLINTR/SPLINTR_preprocessing/run_splintr_pipeline.sh -i <input_directory> [-o <output_directory>] [-t <threads>] [-h]

Processes SPLINTR R1 FASTQ files found in subdirectories of the input directory.

Options:
  -i <input_directory>   : Directory containing sample subdirectories (e.g., 1/, 2/, 10/)
                           Each subdirectory must contain a corresponding _1.fq.gz file.
                           (Mandatory)
  -o <output_directory>  : Directory where processed files will be saved.
                           (Optional, default: ./Processed_SPLINTR_Output)
  -t <threads>           : Number of threads for starcode.
                           (Optional, default: 8)
  -h                     : Display this help message.

Example:
  /home/YangZongmian/SPLINTR/SPLINTR_preprocessing/run_splintr_pipeline.sh -i /path/to/01.RawData -o /path/to/results
  /home/YangZongmian/SPLINTR/SPLINTR_preprocessing/run_splintr_pipeline.sh -i /path/to/01.RawData  # Uses default output directory
```

## 以上步骤均可以忽略
> 下面已经封装好所有的步骤，一行代码即可


```
(base) [YangZongmian@master01 ~]$ python --version
Python 2.7.5
(base) [YangZongmian@master01 ~]$ conda deactivate
[YangZongmian@master01 ~]$ 



```

```
(base) [YangZongmian@master01 X101SC24127971-Z01-J027]$ run_splintr -i 输入文件夹名字 -o 输出文件夹名字
```

```
cd /home/YangZongmian/SPLINTR/X101SC24127971-Z01-J027 ### 进入文件夹

(base) [YangZongmian@master01 X101SC24127971-Z01-J027]$ run_splintr -i 01.RawData/ -o Processed_SPLINTR_Output_Parallel
```


## 或者绝对路径

```
run_splintr -i /home/YangZongmian/SPLINTR/X101SC24127971-Z01-J027/01.RawData/ -o  /home/YangZongmian/SPLINTR/X101SC24127971-Z01-J027/Processed_SPLINTR_Output_Parallel2
```




