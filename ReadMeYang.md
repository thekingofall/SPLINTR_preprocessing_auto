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



