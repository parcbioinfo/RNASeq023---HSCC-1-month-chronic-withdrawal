RNASeq023 - HSCC after chronic 1-month EtOH

The analysis related files are distributed amond 3 directories:
  data: input data and intermediate files
    data/Alignment_mm10 --> bam and bai files
    data/Normalization --> all input files needed to generate normalized read counts
  scripts: all scripts used for analysis
  analysis: final analysis results files 

Alignment:
  Since Sample file numbers had some breaks in sequence, I used 3 scripts to run the alignment simultaneously for different batches of samples.
    In lawrence:
      nohup ./run_star_alignment_1.sh > ./run_star_alignment_1.log 2>&1 </dev/null &
      nohup ./run_star_alignment_1b.sh > ./run_star_alignment_1b.log 2>&1 </dev/null &
      nohup ./sort_and_index_1b.sh > ./sort_and_index_1b.log 2>&1 </dev/null &

    In exacloud:
      condor-submit align.sub
      align.sub:
        executable              = run_star_alignment_2.sh
        log                     = submit_run_star_alignment.log
        output                  = submit_run_star_alignment.out
        error                   = submit_run_star_alignment.err
        request_cpus            = 20
        request_memory          = 64 GB
        notification            =Complete
        notify_user             =darakjia@ohsu.edu
        queue
    
Counts:
  nohup ../../scripts/run_bedtools_coverage.sh > ./run_bedtools_coverage.log 2>&1 </dev/null &
  mv RNASeq023_mm10_coverage_splitoption.txt ../../analysis/
  ls -d Sample_RNA150217RH_* > samples1.txt
  cd ../150423_D00735_0035_BC7D48ACXX/
  ls -d Sample_RNA150217RH_* >> ../150423_D00735_0034_AC791EACXX/samples1.txt
  cd ../150423_D00735_0034_AC791EACXX/
  #### After removing unnecessary substrings from sample names and tab separating core id from lab id:
  mv samples1.txt sample_key.txt
  mv sample_key.txt ../data/
  cd ../data/
  sort -n sample_key.txt > sample_key_sorted
  mv sample_key_sorted sample_key_sorted.txt
  cd Alignment_mm10/
  ls Alignment_mm10/*sorted*bam > bam_files_counts_header_order.txt
  #### Removed unnecessary substrings from sample names and tab separating core id from lab id
  
Normalization:
  used --> scripts/selectNormalizeGeneExonCounts_RNASeq023.R
  
  