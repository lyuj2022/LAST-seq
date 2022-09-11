1. mac use ruby install home-brew 

   ```shell
   ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
   ```

2. use homebrew install wget

   ```shell
   brew install wget
   ```

3. use wget install sratoolkit

   ```shell
   wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6-1/sratoolkit.2.9.6-1-mac64.tar.gz
   ```

4. decompress sratoolkit

5. ```
   tar zvxf sratoolkit.2.9.6-1-mac64.tar.gz
   ```

5. Test installation

7. ```shell
   ../bin/fastq-dump
   ```

8. download batch SRA data by prefetch, SraAccList.txt can be export from SRA website

9. ```shell
   ./prefetch -O ./ --option-file /Users/lyuj2/Downloads/SraAccList.txt
   ```

10. extract fastq file from sra file

11. ```shell
    fastq-dump --gzip --split-files SRR6232298.sra
    
    #--split-files is to split one SRA into two fastq, --gzip generates fastq.gz.
    ```

12. for batch work, you should write a script

13. ```shell
    !/bin/sh
    module load sratoolkit
    for i in *sra
    do
    echo $i
    fastq-dump -I --gzip --split-files $i
    done
    
    # -I will use the suffix _1.fastq and _2.fastq
    
    #if work on biowulf, please consult on https://hpc.nih.gov/apps/sratoolkit.html
    #but note that, sratoolkit is a single thread app, it will cost 10min to encode one SRR file, so I think it's better to use batch array (i.e. Swarm) instead of circular script.
    ```

    


