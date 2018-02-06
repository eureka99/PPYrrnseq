# RNAseq
differentially expression genes of plants (tissues and growth stage)

[zhangxuan@c0102 test_Trinity_Assembly]$ ./runMe.sh 
module () {  eval `/usr/bin/modulecmd bash $*`
}
#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

../../Trinity --seqType fq --max_memory 2G \
              --left reads.left.fq.gz \
              --right reads.right.fq.gz \
              --SS_lib_type RF \
              --CPU 4 


     ______  ____   ____  ____   ____  ______  __ __
    |      ||    \ |    ||    \ |    ||      ||  |  |
    |      ||  D  ) |  | |  _  | |  | |      ||  |  |
    |_|  |_||    /  |  | |  |  | |  | |_|  |_||  ~  |
      |  |  |    \  |  | |  |  | |  |   |  |  |___, |
      |  |  |  .  \ |  | |  |  | |  |   |  |  |     |
      |__|  |__|\_||____||__|__||____|  |__|  |____/


Left read files: $VAR1 = [
          'reads.left.fq.gz'
        ];
Right read files: $VAR1 = [
          'reads.right.fq.gz'
        ];
Trinity version: Trinity-v2.5.1
-currently using the latest production release of Trinity.

Thursday, January 4, 2018: 20:19:19	CMD: java -Xmx64m -XX:ParallelGCThreads=2  -jar /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/support_scripts/ExitTester.jar 0
Thursday, January 4, 2018: 20:19:20	CMD: java -Xmx64m -XX:ParallelGCThreads=2  -jar /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/support_scripts/ExitTester.jar 1
Thursday, January 4, 2018: 20:19:20	CMD: mkdir -p /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir
Thursday, January 4, 2018: 20:19:20	CMD: mkdir -p /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis


----------------------------------------------------------------------------------
-------------- Trinity Phase 1: Clustering of RNA-Seq Reads  ---------------------
----------------------------------------------------------------------------------

---------------------------------------------------------------
------------ In silico Read Normalization ---------------------
-- (Removing Excess Reads Beyond 50 Coverage --
---------------------------------------------------------------

# running normalization on reads: $VAR1 = [
          [
            '/lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/reads.left.fq.gz'
          ],
          [
            '/lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/reads.right.fq.gz'
          ]
        ];


Thursday, January 4, 2018: 20:19:20	CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/insilico_read_normalization.pl --seqType fq --JM 2G  --max_cov 50 --CPU 4 --output /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization   --max_pct_stdev 10000  --SS_lib_type RF  --left /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/reads.left.fq.gz --right /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/reads.right.fq.gz --pairs_together --PARALLEL_STATS  
Converting input files. (both directions in parallel)CMD: seqtk-trinity seq -r -A <(gunzip -c /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/reads.left.fq.gz) >> left.fa
CMD: seqtk-trinity seq -A <(gunzip -c /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/reads.right.fq.gz) >> right.fa
CMD finished (3 seconds)
CMD finished (3 seconds)
CMD: touch left.fa.ok
CMD finished (0 seconds)
CMD: touch right.fa.ok
CMD finished (0 seconds)
Done converting input files.CMD: cat left.fa right.fa > both.fa
CMD finished (0 seconds)
CMD: touch both.fa.ok
CMD finished (0 seconds)
-------------------------------------------
----------- Jellyfish  --------------------
-- (building a k-mer catalog from reads) --
-------------------------------------------

CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/..//trinity-plugins/jellyfish/bin/jellyfish count -t 4 -m 25 -s 100000000  both.fa
CMD finished (2 seconds)
CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/..//trinity-plugins/jellyfish/bin/jellyfish histo -t 4 -o jellyfish.K25.min2.kmers.fa.histo mer_counts.jf
CMD finished (1 seconds)
CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/..//trinity-plugins/jellyfish/bin/jellyfish dump -L 2 mer_counts.jf > jellyfish.K25.min2.kmers.fa
CMD finished (0 seconds)
CMD: touch jellyfish.K25.min2.kmers.fa.success
CMD finished (0 seconds)
CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/..//Inchworm/bin/fastaToKmerCoverageStats --reads left.fa --kmers jellyfish.K25.min2.kmers.fa --kmer_size 25  --num_threads 4  > left.fa.K25.stats
CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/..//Inchworm/bin/fastaToKmerCoverageStats --reads right.fa --kmers jellyfish.K25.min2.kmers.fa --kmer_size 25  --num_threads 4  > right.fa.K25.stats
-reading Kmer occurrences...
-reading Kmer occurrences...

 done parsing 100973 Kmers, 100873 added, taking 0 seconds.

 done parsing 100973 Kmers, 100873 added, taking 0 seconds.
STATS_GENERATION_TIME: 0 seconds.
CMD finished (0 seconds)
STATS_GENERATION_TIME: 0 seconds.
CMD finished (1 seconds)
CMD: touch left.fa.K25.stats.ok
CMD finished (0 seconds)
CMD: touch right.fa.K25.stats.ok
CMD finished (0 seconds)
-sorting each stats file by read name.
CMD: /bin/sort -k5,5 -T . -S 1G left.fa.K25.stats > left.fa.K25.stats.sort
CMD: /bin/sort -k5,5 -T . -S 1G right.fa.K25.stats > right.fa.K25.stats.sort
CMD finished (0 seconds)
CMD finished (0 seconds)
CMD: touch left.fa.K25.stats.sort.ok
CMD finished (0 seconds)
CMD: touch right.fa.K25.stats.sort.ok
CMD finished (0 seconds)
CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/..//util/support_scripts//nbkc_merge_left_right_stats.pl --left left.fa.K25.stats.sort --right right.fa.K25.stats.sort --sorted > pairs.K25.stats
-opening left.fa.K25.stats.sort
-opening right.fa.K25.stats.sort
-done opening files.
CMD finished (0 seconds)
CMD: touch pairs.K25.stats.ok
CMD finished (0 seconds)
CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/..//util/support_scripts//nbkc_normalize.pl pairs.K25.stats 50 10000 > pairs.K25.stats.C50.pctSD10000.accs
23031 / 30575 = 75.33% reads selected during normalization.
0 / 30575 = 0.00% reads discarded as likely aberrant based on coverage profiles.
0 / 30575 = 0.00% reads missing kmer coverage (N chars included?).
CMD finished (0 seconds)
CMD: touch pairs.K25.stats.C50.pctSD10000.accs.ok
CMD finished (0 seconds)
CMD: touch /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/reads.left.fq.gz.normalized_K25_C50_pctSD10000.fq.ok
CMD finished (0 seconds)
CMD: touch /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/reads.right.fq.gz.normalized_K25_C50_pctSD10000.fq.ok
CMD finished (0 seconds)
CMD: ln -sf /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/reads.left.fq.gz.normalized_K25_C50_pctSD10000.fq left.norm.fq
CMD finished (0 seconds)
CMD: ln -sf /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/reads.right.fq.gz.normalized_K25_C50_pctSD10000.fq right.norm.fq
CMD finished (0 seconds)
-removing tmp dir /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/tmp_normalized_reads


Normalization complete. See outputs: 
	/lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/reads.left.fq.gz.normalized_K25_C50_pctSD10000.fq
	/lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/reads.right.fq.gz.normalized_K25_C50_pctSD10000.fq
Thursday, January 4, 2018: 20:19:30	CMD: touch /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/normalization.ok
Converting input files. (in parallel)Thursday, January 4, 2018: 20:19:30	CMD: cat /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/left.norm.fq | seqtk-trinity seq -r -A - >> left.fa
Thursday, January 4, 2018: 20:19:30	CMD: cat /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/insilico_read_normalization/right.norm.fq | seqtk-trinity seq -A - >> right.fa
Thursday, January 4, 2018: 20:19:30	CMD: touch left.fa.ok
Thursday, January 4, 2018: 20:19:30	CMD: touch right.fa.ok
Thursday, January 4, 2018: 20:19:30	CMD: touch left.fa.ok right.fa.ok
Thursday, January 4, 2018: 20:19:30	CMD: cat left.fa right.fa > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa
Thursday, January 4, 2018: 20:19:30	CMD: touch /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa.ok
-------------------------------------------
----------- Jellyfish  --------------------
-- (building a k-mer catalog from reads) --
-------------------------------------------

* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/trinity-plugins/jellyfish/bin/jellyfish count -t 4 -m 25 -s 100000000  /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/trinity-plugins/jellyfish/bin/jellyfish dump -L 1 mer_counts.jf > jellyfish.kmers.fa
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/trinity-plugins/jellyfish/bin/jellyfish histo -t 4 -o jellyfish.kmers.fa.histo mer_counts.jf
----------------------------------------------
--------------- Inchworm ---------------------
-- (Linear contig construction from k-mers) --
----------------------------------------------



* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/Inchworm/bin//inchworm --kmers jellyfish.kmers.fa --run_inchworm -K 25 -L 25 --monitor 1   --num_threads 4  --PARALLEL_IWORM  > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa.tmp
* Running CMD: mv /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa.tmp /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa
Thursday, January 4, 2018: 20:19:33	CMD: touch /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa.finished
--------------------------------------------------------
-------------------- Chrysalis -------------------------
-- (Contig Clustering & de Bruijn Graph Construction) --
--------------------------------------------------------

inchworm_target: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa
bowite_reads_fa: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa
chrysalis_reads_fa: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/support_scripts/filter_iworm_by_min_length_or_cov.pl /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa 100 10 > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/inchworm.K25.L25.fa.min100
* Running CMD: bowtie2-build --threads 4 -o 3 /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/inchworm.K25.L25.fa.min100 /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/inchworm.K25.L25.fa.min100 1>/dev/null
* Running CMD: bash -c " set -o pipefail;bowtie2 --local -k 2 --threads 4 -f --score-min G,46,0 -x /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/inchworm.K25.L25.fa.min100 /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa  | samtools view -@ 4 -F4 -Sb - | samtools sort -m 268435456 -@ 4 -no - - > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/iworm.bowtie.nameSorted.bam" 
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/support_scripts/scaffold_iworm_contigs.pl /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/iworm.bowtie.nameSorted.bam /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/iworm_scaffolds.txt
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/Chrysalis/GraphFromFasta -i /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa -r /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa -min_contig_length 200 -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05 -t 4 -k 24 -kk 48  -strand  -scaffolding /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/iworm_scaffolds.txt  > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/iworm_cluster_welds_graph.txt
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/Chrysalis/BubbleUpClustering -i /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/inchworm.K25.L25.fa  -weld_graph /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/iworm_cluster_welds_graph.txt -min_contig_length 200  > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/GraphFromIwormFasta.out
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/Chrysalis/CreateIwormFastaBundle -i /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/GraphFromIwormFasta.out -o /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/bundled_iworm_contigs.fasta -min 200
* Running CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/Chrysalis/ReadsToTranscripts -i /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/both.fa -f /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/bundled_iworm_contigs.fasta -o /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/readsToComponents.out -t 4 -max_mem_reads 50000000  -strand 
* Running CMD: /bin/sort -T . -S 2G -k 1,1n /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/readsToComponents.out > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/chrysalis/readsToComponents.out.sort
Thursday, January 4, 2018: 20:19:43	CMD: mkdir -p /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/read_partitions/Fb_0/CBin_0
Thursday, January 4, 2018: 20:19:44	CMD: touch /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/partitioned_reads.files.list.ok
Thursday, January 4, 2018: 20:19:44	CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/support_scripts/write_partitioned_trinity_cmds.pl --reads_list_file /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/partitioned_reads.files.list --CPU 1 --max_memory 1G  --run_as_paired  --SS_lib_type F  --seqType fa --trinity_complete --full_cleanup  > recursive_trinity.cmds
Thursday, January 4, 2018: 20:19:44	CMD: touch recursive_trinity.cmds.ok
Thursday, January 4, 2018: 20:19:44	CMD: touch recursive_trinity.cmds.ok


--------------------------------------------------------------------------------
------------ Trinity Phase 2: Assembling Clusters of Reads ---------------------
--------------------------------------------------------------------------------

Thursday, January 4, 2018: 20:19:44	CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/trinity-plugins/BIN/ParaFly -c recursive_trinity.cmds -CPU 4 -v -shuffle 
Number of Commands: 32
succeeded(32)   100% completed.       

All commands completed successfully. :-)



** Harvesting all assembled transcripts into a single multi-fasta file...

Thursday, January 4, 2018: 20:20:46	CMD: find /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/read_partitions/ -name '*inity.fasta'  | /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/support_scripts/partitioned_trinity_aggregator.pl TRINITY_DN > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta.tmp
-relocating /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta.tmp to /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta
Thursday, January 4, 2018: 20:20:46	CMD: mv /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta.tmp /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta


###################################################################
Trinity assemblies are written to /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta
###################################################################


Thursday, January 4, 2018: 20:20:46	CMD: /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/util/support_scripts/get_Trinity_gene_to_trans_map.pl /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta > /lustre/home/zhangxuan/gxy/software/trinityrnaseq-Trinity-v2.5.1/sample_data/test_Trinity_Assembly/trinity_out_dir/Trinity.fasta.gene_trans_map

##### Done Running Trinity #####

if [ $* ]; then
    # check full-length reconstruction stats:
    ./test_FL.sh --query trinity_out_dir/Trinity.fasta --target __indiv_ex_sample_derived/refSeqs.fa --no_reuse
Fi

