import os
import numpy as np
import pandas as pd

path_fastp_files = "/shared/projects/origamis/wgs/sphegodes/HiSeq/"
path_github = "/shared/projects/origamis/wgs/"
path_output = "/shared/projects/origamis/wgs/"
path_ref_genome = "/shared/projects/origamis/ref_genome/sphegodes/GCA_040285675.1_Osph-v2.3_genomic.fna"

raw_samples = os.popen("ls /shared/projects/origamis/wgs/sphegodes/HiSeq/ | grep '_1.fastq.gz'").read().split("\n")[:-1]
SAMPLES = [x.split("_1.fastq.gz")[0] for x in raw_samples]
print(SAMPLES)
REGIONS = pd.read_csv("/shared/projects/origamis/ref_genome/sphegodes/chr_regions.txt",header=None)[0].tolist()
#REGIONS = REGIONS[0]

INDIV = ["sphegodes_hiseq"]
chrom = os.popen("grep '>' /shared/projects/origamis/ref_genome/sphegodes/GCA_040285675.1_Osph-v2.3_genomic.fna | grep 'CM' | awk '{print $1}'").read().split("\n")[:-1]
CHROM = [x.split(">")[1] for x in chrom]
CHROM = CHROM[:-2]
print(CHROM)

rule all:
        input:
                #fastp_R1 = expand(path_output+"fastp/{sample}_R1_fastp.fastq.gz",sample=SAMPLES),
                #fastp_R2 = expand(path_output+"fastp/{sample}_R2_fastp.fastq.gz",sample=SAMPLES),
                #report_html = expand(path_output+"fastp/fastp_report_{sample}.html",sample=SAMPLES),
                #report_json = expand(path_output+"fastp/fastp_report_{sample}.json",sample=SAMPLES)
                #raw_cram = expand(path_output+"mapping/{sample}.cram",sample=SAMPLES)
                #raw_bam = expand(path_output+"mapping/{sample}.bam",sample=SAMPLES)
                #markdup_cram = expand(path_output+"mapping/{sample}_markdup.cram",sample=SAMPLES)
                #markdup_bam = expand(path_output+"mapping/{sample}_markdup.bam",sample=SAMPLES)
                #samtools_stats = expand(path_output+"mapping/{sample}_samtools_stats.stats",sample=SAMPLES),
                #flagstat = expand(path_output+"mapping/{sample}_flagstat.txt",sample=SAMPLES)
                #mosepth_output = expand(path_output + "mosdepth/{sample}/dist.html",sample=SAMPLES), 
                #ngsParalog = path_output + "ngsParalog/ERR5100914.lr",
                #mito_plasto = expand(path_output + "getorganelle/{sample}_plastome_mitogenome/extended_K115.assembly_graph.fastg",sample=SAMPLES),
                #vcf_regions =  expand(path_output + "vcf/{regions}.vcf.gz",regions=REGIONS),
                #vcf_filtered_only_variant = path_output + "ophrys_variant_only.vcf.gz"
                #vcf_mark_filter = path_output + "ophrys_mark_filter.vcf.gz",
                #vcf_filtered = path_output + "ophrys_filtered.vcf.gz"
                bpp = expand(path_github+"output/iSMC/{indiv}/{chrom}/{chrom}.bpp",indiv=INDIV,chrom=CHROM),
                tab = expand(path_github+"output/iSMC/{indiv}/{chrom}/{chrom}.tab",indiv=INDIV,chrom=CHROM),
                decoding_label = expand(path_github+"output/iSMC/{indiv}/{chrom}/my_dataset_diploid_decoding_labels.txt",indiv=INDIV,chrom=CHROM)

rule fastp:
        input:
                raw_R1 = path_fastp_files + "{sample}_1.fastq.gz",
                raw_R2 = path_fastp_files + "{sample}_2.fastq.gz"
        output:
                fastp_R1 = path_output + "fastp/{sample}_R1_fastp.fastq.gz",
                fastp_R2 = path_output + "fastp/{sample}_R2_fastp.fastq.gz",
                report_html = path_output + "fastp/fastp_report_{sample}.html",
                report_json = path_output + "fastp/fastp_report_{sample}.json"
        message:
                "Fastp processing : {wildcards.sample}"
        log:
                stdout = path_github + "log/{sample}_fastp.log"
        benchmark:
                path_github + "benchmark/{sample}_fastp_benchmark.txt"
        params:
                thread = 1
        shell :
               "(fastp "
               "-i {input.raw_R1} "
               "-I {input.raw_R2} "
               "-o {output.fastp_R1} "
               "-O {output.fastp_R2} "
               "--trim_poly_g "
               "--correction "
               "--low_complexity_filter "
               "--html {output.report_html} "
               "--json {output.report_json} "
               "--report_title {wildcards.sample} "
               "--thread {params.thread} "
               "--dont_overwrite) 2> {log.stdout}"

rule reference_mapping:
        input:
                fastp_R1 = path_output + "fastp/{sample}_R1_fastp.fastq.gz",
                fastp_R2 = path_output + "fastp/{sample}_R2_fastp.fastq.gz",
                reference_genome = path_ref_genome
        output:
                raw_cram = path_output + "mapping/{sample}.cram"
                #raw_bam = path_output + "mapping/{sample}.bam"
        message:
                "Reference mapping: {wildcards.sample}"
        log:
                stdout = path_github+"log/{sample}_reference_mapping.log"
        benchmark:
                path_github+"benchmark/{sample}_reference_mapping_benchmark.txt"
        params:
                threads = 16
        shell:
                "(bwa mem "
                "-M "
                "-t {params.threads} "
                "{input.reference_genome} "
                "{input.fastp_R1} "
                "{input.fastp_R2} | "
                #"samtools view -b > {output.raw_bam}) 2> {log.stdout}"
                "samtools view -C "
                "-T {input.reference_genome} "
                "-o {output.raw_cram} "
                "-) 2> {log.stdout}"

rule sort_readgroup_markdup_bam:
        input:
                raw_cram = path_output + "mapping/{sample}.cram",
                #raw_bam = path_output + "mapping/{sample}.bam",
                reference_genome = path_ref_genome
        output:
                markdup_cram = path_output + "mapping/{sample}_markdup.cram",
                #markdup_bam = path_output + "mapping/{sample}_markdup.bam"
                duplicate_metrics = path_output + "mapping/{sample}_duplicate_metrics.txt"
        message:
                "Sort, add read group and mark duplicates of bam files: {wildcards.sample}"
        benchmark:
                path_github + "benchmark/{sample}_markdup_benchmark.txt"
        log:
                stdout = path_github + "log/{sample}_markdup.log"
        shell:
                "(cd /shared/projects/origamis/wgs/ && "
                "samtools sort -n -O BAM {input.raw_cram} | "
                "samtools fixmate -m - - | "
                "samtools sort -O BAM - | "
                "samtools markdup -s - - | "
                "samtools addreplacerg -r '@RG\tID:1\tSM:{wildcards.sample}\tPL:ILLUMINA\tLB:lib\tPU:ophrys' "
                 "-O BAM - - | "
                #"samtools view -b > {output.markdup_bam}) 2> {log.stdout}"
                "samtools view -C "
                "-T {input.reference_genome} "
                "-o {output.markdup_cram} "
                "-) 2> {output.duplicate_metrics}"

rule samtools_stats:
        input:
                markdup_cram = path_output + "mapping/{sample}_markdup.cram"
        output:
                samtools_stats = path_output + "mapping/{sample}_samtools_stats.stats"
        message:
                "Samtools stats : {wildcards.sample}"
        log:
                stdout = path_github+"log/{sample}_samtools_stats.log"
        benchmark:
                path_github+"benchmark/{sample}_samtools_stats_benchmark.txt"
        shell:
                "(samtools index "
                "{input.markdup_cram} && "
                "samtools stats {input.markdup_cram} "
                "> {output.samtools_stats}) 2> {log.stdout}"

rule samtools_flagstats:
       input:
               markdup_cram = path_output + "mapping/{sample}_markdup.cram"
       output:
               flagstat = path_output + "mapping/{sample}_flagstat.txt"
       message:
               "Flagstat : {wildcards.sample}"
       log:
               stdout = path_github + "log/{sample}_flagstat.log"
       benchmark:
               path_github + "benchmark/{sample}_flagstat_benchmark.txt"
       shell:
              "(samtools flagstat "
              "{input.markdup_cram} > "
              "{output.flagstat}) 2> {log.stdout}"

#RegSize=10000000
#cut -f1-2 /shared/projects/origamis/ref_genome/sphegodes/GCA_040285675.1_Osph-v2.3_genomic.fna.fai > /shared/projects/origamis/ref_genome/sphegodes/chrSize.txt 
#for i in "CM079930.1" "CM079931.1" "CM079932.1" "CM079933.1" "CM079934.1" "CM079935.1" "CM079936.1" "CM079937.1" "CM079938.1" "CM079939.1" "CM079940.1" "CM079941.1" "CM079942.1" "CM079943.1" "CM079944.1" "CM079945.1" "CM079946.1" "CM079947.1" "CM079948.1"; 
#  do 
#	chr=$i
#	chrL=$(grep -w "$chr" /shared/projects/origamis/ref_genome/sphegodes/chrSize.txt | cut -f2)
#	# List regions (e.g., “1000001-2000000”)
#	paste -d '-' \
#	<(seq 1 $RegSize $chrL) \
#	<(printf "$(seq $RegSize $RegSize $chrL)\n$chrL") \
#	> /shared/projects/origamis/ref_genome/sphegodes/chr_regions_$chr.txt
#	# Add chromosome prefix (e.g., “chr1:1000001-2000000”)
#	sed -i -e 's/^/'"$chr"':/' /shared/projects/origamis/ref_genome/sphegodes/chr_regions_$chr.txt
#  done
#cat /shared/projects/origamis/ref_genome/sphegodes/chr_regions_* > /shared/projects/origamis/ref_genome/sphegodes/chr_regions.txt

rule variant_calling:
       input:
               #chr_regions = "/shared/projects/origamis/ref_genome/sphegodes/chr_regions.txt",
               reference_genome = path_ref_genome,
               bam_list = "/shared/projects/origamis/wgs/bam_list"
       output:
               vcf_regions =  path_output + "vcf/{regions}.vcf.gz"
       message:
               "Variant calling for {wildcards.regions}"
       log:
               stdout = path_github + "log/{regions}_vcf.log"
       benchmark:
               path_github + "benchmark/{regions}_vcf.txt"
       shell:
               "(bcftools mpileup -Ou -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD "
               "-b {input.bam_list} "
               "-f {input.reference_genome} "
               "-r {wildcards.regions} | "
               "bcftools call -mOz -f GQ,GP -o {output.vcf_regions})"

#os.system("ls " + path_output+"/vcf/*.vcf.gz > " + path_output + "regions.vcf.text && bcftools concat -f " + path_output + "regions.vcf.text | bcftools sort -Oz - > " + path_output + "ophrys.vcf.gz")

rule filtered_only_variant:
       input:
               vcf = path_output + "ophrys.vcf.gz"
       output:
               vcf_filtered_only_variant = path_output + "ophrys_variant_only.vcf.gz"
       message:
               "Filter only variant"
       log:
               stdout = path_github + "log/filter_only_variant_vcf.log"
       benchmark:
               path_github + "benchmark/filter_only_variant_vcf.txt"
       shell:
               "(bcftools filter -g 5 {input.vcf} | "
               "bcftools view --types snps -m 2 -M 2 | "
               "bcftools view -Oz -o {output.vcf_filtered_only_variant}) 2> {log.stdout}"

rule mark_filter:
        input:
               vcf_filtered_only_variant = path_output + "ophrys_variant_only.vcf.gz"
        output:
               vcf_mark_filter = path_output + "ophrys_mark_filter.vcf.gz"
        message:
               "Mark filter"
        log:
               stdout = path_github+"log/mark_filter.log"
        benchmark:
               path_github+"benchmark/mark_filter_benchmark.txt"
        shell:
               "(bcftools index -f -t {input.vcf_filtered_only_variant} && "
               "vcftools --gzvcf {input.vcf_filtered_only_variant} --minGQ 10 --minDP 5 --stdout --recode --recode-INFO-all | "
               "bcftools filter -e 'QUAL<20.0' -s 'QUAL20' -m + | "
               "bcftools filter -e 'MQ<40.0' -s 'MQ40' -m + | "
               "bcftools filter -e 'INFO/DP>2.5*AVG(INFO/DP)' -s 'DP_high' -m + | "
               "bcftools filter -e 'F_MISSING>0.1' -s 'MISSING10' -m + | "
               "bcftools view -O z -o {output.vcf_mark_filter}) 2> {log.stdout}"

rule filtered:
        input:
               vcf_mark_filter = path_output + "ophrys_mark_filter.vcf.gz"
        output:
               vcf_filtered = path_output + "ophrys_filtered.vcf.gz"
        message:
                "Filtered"
        log:
                stdout = path_github+"log/filtered.log"
        benchmark:
                path_github+"benchmark/filtered_benchmark.txt"
        shell:
                "(bcftools filter -g 5 {input.vcf_mark_filter} | "
                "bcftools view --types snps --max-alleles 2 | "
                "bcftools view -f 'PASS,.' | "
                "bcftools view -O z -o {output.vcf_filtered}) 2> {log.stdout}"

rule get_input_iSMC:
       input:
               vcf_filtered = path_output + "ophrys_filtered.vcf.gz",
               script_make_input_iSMC = path_github + "script/iSMC/make_input_iSMC.R",
               script_create_tab = path_github + "script/iSMC/create_tab.py"
       output:
               bpp = path_github+"output/iSMC/{indiv}/{chrom}/{chrom}.bpp",
               tab = path_github+"output/iSMC/{indiv}/{chrom}/{chrom}.tab"
       message:
               "Get input for iSMC for {wildcards.indiv} {wildcards.chrom}"
       log:
               stdout = path_github+"logs/{indiv}_{chrom}_get_input_iSMC.log"
       benchmark:
               path_github+"benchmarks/{indiv}_{chrom}_get_input_iSMC.txt"
       params:
               tmp_vcf = path_github+"output/iSMC/{indiv}/{chrom}/{indiv}_filtered_passed.vcf.gz"
       shell:
               "(bcftools index -f {input.vcf_filtered} && "
               "mkdir -p "+path_github+"output/iSMC/{wildcards.indiv}/{wildcards.chrom} && "
               "bcftools view -s {wildcards.indiv} -Oz -o {params.tmp_vcf} {input.vcf_filtered} && "
               "python {input.script_create_tab} {params.tmp_vcf} {output.tab} {wildcards.chrom} && "
               "R --vanilla --slave --args {wildcards.indiv} {output.bpp} {params.tmp_vcf} {output.tab} < {input.script_make_input_iSMC}) && "
               "rm -f {params.tmp_vcf}"

rule iSMC:
       input:
               vcf_filtered = path_output + "ophrys_filtered.vcf.gz",
               bpp = path_github+"output/iSMC/{indiv}/{chrom}/{chrom}.bpp"
       output:
               decoding_label = path_github+"output/iSMC/{indiv}/{chrom}/my_dataset_diploid_decoding_labels.txt"
       message:
               "iSMC for {wildcards.indiv} and {wildcards.chrom}"
       log:
               stdout = path_github+"logs/{indiv}_{chrom}_iSMC.log"
       benchmark:
               path_github+"benchmarks/{indiv}_{chrom}_iSMC.txt"
       params:
               tmp_vcf = path_github+"output/iSMC/{indiv}/{chrom}/{indiv}_filtered_passed.vcf.gz"
       shell:
              "bcftools view -s {wildcards.indiv} -Oz -o {params.tmp_vcf} {input.vcf_filtered} && "
              "cd "+path_github+"output/iSMC/{wildcards.indiv}/{wildcards.chrom} && "
              "ismc params={input.bpp} && "
              "rm -f {params.tmp_vcf}"


rule mosdepth:
       input:
               markdup_cram = path_output + "mapping/{sample}_markdup.cram",
               reference_genome = path_ref_genome
       output:
               mosepth_output = path_output + "mosdepth/{sample}/dist.html"
       message:
               "mosdepth : {wildcards.sample}"
       log:
               stdout = path_github + "log/{sample}_mosdepth.log"
       benchmark:
               path_github + "benchmark/{sample}_mosdepth_benchmark.txt"
       shell:
              "(mkdir " + path_output + "mosdepth/{wildcards.sample} && "
              "cd " + path_output + "mosdepth/{wildcards.sample} && "
              "mosdepth -n -T 1,10,20,30 "
              "--fast-mode --by 500 "
              "-f {input.reference_genome} "
              "/shared/projects/origamis/wgs/mosdepth/{wildcards.sample}.wgs "
              "{input.markdup_cram}) 2> {log.stdout}"

rule ngsParalog:
       input:
               reference_genome = path_ref_genome,
               ngsParalog_script = "/shared/projects/origamis/software/ngsParalog/ngsParalog"
       output:
               cram_list = path_output + "mapping/cram_list",
               ngsParalog = path_output + "ngsParalog/ERR5100914.lr"
       message:
               "ngsParalog"
       log:
               stdout = path_github + "log/ngsParalog.log"
       benchmark:
               path_github + "benchmark/ngsParalog_benchmark.txt"
       shell:
               "find '/shared/projects/origamis/wgs/mapping/' "
               "-type f -name 'markdup_*.cram' > {output.cram_list} && "
               "samtools mpileup "
               "-b {output.cram_list} "
               "--reference {input.reference_genome} "
               "-q 0 -Q 0--ff UNMAP,DUP | "
               "{input.ngsParalog_script} calcLR "
               "-infile - -outfile {output.ngsParalog} "
               "-minQ 20 -minind 25 -mincov 1 -minind 1) 2> {log.stdout}"

rule getorganelle:
       input:
               fastp_R1 = path_output + "fastp/{sample}_R1_fastp.fastq.gz",
               fastp_R2 = path_output + "fastp/{sample}_R2_fastp.fastq.gz",
       output:
               mito_plasto = path_output + "getorganelle/{sample}_plastome_mitogenome/extended_K115.assembly_graph.fastg"
       message:
               "getorganelle : {wildcards.sample}"
       log:
               stdout = path_github + "log/{sample}_getorganelle.log"
       benchmark:
               path_github + "benchmark/{sample}_getorganelle_benchmark.txt"
       shell:
               "cd " + path_output + "getorganelle/ && "
               "get_organelle_from_reads.py "
               "-1 {input.fastp_R1} "
               "-2 {input.fastp_R2} "
               "-t 1 "
               "-o {wildcards.sample}_plastome_mitogenome "
               "-F embplant_pt,embplant_mt "
               "-R 10 "
               "--continue"
