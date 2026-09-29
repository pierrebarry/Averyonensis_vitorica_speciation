import os
import sys

#vcf= "/shared/ifbstor1/home/pbarry/cogediv_genome_architecture/output/iSMC/Dlabr/DlabrLi1/DlabrLi1_filtered_passed.vcf.gz"
vcf = sys.argv[1]
tab = sys.argv[2]
chrom = sys.argv[3]
header = os.popen("bcftools view -h "+vcf).read().split("\n")[:-1]

for i in header:
  if i.split("=")[0]=="##contig":
    length = i.split("=")[-1:][0].split(">")[0]
    chr = i.split("=")[2].split(",")[0]
    if chr==chrom:
      with open(tab, 'a') as f:
        _ = f.write(str(chr) + "\t" + str(1) + "\t" + str(length) + "\t" + str(0) + "\t" + str(int(length) - 1) + "\t" + str(50000) + "\t" + str(int(length) - 50000) + "\n")
