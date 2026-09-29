import numpy as np
import json
import sys
import os

#python /shared/projects/origamis/wgs/script/extract_fastp.py /shared/projects/origamis/wgs/fastp
#sys.argv=["","/shared/projects/origamis/wgs/fastp"]
files = os.listdir(sys.argv[1])
# Filtering only the files.
files = [f for f in files if os.path.isfile(sys.argv[1]+'/'+f)]

files = [x for x in files if  "json" in x]

result = ["sample",
"total_reads_before",
"total_bases_before",
"gc_content_before",
"q20_rate_before",
"q30_rate_before",
"passed_filter_reads",
"corrected_reads",
"low_quality_reads",
"too_many_N_reads",
"too_short_reads",
"low_complexity_reads",
"total_reads_after",
"total_bases_after",
"gc_content_after",
"q20_rate_after",
"q30_rate_after",
"duplication_rate"
]
with open(sys.argv[1]+'/Summary_fastp.txt',"w") as f:
	np.savetxt(f, [result],newline='\n', fmt='%s', delimiter=";")

for i in files:
	with open(sys.argv[1]+"/"+i) as json_fastp:
		a = json.load(json_fastp)
	
	result=[]
	
	result+=[i.split("fastp_report_")[1].split(".json")[0]]
	## total reads_before_filtering
	result+=[a['summary']['before_filtering']['total_reads']]
	result+=[a['summary']['before_filtering']['total_bases']]
	result+=[a['summary']['before_filtering']['gc_content']]
	result+=[a['summary']['before_filtering']['q20_rate']]
	result+=[a['summary']['before_filtering']['q30_rate']]
	
	## filtering result
	
	reads_inital=a['summary']['before_filtering']['total_reads']
	base_inital=a['summary']['before_filtering']['total_bases']
	result+=[a['filtering_result']['passed_filter_reads']/reads_inital*100]
	result+=[a['filtering_result']['corrected_reads']/reads_inital*100]
	result+=[a['filtering_result']['corrected_bases']/base_inital*100]
	result+=[a['filtering_result']['low_quality_reads']/reads_inital*100]
	result+=[a['filtering_result']['too_many_N_reads']/reads_inital*100]
	result+=[a['filtering_result']['too_short_reads']/reads_inital*100]
	result+=[a['filtering_result']['low_complexity_reads']/reads_inital*100]
	
	## total after_filtering
	
	result+=[a['summary']['after_filtering']['total_reads']]
	result+=[a['summary']['after_filtering']['total_bases']]
	result+=[a['summary']['after_filtering']['gc_content']]
	result+=[a['summary']['after_filtering']['q20_rate']]
	result+=[a['summary']['after_filtering']['q30_rate']]
	result+=[a['duplication']['rate']]
	
	with open(sys.argv[1]+'/Summary_fastp.txt',"a") as f:
		np.savetxt(f, [result],newline='\n', fmt='%s', delimiter=";")

