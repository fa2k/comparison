import pandas as pd

samples = list(set(pd.read_table("sampleList.txt")['sampleName']))
readCounts = []

data = []

print("Computing downsample factors based on raw reads")
for sampleName in samples:
    with open("20_piccard/{}.AlignmentSummaryMetrics.txt".format(sampleName)) as f:
        for line in f:
            if line.startswith("FIRST_OF_PAIR\t"):
                readCount = float(line.split("\t")[1])
    readCounts.append(readCount)

minReadCount = min(readCounts)
print("Minimum number of reads: ", minReadCount)
normalizer = (minReadCount // 10000000) * 10000000
normalizer = int(90e6)
print("Normalized number of reads to use: ", normalizer)

for sampleName, readCount in zip(samples, readCounts):
    bamfile = "/data/runScratch.boston/analysis/projects/Project_Ribarska-DNAlibs-2019-04-01/20_piccard/{sampleName}.bam".format(sampleName=sampleName)
    data.append({'sampleName': sampleName + '-aa', 'fileBAM': bamfile, 'downSampleFactor': min(1.0, normalizer / readCount)})
    if normalizer / readCount > 1: print("Warning: sample {} only has {} reads, set to 1.0.".format(sampleName, readCount))

df = pd.DataFrame(data)
output_file = "downsample_factors.txt"
df.to_csv(output_file, index=False, sep="\t")
print("Results written to", output_file)
