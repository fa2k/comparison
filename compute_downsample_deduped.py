import pandas as pd

samples = list(set(pd.read_table("sampleList.txt")['sampleName']))
readCounts = []
for sampleName in samples:
    with open("20_piccard_dd/{}.readCount.txt".format(sampleName)) as f:
        readCount = float(f.read().strip())
    readCounts.append(readCount)

minReadCount = min(readCounts)
print("Minimum number of deduplicated reads: ", minReadCount)
normalizer = (minReadCount // 5000000) * 5000000
print("Normalized number of reads to use: ", normalizer)

# Check this read count against fastq file content: 24399976

data = []
for sampleName, readCount in zip(samples, readCounts):
    bamfile = "20_piccard_dd/{sampleName}.filter.bam".format(sampleName=sampleName)
    data.append({'sampleName': sampleName + '-dd', 'fileBAM': bamfile, 'downSampleFactor': normalizer / readCount})

df = pd.DataFrame(data)
df.to_csv("downsample_factors_dd_fix.txt", index=False, sep="\t")

