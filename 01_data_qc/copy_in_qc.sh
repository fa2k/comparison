for dir in \
    ../../../../demultiplexed/completed/1906/190606_E00426_0063_AHYK7TCCXY/Data/Intensities/BaseCalls/QualityControl/Ribarska-DNAlibs-2019-04-01 \
    ../../../../demultiplexed/completed/1906/190613_E00426_0064_AHYKCHCCXY/Data/Intensities/BaseCalls/QualityControl/Ribarska-DNAlibs-2019-04-01 \
    ../../../../demultiplexed/completed/1906/190626_E00396_0021_BHYNTHCCXY/Data/Intensities/BaseCalls/QualityControl/Ribarska-DNAlibs-2019-04-01 \
    ../../../../demultiplexed/completed/1912/191203_J00146_0133_AHFV2HBBXY/Data/Intensities/BaseCalls/QualityControl/Ribarska-SwiftLibs2-2019-11-22
do
    cp -r $dir/*/ QualityControl/
done


