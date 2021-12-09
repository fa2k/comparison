library(tidyverse)


library(paletteer)
paletteer_d("basetheme::minimal")


df <- read_csv('summary-tables/main-library-metrics.csv')

# Tapestation
tape <- read_csv('summary-tables/tapestation-insert-sizes.csv')
df <- left_join(df, tape, by=c('Kit', 'Conc', 'Replicate'))

#load.ins.size.hist <- function(..., ins.path) {
load.ins.size.hist <- function(..., ins.path, Kit, Replicate, Conc) {
  dp <- read_delim(ins.path, skip=10) %>%
    rename(reads = `All_Reads.fr_count`) %>%
    mutate(
      Kit=Kit, Conc=Conc, Replicate=Replicate
    )
}

df %>%
  mutate(ins.path = str_glue("d30_downsample/{LIBRARY}_DS_MD.InsertSizeMetrics.txt")) %>%
  pmap_dfr(load.ins.size.hist) ->
    insert.sizes

in1 <- insert.sizes %>% filter(Conc == '100ng') %>% filter(insert_size >= 25)

ggplot(in1, aes(x=insert_size, y=reads, color=Kit)) +
  geom_line()


ggplot(df, aes(x=MEAN_INSERT_SIZE, y=TapeInsertSize, colour=Kit)) + geom_point()


# *** Variant calling performance vs Insert size ***

vc_ins <- read_csv('summary-tables/b_analysis-merged_f1_df-dataset_a.csv')

ggplot(vc_ins, aes(x=MEAN_INSERT_SIZE, y=F1_score, colour=Kit, shape=Type)) +
  geom_point() + 
  scale_color_paletteer_d("basetheme::minimal")


#stat_summary(fun.data= mean_cl_normal) + 
  #geom_smooth(method='lm')



