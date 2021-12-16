library(tidyverse)


library(paletteer)
paletteer_d("basetheme::minimal")


df <- read_csv('summary-tables/main-library-metrics.csv')

# Tapestation
tape <- read_csv('summary-tables/tapestation-insert-sizes.csv')
df <- left_join(df, tape, by=c('Kit', 'Conc', 'Replicate'))

#load.ins.size.hist <- function(..., ins.path) {
load.ins.size.hist <- function(..., ins.path, Kit, Replicate, Conc) {
  dp <- read_delim(ins.path, skip=10, delim='\t') %>%
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





# Compare mean insert size Tapestation vs by sequencing reads
model <- lm(MEAN_INSERT_SIZE - TapeInsertSize ~ 1, data=df)
print(summary(model))

ggplot(df, aes(x=TapeInsertSize, y=MEAN_INSERT_SIZE, colour=Kit, shape=Conc)) +
  geom_point(size=3) +
  scale_color_paletteer_d("basetheme::minimal") +
  geom_abline(slope=1,
              intercept=coef(model)[["(Intercept)"]],
              linetype='dashed') +
  ylab('Mean insert size by sequencing reads') +
  xlab('Mean insert size by Tapestation') +
  theme_classic()
ggsave('plots-suppl/insert-size-tape-vs-reads.pdf')



# Pairwise comparison MEAN_INSERT_SIZE between kits / conc condition
get.pairwise.pvalue <- function(group1, group2) {
  t.test(group1, group2, var.equal = FALSE)
}

grouped <- df %>% group_by(Kit, Conc)

get.get.comparisons <- function(concs) {
  get.comparisons <- function(group.1, key) {
      pvalues <- grouped %>% 
                 filter(Kit != key$Kit && Conc == concs) %>% 
                 group_modify(
                    ~tibble(Other = .y$Kit, p.value=
                                   t.test(
                                      group.1$MEAN_INSERT_SIZE, 
                                      .x$MEAN_INSERT_SIZE, 
                                      equal.var=FALSE)$p.value
                  )
                 )
        ungroup(pvalues) %>% slice_max(p.value, n=1) %>% select(p.value, Other)
    }
}
crossp <- grouped %>% filter(Conc=='10ng') %>% group_modify(get.get.comparisons('10ng'))
print(crossp)
crossp <- grouped %>% filter(Conc=='100ng') %>% group_modify(get.get.comparisons('100ng'))
print(crossp)


# *** Variant calling performance vs Insert size ***

vc_ins.a <- read_csv('summary-tables/b_analysis-merged_f1_df-dataset_a.csv')

print('--- SNP MODEL ---')
f1.model.snp <- lm(F1_score ~ MEAN_INSERT_SIZE, data=vc_ins.a %>% filter(Type=='SNP'))
print(summary(f1.model.snp))
print('--- INDEL MODEL ---')
f1.model.indel <- lm(F1_score ~ MEAN_INSERT_SIZE, data=vc_ins.a %>% filter(Type=='INDEL'))
print(summary(f1.model.indel))

ggplot(vc_ins.a, aes(x=MEAN_INSERT_SIZE, y=F1_score, colour=Kit, shape=Conc)) +
  geom_point(size=5) + 
  scale_color_paletteer_d("basetheme::minimal") +
  geom_abline(slope=coef(f1.model.snp)[["MEAN_INSERT_SIZE"]], intercept=coef(f1.model.snp)[["(Intercept)"]]) +
  facet_grid(row=vars(Type)) + 
  theme_classic()


vc_ins.d <- read_csv('summary-tables/b_analysis-merged_f1_df-dataset_d.csv')

print('--- SNP MODEL ---')
f1.model.snp <- lm(F1_score ~ MEAN_INSERT_SIZE, data=vc_ins.d %>% filter(Type=='SNP'))
print(summary(f1.model.snp))
print('--- INDEL MODEL ---')
f1.model.indel <- lm(F1_score ~ MEAN_INSERT_SIZE, data=vc_ins.d %>% filter(Type=='INDEL'))
print(summary(f1.model.indel))

ggplot(vc_ins.d, aes(x=MEAN_INSERT_SIZE, y=F1_score, colour=Kit, shape=Conc)) +
  geom_point(size=5) + 
  scale_color_paletteer_d("basetheme::minimal") +
  geom_abline(slope=coef(f1.model.snp)[["MEAN_INSERT_SIZE"]], intercept=coef(f1.model.snp)[["(Intercept)"]]) +
  facet_grid(row=vars(Type)) + 
  theme_classic()


# *** Mean coverage significant differences ***

get.get.comparisons.cov <- function(concs) {
  get.comparisons.cov <- function(group.1, key) {
    pvalues <- grouped %>% 
      filter(Kit != key$Kit && Conc == concs) %>% 
      group_modify(
        ~tibble(Other = .y$Kit, p.value=
                  t.test(
                    group.1$MEAN_COVERAGE, 
                    .x$MEAN_COVERAGE, 
                    equal.var=FALSE)$p.value
        )
      )
    ungroup(pvalues) %>% slice_max(p.value, n=1) %>% select(p.value, Other)
  }
}
crossp <- grouped %>% filter(Conc=='10ng') %>% group_modify(get.get.comparisons.cov('10ng'))
print(crossp)
crossp <- grouped %>% filter(Conc=='100ng') %>% group_modify(get.get.comparisons.cov('100ng'))
print(crossp)
t.test(df %>% filter(Conc=='100ng' & Kit=='Nextera') %>% select(MEAN_COVERAGE), 
       df %>% filter(Conc=='100ng' & Kit!='Nextera') %>% select(MEAN_COVERAGE), equal.var=FALSE)

t.test(df %>% filter(Conc=='10ng' & Kit %in% c('Swift2S', 'Nextera')) %>% select(MEAN_COVERAGE), 
       df %>% filter(Conc=='10ng' & !(Kit %in% c('Swift2S', 'Nextera'))) %>% select(MEAN_COVERAGE), equal.var=FALSE)
