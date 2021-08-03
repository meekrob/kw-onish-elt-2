
# plot cebp-1 RNA-seq read counts
```{r fig.width=3, fig.height=3}
# cebp1_rnaseq_counts <- dineen_nishimura_counts %>% 
#   filter(WBGeneID == "WBGene00016997") %>% 
#   pivot_longer(colnames(dynamic_counts_matrix_scaled_ascend)) %>% 
#   separate(name, sep = "_", into = c("sample", "sorted", "rep")) %>% 
#   mutate(value = as.numeric(value), sample = fct_relevel(sample, c("wt", "elt7D", "elt2D", "elt2Delt7D"))) %>%
#   ggplot(aes(x = sample, y = value)) + 
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(stat = "identity", shape = 21, size = 3, fill = "grey", width = 0.3, height = 0) + 
#   theme_classic() +
#   xlab("genotype")+
#   ylab("normalized count")

cebp1_rnaseq_counts <- dineen_nishimura_counts %>%
  filter(WBGeneID == "WBGene00016997") %>%
  pivot_longer(colnames(dynamic_counts_matrix_scaled_ascend)) %>%
  separate(name, sep = "_", into = c("sample", "sorted", "rep")) %>%
  filter(sample %in% c("wt", "elt2D")) %>%
  mutate(value = as.numeric(value), sample = fct_relevel(sample, c("wt", "elt2D"))) %>%
  ggplot(aes(x = sample, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(stat = "identity", shape = 21, size = 3, fill = "grey", width = 0.25, height = 0) +
  theme_classic() +
  xlab("genotype")+
  ylab("normalized count")
cebp1_rnaseq_counts
if (do_plot == TRUE){
  myggsave(plot = cebp1_rnaseq_counts, 
           name = "Figure5B_CEBP1_RNAseq_ReadCounts",
           height = 3,
           width = 3,
           plotdir = plotdir)
}
```
```{r}
dineen_nishimura_counts %>% filter(WBGeneID == "WBGene00017687")%>% pivot_longer(colnames(dynamic_counts_matrix_scaled_ascend)) %>% separate(name, sep = "_", into = c("sample", "sorted", "rep")) %>% mutate(value = as.numeric(value), sample = fct_relevel(sample, c("wt", "elt7D", "elt2D", "elt2Delt7D"))) %>%
  ggplot(aes(x = sample, y = value)) + geom_point(stat = "identity")

```
