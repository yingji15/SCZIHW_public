# IHW diagnostic plots
# 4/5/21



# input for ihw: p-values and covariates
ihw01 <- ihw(p=tab32$P.meta,covariates=tab32$bsp_fan_dep_lake_scaleavg,alpha=0.1)

# look at weights at different strata and folds
weights(ihw01, levels_only = TRUE)

library(ggplot2)



# Stratified p-value histograms
tab32$Group <- groups_by_filter(tab32$bsp_fan_dep_lake_scaleavg, 11)

ggplot(tab32, aes(x=P.meta)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ Group, nrow = 3)   
  
ggplot(tab32, aes(x = P.meta, col = Group)) + stat_ecdf(geom = "step") 
dev.off()





