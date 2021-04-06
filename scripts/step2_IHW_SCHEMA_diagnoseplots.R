# IHW diagnostic plots
# 4/5/21
# /scratch/cgg/jiy1/RF2021/IHW_2021

# "IHW_sczSCHEMA.040121.RData"

ihw01 <- ihw(p=tab32$P.meta,covariates=tab32$bsp_fan_dep_lake_scaleavg,alpha=0.1)

weights(ihw01, levels_only = TRUE)
#           [,1]       [,2]      [,3]      [,4]      [,5]
 [1,] 0.2609468 0.24437260 0.1823477 0.2271690 0.3921355
 [2,] 0.2609468 0.19976490 0.1823477 0.2271690 0.3921355
 [3,] 0.2609468 0.00000000 0.0000000 0.2271690 0.3921355
 [4,] 0.2609468 0.24017042 0.2192304 0.2271690 0.3921355
 [5,] 0.2609468 0.00000000 0.2192304 0.2271690 0.3921355
 [6,] 0.2609468 0.00000000 0.2192304 0.2271690 0.3921355
 [7,] 0.7662724 0.00238231 0.6444441 0.6670835 0.7773004
 [8,] 1.3071289 0.26602996 0.2428353 1.4252342 0.7773004
 [9,] 1.3071289 0.22206875 2.2728618 1.4252342 0.7773004
[10,] 3.0109468 2.77098312 2.7057101 2.9771690 3.1421355
[11,] 3.0109468 7.69321136 4.5144345 2.9771690 3.1421355

library(ggplot2)
pdf("ihw_schema.alpha01.pdf")
plot(ihw01,what="decisioboundary")
gg <- ggplot(as.data.frame(ihw01), aes(x = pvalue, y = adj_pvalue, col = group)) + 
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
plot( gg%+% subset(as.data.frame(ihw01), adj_pvalue <= 0.1))

rbind(data.frame(pvalue = tab32$P.meta, covariate = rank(tab32$bsp_fan_dep_lake_scaleavg)/nrow(tab32), 
                 covariate_type="avg 4pred")) %>%
ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) + 
   ylab(expression(-log[10]~p))



tab32$Group <- groups_by_filter(tab32$bsp_fan_dep_lake_scaleavg, 11)

ggplot(tab32, aes(x=P.meta)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ Group, nrow = 3)   
  
ggplot(tab32, aes(x = P.meta, col = Group)) + stat_ecdf(geom = "step") 
dev.off()



library(ggpubr)

pdf("ihw_schema.alpha01.scatter.pdf",width = 8,height = 4)

a = rbind(data.frame(pvalue = tab32$P.meta, covariate = rank(tab32$bsp_fan_dep_lake_scaleavg)/nrow(tab32), 
                 covariate_type="avg 4pred")) %>%
ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) + 
   ylab(expression(-log[10]~p))

tab32$Group2 <- groups_by_filter(tab32$bsp_fan_dep_lake_scaleavg, 3)

b = ggplot(tab32, aes(x = P.meta, col = Group2)) + stat_ecdf(geom = "step") 

ggarrange(a,b, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)
dev.off()



#tab32$Group <- groups_by_filter(tab32$bsp_fan_dep_lake_scaleavg, 11)
pdf("ihw_schema.alpha01.hist11.pdf")
ggplot(tab32, aes(x=P.meta)) + 
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ Group, nrow = 3)   
dev.off()
 

dev.off()


