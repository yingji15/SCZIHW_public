The raw data for features (BRAINSPAN, DEPICT, FANTOM5, LAKE) are from the following sources
The processed data is available in Zenodo repo

* BRAINSPAN
https://www.brainspan.org/static/download.html
We added column names, and log2(value+1) transformed data
```
row.meta <- read.table("rows_metadata",sep="\t",head=TRUE) %>% tbl_df
col.meta <- read.table("columns_metadata",sep="\t",head=TRUE) %>% tbl_df
exprs.matrix <- read.table("expression_matrix",sep="\t") %>% tbl_df

names(exprs.matrix)[1]='row_num'
exprs.matrix[2:ncol(exprs.matrix)]=log(exprs.matrix[2:ncol(exprs.matrix)]+1,2)
```

* DEPICT
https://data.broadinstitute.org/mpg/depict/depict_download/reconstituted_genesets/GPL570-GPL96-GPL1261-GPL1355TermGeneZScores-MGI_MF_CC_RT_IW_BP_KEGG_z_z.txt
Gene symbols are added, and the dataframe was saved in RDS format, no other treatments


* FANTOM5
https://fantom.gsc.riken.jp/5/datafiles/latest/extra/gene_level_expression/hg19.gene_phase1and2combined_tpm.osc.txt.gz
Here is the R script to process, mostly just added gene symbols
```
human<-read.delim("hg19.gene_phase1and2combined_tpm.osc.txt.gz",skip=1832,as.is=T)

colnames(human)[-1]<-substr(colnames(human)[-1],nchar(colnames(human)[-1])-20,nchar(colnames(human)[-1])-12)
colnames(human)[1]<-"gene_symbol"
human<-human[,apply(human,2,function(x) !all(is.na(x)))]
human[,-1]<-round(human[,-1],2)

write.table(human,"human_fantom5_expression.txt",quote=F,row.names=F,sep="\t")
```


* LAKE
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930
To reduce the dimension of this dataset, we take an average of all the expression in each cell type and state (with labels provided in the dataset), resulting in a matrix with 61 columns.
Here are some bash script from Rui Chen using Perl to process the data.
```
# for each region, remove the double quote and add the gene column in header
zcat GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz | sed 's/"//g' | perl -lane 'if($.==1){print "gene\t$_";}else{print}'> Cerebellar.count

#Identify the common genes (intersection of the three sets) in all three regions.
cat <(cut -f 1 Cerebellar.count | sed '1d') <(cut -f 1 FrontalCortex.count | sed '1d') <(cut -f 1 VisualCortex.count | sed '1d') | sort | uniq -c | perl -lane 'print $F[1] if $F[0]==3' | sed '1i\gene' > allcommgeneid

#Identify the union of the three sets
cat <(cut -f 1 Cerebellar.count | sed '1d') <(cut -f 1 FrontalCortex.count | sed '1d') <(cut -f 1 VisualCortex.count | sed '1d') | sort | uniq | sed '1i\gene'> allgeneid

##combine the single cell samples/columns to cell type clusters.
parallel cat {} \| perl -MList::Util=sum -lane \''@o=();if($.==1){ map { @tmp=split/_/,$F[$_];push @{$a{$tmp[0]}},$_; } 1..$#F;@sorted_key=sort keys %a;print "$F[0]\t".join ("\t",@sorted_key);next} for $key(@sorted_key){push @o,sum(@F[@{$a{$key}}]);}print $F[0],"\t",join ("\t",@o)'\' \> {}.region ::: *count

parallel cat {} \| perl -MList::Util=sum -slane \''@o=();if($.==1){ map { @tmp=split/_/,$F[$_];push @{$a{$tmp[0]}},"$_"; } 1..$#F;@sorted_key=map {$x."_".$_} sort keys %a;print "$F[0]\t".join ("\t",@sorted_key);next} for $key(sort keys %a){push @o,sum(@F[@{$a{$key}}])/scalar(@{$a{$key}});}print $F[0],"\t",join ("\t",@o)'\' -- -x={} \> {}.region.mean ::: *count

##Coordinate all the files using intersection or union of gene id.
parallel comp2line.hash.pl -q {2} -c 1 -d 1 -db {1} -e \| cut -f 1,3- \> {1}.{2} :::  *region.mean ::: allcommgeneid allgeneid

##combine three tissues into  one single file
##The reads counts in each single cell (column) of intersection of gene id (row) of all three tissues
paste *.count.region.mean.allcommgeneid  | perl -lane 'if($.==1){do {push @h,$_ if $F[$_]!~/gene/} for 1..$#F;print join "\t",@F[(0,@h)];next}print join "\t",@F[(0,@h)]' > all.count.region.mean.allcommgeneid


```

