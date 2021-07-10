GRN

# filter gmat

```
gmat=FilterGmat(art,hvg=4000)
```

# GRNboost

#path.to.pyscenic: 系统中pyscenic的位置，如/home/lijingyu/anaconda3/bin/pyscenic

#tf_file ：tf_file   的位置'hs_hgnc_tfs.txt'

#only_pos 只保留正调控

#min_importance： 最小的阈值

```
regulon=runGRNboost(gmat,path.to.pyscenic,tf_file=' ',min_importance=2,only_pos=TRUE)
```





# GRN.R



