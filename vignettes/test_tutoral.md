# Load the package and tutorial data

```{r}
source('function.R')
source('scART.R')
load('test_step.Rdata')

```

# Create the scART object, impute missing value and filter data
```{r message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE}
art<-CreatescART(ncounts,metadata = anno)  #metadata = annotation)
art<-RunImputation(art,k=1)
art<-SparseFilter(art, ncell=2, ncell2=0.8, ncell3=2, nbin=10)
```
![image-20210423151512646](../image/image-20210423151512646-1620016637304.png)

# Reduce Dimensionality 

```{r include=FALSE}
art<-RunSim(art)
art<-DimReduce(art)
```

![image-20210503123815138](../image/image-20210503123815138.png)

```{r}
art<-RunTSNE(art, nSV=40, ndims=2, perplexity=30,seed.use = 10)
art<-RunUMAP(art ,seed.use = 10,nSV=40)
```
# Group cells into clusters, you can take a good look at the output pdf to adjust 'rho_cutoff' and 'delta_cutoff'
```{r message=FALSE, warning=FALSE, include=FALSE, paged.print=FALSE}
art<-RunCluster(art,delta_cutoff = 4,rho_cutoff = 8)
```
![image-20210503123838311](../image/image-20210503123838311.png)

# Visualize an Embedding

```{r}
Visualization_2D(art,reductions = 'TSNE',anno='type2',size = c(1,1))
```

![image-20210503123859393](../image/image-20210503123859393.png)

```{r}
Visualization_2D(art,reductions = 'UMAP',anno='type2',size = c(1,1))
```

![image-20210503123904486](../image/image-20210503123904486.png)

# Run and plot trajectory

```{r}
art<-RunTrajectory(art,nSV = 20,ndim =10,gamma = 10)
plotTrajectory(art,anno='type2')
```

![image-20210503123950798](../image/image-20210503123950798.png)

# Create cell-gene matrix and explore gene chromtin accessibility

```{r}
art<-MapBin2Gene(art,Org = 'hg19')
PlotSelectGenesATAC(art,gene2plot = c('OTOAP1','PNPLA5'))  

```
![image-20210503123959699](../image/image-20210503123959699.png)



```{r}
art<-RunChromVAR(art,Org=c('hg19'),species = c("Homo sapiens") ,min.count=10)
### min.count:the threshold of a peaks found at at least 10 cells
PlotSelectTF(art,)
```

![image-20210503124005433](../image/image-20210503124005433.png)




```{r}
seurat=art2seurat(art)
snap=art2snap(art)
```