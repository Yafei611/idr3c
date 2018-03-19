# 3-components IDR

*idr3c* is a semi-parametric model to determine differential expressions using information across different sources. It can be applied to the integration of RNA-seq and microarray data. By incorporating both the significance of differential expression and the consistency across platforms, our method effectively detects differential expressions with moderate but consistent signals. 

## Installing

The package can be installed from github. R package **devtools** is required.

```
library(devtools)
install_github("Yafei611/idr3c")
library(idr3c)
```

## Example

```
x = expr
plot(x[,1],x[,2],xlim=c(-8,8),ylim=c(-8,8),cex = .4);
idr.out = IDR.3component(x = x)
```

## Citation

Please cite the following paper: A semi-parametric statistical model for integrating gene expression profiles across different platforms: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0847-y

