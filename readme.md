# ddqc - Biology-centered data-driven quality control for single cell/nucleus RNA sequencing
## Required packages
- numpy>1.20
- matplotlib>=3.4.0
- pandas>=1.2.0
- pegasusio
- pegasuspy>=1.3 
- seaborn>=0.11
## installation
1. Clone this repository using the following code:
   
   `git clone https://github.com/ayshwaryas/ddqc.git`
2. In the root directory of the package run the following command:
   
    `pip install .`
3. For the usage instructions refer to tutorials/ddqc_tutorial.ipynb
## FAQs
### What is the best input for ddqc?
ddqc takes pegasus MultimodalData object as an input. The call of ddqc.ddqc_metrics ideally 
should be the next step after reading the data, similar to when regular QC is done in scRNA-seq pipelines. 
ddqc cant work with the normalized matrix, so it should be performed before normalization step.

### What data formats ddqc supports?
ddqc can work with all formats supported by Pegasus. This includes h5ad, h5, mtx, csv, loom. <br>
Here is a sample code for reading different data formats:
```python
import pegasusio as io
import ddqc
# read h5
data1 = io.read_input("path/file.h5ad", genome = 'hg19')
# read csv
data2 = io.read_input("path/file.csv", genome = 'hg19')
# read mtx
data3 = io.read_input("path/file.mtx", genome = 'hg19')
# call ddqc
ddqc.ddqc_metrics(data1)
```

Pegasus can also aggregate multiple files into one object. To do it, first create a CSV file with the information about your data:
```text
Sample,Location
sample1,path/file1.mtx
sample2,path/file2.mtx
sample3,path/file3.mtx
```

Then use the following Python code:
```python
import pegasusio as io
import ddqc
data = io.aggregate_matrices("pegasusio_test_cases/case6/count_matrix.csv")
# call ddqc
ddqc.ddqc_metrics(data)
```
Please refer to [PegasusIO tutorial](https://pegasusio.readthedocs.io/en/stable/_static/tutorials/pegasusio_tutorial.html) for a complete guide on reading files in Pegasus.

### What are the outputs of ddqc?
There are four plots provided for exploratory data analysis:

* Two boxplots:
  * log2(n_genes) by cluster: shows log2 of number of genes for each cluster in the initial clustering. Red line at 200 genes (7.64 in log2 scale) represents the most common fixed threshold cutoff for n_genes.
  * percent_mito by cluster: shows percent_mito for each cluster in the initial clustering. Red line at 10% represents the most common fixed threshold cutoff for percent_mito.
* If you are using "mad" method and didn't use any of metric threshold overwrites (threshold_counts, threshold_genes, threshold_mito, threshold_ribo), ddqc will generate two facet plots that show how number of cells that are filtered out changes depending on the threshold value. These plots will help you to pick a threshold parameter if you want to tune it. Each of these plots is faceted by cluster, and has a threshold parameter (from 1 to 3) on x-axis. First plot will have a number of cells filtered out in cluster on y-axis, second plot will have a percentage of cells filtered out in cluster.


If you requested to return df_qc the function will return a pandas dataframe containing the following info for each cell:
* metric: QC metric number
* cluster_labels: cluster from initial clustering performed by ddqc
* metric_lower_co and metric_upper_co: lower and upper cuttofs for each metric on which ddqc was performed. If ddqc was not performed for upper or lower end of this metric this field will be None
* metric.passed.qc: whether the cell passed qc for a given metric This information is useful if you want to understand based on which metric the cell was filtered out.

The pegasus object will have the following data added to its obs field:
* n_counts: number of counts in each cell
* n_genes: number of genes in each cell
* percent_mito: percentage of mitochondrial reads in each cell
* percent_ribo: percentage of ribosomal reads in each cell
* passed_qc (true or false): whether the cell passed the quality control