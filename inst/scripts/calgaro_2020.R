## Samples for Calgaro 2020 study
library(magrittr)

## Import data from Calgaro repository
load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/16Sdatasets_for_replicability_filtered.RData?raw=true"))

## get rownames and colnames
samples <- phyloseq::sample_names(ps_list_16S[["Subgingival_Supragingival"]])
taxa <- phyloseq::taxa_names(ps_list_16S[["Subgingival_Supragingival"]])

## Import the data from the HMP16SData package
se <- HMP16SData::V35()
se <- se[taxa,samples]
se

tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = SimpleList(counts = SummarizedExperiment::assay(se)),
  colData = SummarizedExperiment::colData(se),
  rowData = SummarizedExperiment::rowData(se),
  rowTree = S4Vectors::metadata(se)[["phylogeneticTree"]]
)

S4Vectors::metadata(se)[["phylogeneticTree"]]
