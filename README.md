# spatial-metabolomics

This repository contains a collection of R scripts for processing and analyzing Mass Spectrometry Imaging (MSI) data, specifically designed for working with imzML format files. The pipeline provides a comprehensive set of functions for data loading, quality control, preprocessing, and visualization.

## Overview

The scripts facilitate the complete workflow for MSI data analysis, from raw data import to advanced statistical analysis and visualization. The whole analysis was based on Cardinal v3.xx, please confirm your package version.

## Notes:

0. Please check your MSI platform at first, record all the processing steps and parameters of MSI. Remember, batch effects, batch effects, batch effects🤷
1. Data Import: it is supposed .imzML files were the default file format, please confirm you have the ## writable permission ##, or it would not be loaded to RAM, which could cause some unpleasant errors.
2. Processing: the whole processing steps were based on Cardinal package, including feature selection (TIC normalization, peak picking, recalibrating, etc.) and representative marker finding (PCA and SSC (Spatial Shrunken Centroid)). Please refer to [`Cardinal`](https://www.bioconductor.org/packages/release/bioc/html/Cardinal.html) and [`CardinalWorkflows`](https://www.bioconductor.org/packages/release/data/experiment/html/CardinalWorkflows.html) documents in Bioconductor if you have any question.
3. Others: the following steps, usually meaning differential analysis & enrichment, were not shown here due to no standard steps acknowledged. For biomarker finding, some processes like `FindAllMarkers` in `Seurat` or likewise could be adopted. and for the other, `mmumichog` may be better for the feature in MSI, which actually `m/z` format, not a single metabolite annotation. of course, if your data has MS2 and even MS3/4 or plus LC/GC, or a quite clear compound candidate, it is much easier, just do what bulk metabolomics does.

