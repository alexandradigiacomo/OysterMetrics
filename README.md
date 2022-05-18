# Oyster Metrics

This repository contains the code to extract structural metrics from UAS- and LiDAR- derived point clouds of oyster reef surfaces. This work was performed in 2020-2021 in the Marine Robotics and Remote Sensing Laboratory at Duke University Marine Laboratory. The 'OysterMetrics.py' script contained in this repository interface with ArcGIS to create user-friendly GUIs that take in point cloud files (.csv, .txt, .las) and export structural metrics (e.g. rugosity) of the reef around the reef edge and reef crest. The computed structural metrics are rugosity, slope, aspect, profile and planform curvature, and fractal dimension. 

The .Rmd files in this repository contain the data analysis scripts specific to this project. Three reefs (1L5, CCA3, and Natural) that represent different reef development stages were examined. The 'smoothedAnalysis.Rmd' file combines the three individual reef metrics and performs statistical comparisons. 

