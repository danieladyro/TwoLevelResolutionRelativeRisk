# TwoLevelResolutionRelativeRisk
This repository contains the data for the paper "Two-level resolution of relative risk of dengue disease in a hyperendemic city of Colombia"
The folder contains the following files:
DengueData_PlosONE.R: R code to reproduce all the results and figures in the paper.
DengueData.RData: binary .RData containing two dataframes, DengueData_INLA (which contains the observed and expected data for the 94 census sectors and 91 epidemiological periods), and cuminc (which contains average annual incidence cases per 100.000 inhabitants and age-group).
Dengue_INLA_TLModelA_TypeIV.RData: binay .RData containing the fitted model with INLA (TL-Model A, Type IV interaction).
comunasGraph.dat: represents the spatial neighborhood matrix the communes.
sectoresGraph.dat: represents the spatial neighborhood matrix the census sectors.
comunaNP: Folder containing the shapefile of the communes.
sector293NP: Folder containing the shapefile of the census sectors.
