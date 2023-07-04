MultiAssayExperiment objects in the R programming language are used to store many assays all relatiing to one phenotype, project, goal, etc. - please see https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html

However, presently the tools used for creating these objects, such as those housed by cBioPortalData are crude and frequently fail because they tend to be based on receiving appropriately formatted files as input.

By contrast, the MaeUtils repository is founded on the tautology that the internal structure and logic of a Mae object provides substantial information as to whether a set of flatfiles is appopriate for instantiation of the same.

Specifically, we demonstrate more robust functionality by devoting a series of heuristics to unambiguous identification of the primaryId column that lays within the file best suited to become colData(MAE). The functions will attempt to identify such a column even in the absence of exact matches to the strings in a candidate column. 

The first functionality (i.e., FlatFilesToMae.R) attempts to automate creation of a Mae objects given only a directory containing (appropriate) flat files as input. Future functonality to enchance ease of use is planned for the future.
