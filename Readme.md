MultiAssayExperiment objecs are objects in the R programming language used to store many assays relation to one phenotype, project, goal, etc. - please see https://bioconductor.org/packages/release/bioc/html/MultiAssayExperiment.html

However, presently the tools used for creating these objects, such as those housed by cBioPortalData are crude and frequently fail because they tend to be based on receiving appropriately formatted files as input.

By Contrast, MaeUtils repository is founded on the idea that the internal structure and logic of a Mae object is a better starting point for creation of these object.

Specifically, we demonstrate more robust functionality by devoting a series of heuristics to unambiguous identificaiton of a primaryId column, even if there is no exact match for the strings found therein among the other flatfiles in a given directory.

The first functionality intended is automation of the creation of a Mae objects given only a directory containing (appropriate) flat files as input.
