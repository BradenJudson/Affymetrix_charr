# Affymetrix_charr
An R script for processing Affymetrix genotyping array data. Data are formatted into .map and .ped files for future manipulation (e.g., via PLINK, PGDSPIDER, etc.). Locations of polymorphic SNP loci are also visualized on their respective linkage groups.


charrSNP.map and cavecharr_polySNP.ped are complementary files ammenable for filtering via PLINK (Purcell et al., 2007). 
For example, PLINK can then be used to filter by minor allele frequency, deviations from HWE, loci in LDE, etc. 


Data from LGstats.csv are retrieved from the available Arctic charr genome assembly (Christensen et al., 2018 and GenBank accession: GCF_002910315.2). Total lengths are retained.
