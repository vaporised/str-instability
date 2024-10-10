# str-instability
Code and data for my Honours project "Genomic Instability of Tandem Repeat Expansion Disorders"


## str_disease_trios.csv
This file contains the extracted genotypes for all 71 pathogenic loci across all related 1000 Genomes Project (1KG) individuals. A '0' in the `FatherID` and `MotherID` columns means the parent's data is not in 1KG. A '1' in `Sex` is male, '2' for female.

## return_transmissions.R
An algorithm that returns the most likely inherited parental STR alleles for a trio. This is returned as a list of two vectors, where each vector represents the transmission of an STR allele to the child. The first vector is the maternal transmission and the second is the paternal transmission. It takes 3 vectors as input that represent the genotypes of the trio members (in repeat units), in the order of mother, father, and child. This supports the different inheritance patterns of loci on the X chromosome. This algorithm optimises based on a least-euclidean distance method to minimise the repeat magnitude changes from parent to child. However, it will always prioritise stable transmissions before a lower average repeat magnitude change.

Takes 3 vectors each containing the 2 alleles (repeat lengths) in the trio member for any locus. 

### Example:
```R
mother_alleles <- c(6, 9)
father_alleles <- c(4, 15)
child_alleles <- c(5, 9)

transmissions <- return_transmissions(mother_alleles, father_alleles, child_alleles)
```

Here, the `transmissions` variable will be a list containing `c(9, 9)` as the first element and `c(4, 5)` as the second element. The first element is the maternal transmission, with the 9 repeat allele in the mother being the allele of origin for the 9 repeat allele in the child. The second element is the paternal transmission, with the 4 repeat allele in the father being the allele of origin for the 5 repeat allele in the child. In this example a paternal expansion occurs from 4 to 5 in the child.
