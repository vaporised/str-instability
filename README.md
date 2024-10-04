# str-instability
Code for my Honours project "Genomic Instability of Tandem Repeat Expansion Disorders"

## return_transmissions.R
An algorithm that returns the most likely transmissions of (parent allele ? -> child allele 1) + (parent allele ? -> child allele 2). This algorithm optimises based on minimum repeat length changes from parent to child. However it will always prioritise stable transmissions.

Takes 3 vectors each containing the 2 alleles (repeat lengths) in the trio member for any locus. 

### Example:
```R
mother_alleles <- c(6, 9)
father_alleles <- c(4, 15)
child_alleles <- c(5, 9)

transmissions <- return_transmissions(mother_alleles, father_alleles, child_alleles)
```

Here, the `transmissions` variable will be a list containing `c(9, 9)` as the first element and `c(4, 5)` as the second element. The first element is the maternal transmission, with the 9 repeat allele in the mother being the allele of origin for the 9 repeat allele in the child. The second element is the paternal transmission, with the 4 repeat allele in the father being the allele of origin for the 5 repeat allele in the child. In this example a paternal expansion occurs from 4 to 5 in the child.
