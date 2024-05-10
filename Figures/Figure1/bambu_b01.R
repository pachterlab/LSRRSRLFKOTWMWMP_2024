library(bambu)

bamb01 <- "sorted_b01_nanopore_13G_spliced.bam"

fa.file <- "/home/rebekah/mm39.fa"

gtf.file <- "/home/rebekah/mm39.gtf"

bambuAnnotations <- prepareAnnotations(gtf.file)

se <- bambu(reads = bamb01, annotations = bambuAnnotations, genome = fa.file, discovery = FALSE)
writeBambuOutput(se, path = "/home/rebekah/b01_nanopore_13G_bambu")
