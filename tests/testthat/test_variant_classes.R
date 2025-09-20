library(testthat)
library(GenomicRanges)

test_gr <- GenomicRanges::GRanges("seq1", IRanges::IRanges(123,126), strand = "*")

test_that("allele.check works correctly",{
  expect_null(allele.check("AT", "ref", exp.length = 2))
  expect_match(allele.check("ATN","alt"), "must contain only A,T,G,C characters")
  expect_match(allele.check("A","alt", exp.length = 2), "must be 2 nucleotides")
  expect_match(allele.check("","alt", allow.empty = FALSE), "must be not empty")
  expect_null(allele.check("","alt", allow.empty = TRUE))
  expect_match(allele.check("A", "ref", min.length = 3), "must be at least 3 nucleotides")
})

test_that("SNV objects are created and validated correctly", {
  snv <- SNV("rs123",test_gr,"A","G","TP53","protein_coding","missense", HGVSp = "p.R381L")
  expect_s4_class(snv, "SNV")
  expect_equal(refAllele(snv), "A")
  expect_identical(snv@geneType, "protein_coding")
  expect_identical(snv@geneSymbol, "TP53")
  expect_identical(snv@altAllele, "G")
  expect_identical(snv@geneType, "protein_coding")
  expect_identical(snv@geneSymbol, "TP53")
  expect_error(SNV("rs123",test_gr,"AT","G","TP53","protein_coding", "missense"),
               "SNV refAllele must be 1 nucleotide")
  expect_error(SNV( "rs123",test_gr,"A", "X","TP53", "protein_coding","missense"),
    "must contain only A,T,G,C characters")
})

test_that("DNV", {
  dnv <- DNV("dnv1",test_gr,"AT","GC","BRCA1","protein_coding","missense")
  expect_s4_class(dnv, "DNV")
  expect_equal(altAllele(dnv), "GC")
  expect_error(DNV("rs1", test_gr , "A", "GC", "BRCA1", "protein_coding", "missense"),
               "DNV refAllele must be 2 nucleotides")
  expect_error(DNV("rs1", test_gr , "AT", "X", "BRCA1", "protein_coding", "missense"),
               "DNV altAllele must contain only A,T,G,C characters")
})

test_that("ONV objects are created and validated correctly", {
  onv <- ONV("onv1",test_gr,"ATG","TAC","KRAS","protein-coding","missense")
  expect_s4_class(onv, "ONV")
  expect_error(ONV("onv1",test_gr,"AT","TA","KRAS","protein-coding","missense"),
               "ONV refAllele must be at least 3 nucleotides")
})

test_that("Insertion", {
  ins <- Insertion("ins1", test_gr , "-", "A", "MYC", "protein_coding", "frameshift")
  expect_s4_class(ins, "Insertion")
  expect_equal(refAllele(ins), "-")
  expect_s4_class(Insertion("ins2", test_gr,"","ATG","EGFR","protein-coding","frameshift"),
                  "Insertion")
  expect_error(Insertion("ins1",test_gr,ref = "A", "ATG", "EGFR","protein-coding","frameshift"),
               "Insertion refAllele must be a '-' or empty")
})

test_that("Deletion", {
  del <- Deletion("del1", test_gr , "T", "-", "BRCA2", "protein_coding", "frameshift")
  expect_s4_class(del, "Deletion")
  expect_equal(altAllele(del), "-")
  expect_error(Deletion("del1",test_gr,"ATG","A","APC", "protein-coding","frameshift"),
               "Deletion altAllele must be a '-' or empty")
})


test_that("Accessors work", {
  snv <- SNV("rs123",test_gr,"A","G", "TP53", "protein-coding","missense","Pathogenic")
  expect_equal(variantID(snv), "rs123")
  expect_equal(position(snv), test_gr)
  expect_equal(geneSymbol(snv), "TP53")
  expect_equal(clinicalRel(snv), "Pathogenic")
})



test_that("Setter methods work with validation", {
  snv <- SNV("rs123",test_gr,"A","G","TP53","protein-coding","missense")
  altAllele(snv) <- "T"
  expect_equal(altAllele(snv), "T")
  expect_error(altAllele(snv) <- "GT", "SNV altAllele must be 1 nucleotide")
  clinicalRel(snv) <- "Benign"
  expect_equal(clinicalRel(snv), "Benign")
})




test_that("Show method displays correctly", {
  snv <- SNV("rs123",test_gr, "A","G","TP53","protein-coding","missense", HGVSp ="p.R381L")
  output <- capture.output(show(snv))
  expect_true(any(grepl("Variant ID: rs123", output)))
  expect_true(any(grepl("Gene: TP53 protein-coding", output)))
  expect_true(any(grepl("Ref -> Alt: A -> G", output)))
  expect_true(any(grepl("HGVS protein: p.R381L", output)))
})



test_that("aamutation shows correct position", {
  snv1 <- SNV("rs123",test_gr, "A","G", "TP53","protein-coding","missense", HGVSp = "p.A381L")
  expect_equal(aamutation(snv1), 381)
  snv2 <- SNV(
    variantID = "rs456",
    pos = test_gr,
    ref = "A",
    alt = "G",
    geneSymb = "TP53",
    geneType = "protein-coding",
    consequence = "silent",
    HGVSp = "p.=")
  expect_equal(aamutation(snv2), "-")
  snv2 <- SNV("rs123", test_gr,"A","T","TP53","protein_coding","silent")
  expect_equal(aamutation(snv2), "-")
})

