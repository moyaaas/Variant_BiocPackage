#' @import methods
#' @import GenomicRanges
#' @importFrom GenomicRanges seqnames start end
#' @import stringr


#' @title Variant Virtual Class
#' @description virtual base class for all mutation types
#' @aliases Variant-class
#' @slot variantID Character. Unique ID for the variant
#' @slot position GRanges. Genomics Coordinates
#' @slot refAllele Character. Reference Allele
#' @slot altAllele Character. Alternative Allele
#' @slot geneSymbol Character. Symbol of the affected gene
#' @slot geneType Character. Type of the mutated gene
#' @slot consequence Character. Type of consequence - effect
#' @slot clinicalRel Character. Clinical Relevance information
#' @slot EnsemblID Character. Ensembl transcript ID
#' @slot HGVSp Character. Amino Acid change in HGVS format
#' @export
setClass("Variant", slots = list(variantID   = "character",
                                 position    = "GRanges",
                                 refAllele   = "character",
                                 altAllele   = "character",
                                 geneSymbol  = "character",
                                 geneType    = "character",
                                 consequence = "character",
                                 clinicalRel = "character",
                                 EnsemblID  = "character",
                                 HGVSp       = "character"),
         contains = "VIRTUAL",
         prototype = list(clinicalRel = NA_character_,
                          EnsemblID  = NA_character_,
                          HGVSp       = NA_character_))




#' SNV class
#' @description Class for Single Nucleotide Variants
#' This class inherits from the \code{\link{Variant-class}} virtual class
#' @export
#' @examples
#' library(GenomicRanges)
#' gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 100))
#' snv <- new("SNV", variantID="rs1", position=gr, refAllele="A", altAllele="T",
#'            geneSymbol="MYC", geneType="protein_coding", consequence="missense",
#'            clinicalRel=NA_character_, EnsemblID=NA_character_, HGVSp=NA_character_)
#' snv
setClass("SNV",
         contains = "Variant")




#' DNV class
#' @description Class for Double Nucleotide Variants
#' This class inherits from the \code{\link{Variant-class}} virtual class
#' @export
setClass("DNV",
         contains = "Variant")




#' ONV class
#' @description Class for Oligo-Nucleotide Variants
#' This class inherits from the \code{\link{Variant-class}} virtual class
#' @export
setClass("ONV",
         contains = "Variant")




#' Insertion class
#' @description Class for insertion mutations
#' This class inherits from the \code{\link{Variant-class}} virtual class
#' @export
setClass("Insertion",
         contains = "Variant")




#' Deletion class
#' @description Class for deletion mutations
#' This class inherits from the \code{\link{Variant-class}} virtual class
#' @export
setClass("Deletion",
         contains = "Variant")





#' Show Variant object
#' Prints summary of objects of Variant class and their subclasses
#' @param object an object of the class Variant or any class that inherits from it
#' @return None
#' It prints a formatted summary of the object to the console
#' @export
#' @name show
#' @aliases show,Variant-method
#' @docType methods
#' @rdname show-methods
setMethod("show",
          "Variant",
          function(object) {
            pos_range <- object@position
            pos_str <- paste0(seqnames(pos_range), ":", 
                              start(pos_range), "-", end(pos_range))
            cat("Variant ID:", object@variantID, "\n")
            cat("Gene:", object@geneSymbol, object@geneType, "\n")
            cat("Position:", pos_str, "\n")
            cat("Ref -> Alt:", object@refAllele, "->", object@altAllele, "\n")
            cat("Consequence:", object@consequence, "\n")
            if (!is.na(object@clinicalRel)) cat("Clinical Relevance:", object@clinicalRel, "\n")
            if (!is.na(object@EnsemblID)) cat("Ensembl ID:", object@EnsemblID, "\n")
            if (!is.na(object@HGVSp)) cat("HGVS protein:", object@HGVSp, "\n")
          })




#' Validate nucleotide allele
#' It checks if A given allele sequence is valid according to several criteria
#' @param allele A character string representing the allele
#' @param allele.type A character string to describe the allele (altAllele, refAllele)
#' @param exp.length Optional integer, the expected required length
#' @param min.length Optional integer, the minimum required length of the allele
#' @param allow.empty Logical, if false, if the allele is empty the function returns an error
#' @return \code{NULL} if the allele is valid, otherwise, a character string with an error
#' @keywords internal
allele.check <- function(allele, allele.type, exp.length = NULL,
                         min.length = NULL, allow.empty = FALSE) {
  msg <- character()
  if (!grepl("^[ATGC]*$", allele)) {
    msg <- paste(allele.type, "must contain only A,T,G,C characters")
    return(msg)
  }
  if (!is.null(exp.length)) {
    if (nchar(allele) != exp.length) {
      msg <- paste(allele.type, "must be", exp.length, "nucleotides")
      return(msg)
    }
  }
  if (!is.null(min.length)) {
    if (nchar(allele) < min.length) {
      msg <- paste(allele.type, "must be at least", min.length, "nucleotides")
      return(msg)
    }
  }
  if (!allow.empty && nchar(allele) == 0) {
    msg <- paste(allele.type, "must be not empty")
    return(msg)
  }
  return(NULL)
}




#' @title Validation method for a DNV object
#' @description It checks if the object is a valid \code{DNV}, it checks both ref and alt allele (exactly 2 nucleotides)
#' @name DNV_validation
#' @param object a \code{DNV} object
#' @return \code{TRUE} if the object is valid, otherwise a message with the error
#' @keywords internal
setValidity("DNV", function(object) {
  ref.msg <- allele.check(object@refAllele, "DNV refAllele", exp.length = 2)
  alt.msg <- allele.check(object@altAllele, "DNV altAllele", exp.length = 2)
  msgs <- c(ref.msg, alt.msg)
  msgs <- msgs[!vapply(msgs, is.null, logical(1))]
  if (length(msgs) == 0)
    TRUE
  else
    msgs
})




#' @title Validation method for a SNV object
#' @description It checks if the object is a valid \code{SNV}, it checks both ref and alt allele (exactly 1 nucleotide)
#' @name SNV_validation
#' @param object a \code{SNV} object
#' @return \code{TRUE} if the object is valid, otherwise a message with the error
#' @keywords internal
setValidity("SNV",
            function(object) {
              msgs <- character()
              ref.msg <- allele.check(object@refAllele, "SNV refAllele", exp.length = 1)
              alt.msg <- allele.check(object@altAllele, "SNV altAllele", exp.length = 1)
              msgs <- c(ref.msg, alt.msg)
              msgs <- msgs[!vapply(msgs, is.null, logical(1))]
              if (length(msgs) == 0)
                return(TRUE)
              else msgs
            })




#' @title Validation method for a ONV object
#' @description It checks if the object is a valid \code{ONV}, it checks both ref and alt allele (at least 3 nucleotides)
#' @name ONV_validation
#' @param object a \code{ONV} object
#' @return \code{TRUE} if the object is valid, otherwise a message with the error
#' @keywords internal
setValidity("ONV", function(object) {
  msgs <- character()
  msgs <- character()
  ref.msg <- allele.check(object@refAllele, "ONV refAllele", min.length = 3)
  alt.msg <- allele.check(object@altAllele, "ONV altAllele", min.length = 3)
  msgs <- c(ref.msg, alt.msg)
  msgs <- msgs[!vapply(msgs, is.null, logical(1))]
  if (length(msgs)==0)
    TRUE
  else
    msgs
})




#' @title Validation method for Insertion object
#' @description It checks if the object is a valid \code{Insertion}, it checks both ref and alt allele
#' @name Insertion_validation
#' @param object an \code{Insertion} object
#' @return \code{TRUE} if the object is valid, otherwise a message with the error
#' @keywords internal
setValidity("Insertion",
            function(object) {
              msgs <- character()
              if (!(object@refAllele %in% c("-",""))) {
                msgs <- c(msgs, "Insertion refAllele must be a '-' or empty")
              }
              alt.msg<- allele.check(object@altAllele, "Insertion altAllele",
                                     min.length = 1)
              if (!is.null(alt.msg))
                msgs <- c(msgs, alt.msg)
              if (length(msgs) == 0)
                return(TRUE)
              else msgs
            })




#' @title Validation method for Deletion object
#' @description It checks if the object is a valid \code{Deletion}, it checks both ref and alt allele
#' @name Deletion_validation
#' @param object an \code{Deletion} object
#' @return \code{TRUE} if the object is valid, otherwise a message with the error
#' @keywords internal
setValidity("Deletion",
            function(object) {
              msgs <- character()
              ref.msg<- allele.check(object@refAllele, "Deletion refAllele",
                                     min.length = 1)
              if (!is.null(ref.msg))
                msgs <- c(msgs, ref.msg)
              if (!(object@altAllele %in% c("-",""))) {
                msgs <- c(msgs, "Deletion altAllele must be a '-' or empty")
              }
              if (length(msgs) == 0)
                return(TRUE)
              else msgs
            })






#' @title SNV constructor
#' @description Creates an object of the class SNV, Single Nucleotide Variant
#' @param variantID Character, unique identifier for the variant
#' @param pos A Granges object indicating the genomic coordinates
#' @param ref Character, reference allele at the specified position
#' @param alt Character, alternative allele
#' @param geneSymb Character, gene symbol of the mutated gene
#' @param geneType Character, type of gene (e.g, protein-coding, lncRNA, etc.)
#' @param consequence Character, Functional consequence (e.g. missense, nonsense, silent, frameshift, etc.)
#' @param clinicalRel Character, clinical relevance (if available, Clinvar of COSMIC info). Optional
#' @param EnsID Character, Ensemble ID of the transcript affected by the variant. Optional
#' @param HGVSp HGVS protein level annotation (e.g. p.R381L). Optional
#' @return An object of the class SNV
#' @export
#' @examples
#' #library(GenomicRanges)
#' #gr <- GenomicRanges::GRanges(seqnames = 'chr1', ranges = IRanges::IRanges(1234, width = 1))
#' #SNV('rs123', gr, 'A','T','TP53', 'protein_coding', 'missense')
SNV <- function(variantID, pos, ref, alt, geneSymb, geneType,
                consequence, clinicalRel = NA_character_,
                EnsID = NA_character_, HGVSp = NA_character_ ) {
  object <- new("SNV",
                variantID = variantID,
                position  = pos,
                refAllele = ref,
                altAllele = alt,
                geneSymbol = geneSymb,
                geneType = geneType,
                consequence = consequence,
                clinicalRel = clinicalRel,
                EnsemblID = EnsID,
                HGVSp = HGVSp
  )
  validObject(object)
  return(object)
}




#' @title DNV constructor
#' @description Creates an object of the class DNV, Double Nucleotide Variant
#' @param variantID Character, unique identifier for the variant
#' @param pos A Granges object indicating the genomic coordinates
#' @param ref Character, reference allele at the specified position
#' @param alt Character, alternative allele
#' @param geneSymb Character, gene symbol of the mutated gene
#' @param geneType Character, type of gene (e.g, protein-coding, lncRNA, etc.)
#' @param consequence Character, Functional consequence (e.g. missense, nonsense, silent, frameshift, etc.)
#' @param clinicalRel Character, clinical relevance (if available, Clinvar of COSMIC info). Optional
#' @param EnsID Character, Ensemble ID of the transcript affected by the variant. Optional
#' @param HGVSp HGVS protein level annotation (e.g. p.R381L). Optional
#' @return An object of the class DNV
#' @export
#' @examples
#' #library(GenomicRanges)
#' #gr <- GenomicRanges::GRanges(seqnames = 'chr1', ranges = IRanges::IRanges(1234, width = 2))
#' #DNV("rs456", gr, "AA","TT","BRCA1","protein_coding","nonsense")
DNV <- function(variantID, pos, ref, alt, geneSymb, geneType,
                consequence, clinicalRel = NA_character_,
                EnsID = NA_character_, HGVSp = NA_character_ ) {
  object <- new("DNV",
                variantID = variantID,
                position  = pos,
                refAllele = ref,
                altAllele = alt,
                geneSymbol = geneSymb,
                geneType = geneType,
                consequence = consequence,
                clinicalRel = clinicalRel,
                EnsemblID = EnsID,
                HGVSp = HGVSp
  )
  validObject(object)
  return(object)
}




#' @title ONV constructor
#' @description Creates an object of the class ONV, Oligo-Nucleotide Variant
#' @param variantID Character, unique identifier for the variant
#' @param pos A Granges object indicating the genomic coordinates
#' @param ref Character, reference allele at the specified position
#' @param alt Character, alternative allele
#' @param geneSymb Character, gene symbol of the mutated gene
#' @param geneType Character, type of gene (e.g, protein-coding, lncRNA, etc.)
#' @param consequence Character, Functional consequence (e.g. missense, nonsense, silent, frameshift, etc.)
#' @param clinicalRel Character, clinical relevance (if available, Clinvar of COSMIC info). Optional
#' @param EnsID Character, Ensemble ID of the transcript affected by the variant. Optional
#' @param HGVSp HGVS protein level annotation (e.g. p.R381L). Optional
#' @return An object of the class ONV
#' @export
#' @examples
#' #library(GenomicRanges)
#' #gr <- GenomicRanges::GRanges(seqnames = 'chr1', ranges = IRanges::IRanges(1234, width = 3))
#' #ONV("rs789", gr, "TAC", "GGA", "EGFR", "protein-coding", "frameshift")
ONV <- function(variantID, pos, ref, alt, geneSymb, geneType,
                consequence, clinicalRel = NA_character_,
                EnsID = NA_character_, HGVSp = NA_character_ ) {
  object <- new("ONV",
                variantID = variantID,
                position  = pos,
                refAllele = ref,
                altAllele = alt,
                geneSymbol = geneSymb,
                geneType = geneType,
                consequence = consequence,
                clinicalRel = clinicalRel,
                EnsemblID = EnsID,
                HGVSp = HGVSp
  )
  validObject(object)
  return(object)
}





#' @title Insertion constructor
#' @description Creates an object of the class Insertion
#' @param variantID Character, unique identifier for the variant
#' @param pos A Granges object indicating the genomic coordinates
#' @param ref Character, reference allele at the specified position
#' @param alt Character, alternative allele
#' @param geneSymb Character, gene symbol of the mutated gene
#' @param geneType Character, type of gene (e.g, protein-coding, lncRNA, etc.)
#' @param consequence Character, Functional consequence (e.g. missense, nonsense, silent, frameshift, etc.)
#' @param clinicalRel Character, clinical relevance (if available, Clinvar of COSMIC info). Optional
#' @param EnsID Character, Ensemble ID of the transcript affected by the variant. Optional
#' @param HGVSp HGVS protein level annotation (e.g. p.R381L). Optional
#' @return An object of the class Insertion
#' @export
#' @examples
#' #library(GenomicRanges)
#' #gr <- GenomicRanges::GRanges(seqnames = 'chr1', ranges = IRanges::IRanges(1234, width = 0))
#' #Insertion("rsIns", gr, "-", "AG", "ALK", "protein-coding", "inframe_insertion")
Insertion <- function(variantID, pos, ref, alt, geneSymb, geneType,
                      consequence, clinicalRel = NA_character_,
                      EnsID = NA_character_, HGVSp = NA_character_ ) {
  object <- new("Insertion",
                variantID = variantID,
                position  = pos,
                refAllele = ref,
                altAllele = alt,
                geneSymbol = geneSymb,
                geneType = geneType,
                consequence = consequence,
                clinicalRel = clinicalRel,
                EnsemblID = EnsID,
                HGVSp = HGVSp
  )
  validObject(object)
  return(object)
}




#' @title Deletion constructor
#' @description Creates an object of the class Deletion
#' @param variantID Character, unique identifier for the variant
#' @param pos A Granges object indicating the genomic coordinates
#' @param ref Character, reference allele at the specified position
#' @param alt Character, alternative allele
#' @param geneSymb Character, gene symbol of the mutated gene
#' @param geneType Character, type of gene (e.g, protein-coding, lncRNA, etc.)
#' @param consequence Character, Functional consequence (e.g. missense, nonsense, silent, frameshift, etc.)
#' @param clinicalRel Character, clinical relevance (if available, Clinvar of COSMIC info). Optional
#' @param EnsID Character, Ensemble ID of the transcript affected by the variant. Optional
#' @param HGVSp HGVS protein level annotation (e.g. p.R381L). Optional
#' @return An object of the class Deletion
#' @export
#' @examples
#' #library(GenomicRanges)
#' #gr <- GenomicRanges::GRanges(seqnames = 'chr1', ranges = IRanges::IRanges(1234, width = 2))
#' #Deletion("rsDel", gr, "CT", "-", "MYC", "protein-coding", "frameshift")
Deletion <- function(variantID, pos, ref, alt, geneSymb, geneType,
                     consequence, clinicalRel = NA_character_,
                     EnsID = NA_character_, HGVSp = NA_character_ ) {
  object <- new("Deletion",
                variantID = variantID,
                position  = pos,
                refAllele = ref,
                altAllele = alt,
                geneSymbol = geneSymb,
                geneType = geneType,
                consequence = consequence,
                clinicalRel = clinicalRel,
                EnsemblID = EnsID,
                HGVSp = HGVSp
  )
  validObject(object)
  return(object)

}








#' @title altAllele accessor
#' @description Get the alternative allele from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # altAllele(variant)
#' @export
setGeneric("altAllele", function(object) standardGeneric("altAllele"))
#' @describeIn altAllele Method for Variant class
#' @export
setMethod(f = "altAllele",
          signature = "Variant",
          definition = function(object) {
            object@altAllele
          })




#' @title refAllele accessor
#' @description Get the reference allele from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # refAllelle(variant)
#' @export
setGeneric("refAllele", function(object) standardGeneric("refAllele"))
#' @describeIn refAllele Method for Variant class
#' @export
setMethod(f = "refAllele",
          signature = "Variant",
          definition = function(object) {
            object@refAllele
          })





#' @title variantID accessor
#' @description Get the variantID from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # variantID(variant)
#' @export
setGeneric("variantID", function(object) standardGeneric("variantID"))
#' @describeIn variantID Method for Variant class
#' @export
setMethod(f = "variantID",
          signature = "Variant",
          definition = function(object) {
            object@variantID
          })




#' @title position accessor
#' @description Get the position from a Variant Object
#' @param object a Variant object
#' @return A GRanges object
#' @examples
#' # position(variant)
#' @export
setGeneric("position", function(object) standardGeneric("position"))
#' @describeIn position Method for Variant class
#' @export
setMethod(f = "position",
          signature = "Variant",
          definition = function(object) {
            object@position
          })




#' @title geneSymbol accessor
#' @description Get the geneSymbol from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # geneSymbol(variant)
#' @export
setGeneric("geneSymbol", function(object) standardGeneric("geneSymbol"))
#' @describeIn geneSymbol Method for Variant class
#' @export
setMethod(f = "geneSymbol",
          signature = "Variant",
          definition = function(object) {
            object@geneSymbol
          })




#' @title geneType accessor
#' @description Get the geneType from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # geneType(variant)
#' @export
setGeneric("geneType", function(object) standardGeneric("geneType"))
#' @describeIn geneType Method for Variant class
#' @export
setMethod(f = "geneType",
          signature = "Variant",
          definition = function(object) {
            object@geneType
          })




#' @title consequence accessor
#' @description Get the consequence from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # consequence(variant)
#' @export
setGeneric("consequence", function(object) standardGeneric("consequence"))
#' @describeIn consequence Method for Variant class
#' @export
setMethod(f = "consequence",
          signature = "Variant",
          definition = function(object) {
            object@consequence
          })




#' @title clinicalRel accessor
#' @description Get the clinicalRel from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # clinicalRel(variant)
#' @export
setGeneric("clinicalRel", function(object) standardGeneric("clinicalRel"))
#' @describeIn clinicalRel Method for Variant class
#' @export
setMethod(f = "clinicalRel",
          signature = "Variant",
          definition = function(object) {
            object@clinicalRel
          })




#' @title EnsemblID accessor
#' @description Get the EnsemblID from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # EnsemblID(variant)
#' @export
setGeneric("EnsemblID", function(object) standardGeneric("EnsemblID"))
#' @describeIn EnsemblID Method for Variant class
#' @export
setMethod(f = "EnsemblID",
          signature = "Variant",
          definition = function(object) {
            object@EnsemblID
          })




#' @title HGVSp accessor
#' @description Get the HGVSp from a Variant Object
#' @param object a Variant object
#' @return A character string
#' @examples
#' # HGVSp(variant)
#' @export
setGeneric("HGVSp", function(object) standardGeneric("HGVSp"))
#' @describeIn HGVSp Method for Variant class
#' @export
setMethod(f = "HGVSp",
          signature = "Variant",
          definition = function(object) {
            object@HGVSp
          })













#' @title altAllele slot setter
#' @description Set a new alternative allele
#' @name altAllele.setter
#' @param object A variant object
#' @param value A character string representing the new alternative allele
#' @return The new Variant object
#' @aliases altAllele<-
#' @export
setGeneric("altAllele<-", function(object, value) {
  standardGeneric("altAllele<-")
})
#' @describeIn altAllele.setter Setter method for the Variant
#' @export
setMethod(f = "altAllele<-",
          signature = "Variant",
          definition = function(object, value) {
            object@altAllele <- value
            validObject(object)
            return(object)
          })




#' @title clinicalRel slot setter
#' @description Set a new clinical Relevance info
#' @name clinicalRel.setter
#' @param object A variant object
#' @param value A character string representing the new clinical Relevance
#' @return The new Variant object
#' @aliases clinicalRel<-
#' @export
setGeneric("clinicalRel<-", function(object, value) {
  standardGeneric("clinicalRel<-")
})
#' @describeIn clinicalRel.setter Setter method for the Variant
#' @export
setMethod(f = "clinicalRel<-",
          signature = "Variant",
          definition = function(object, value) {
            object@clinicalRel <- value
            validObject(object)
            return(object)
          })





#' @title EnsemblID slot setter
#' @description Set a new EnsemblID
#' @name EnsemblID.setter
#' @param object A variant object
#' @param value A character string representing the new EnsemblID
#' @return The new Variant object
#' @aliases EnsemblID<-
#' @export
setGeneric("EnsemblID<-", function(object, value) {
  standardGeneric("EnsemblID<-")
})
#' @describeIn EnsemblID.setter Setter method for the Variant
#' @export
setMethod(f = "EnsemblID<-",
          signature = "Variant",
          definition = function(object, value) {
            object@EnsemblID <- value
            validObject(object)
            return(object)
          })





#' @title HGVSp slot setter
#' @description Set a new HGVSp
#' @name HGVSp.setter
#' @param object A variant object
#' @param value A character string representing the new HGVSp
#' @return The new Variant object
#' @aliases HGVSp<-
#' @export
setGeneric("HGVSp<-", function(object, value) {
  standardGeneric("HGVSp<-")
})
#' @describeIn HGVSp.setter Setter method for the Variant
#' @export
setMethod(f = "HGVSp<-",
          signature = "Variant",
          definition = function(object, value) {
            object@HGVSp <- value
            validObject(object)
            return(object)
          })








#' @title aamutation function
#' @description It extracts the AminoAcid position from the variant. If the mutation is a silent mutation or if the info is missing, it returns a "-"
#' @param object A Variant Object
#' @return An integer position or '-'
#' @export
setGeneric("aamutation", function(object) {
  standardGeneric("aamutation")
})
#' @describeIn aamutation Method for Variant class
#' @export
setMethod(f = "aamutation",
          signature = "Variant",
          definition = function(object) {
            hgvsp <- object@HGVSp
            if (is.na(hgvsp) || hgvsp == "p.=")
              return("-")
            match <- str_extract(hgvsp, "[0-9]+")
            if (is.na(match))
              return("-")
            else
              return(as.integer(match))
          })
