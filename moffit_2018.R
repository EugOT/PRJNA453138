library(SingleCellExperiment)
library(scater)
poa <- readr::read_rds("/home/etretiakov/Dropbox/04_Archives/Projects/Summarise_Hypothalamus/data/resources/moffitt_2018.rds")
assayNames(poa) <- "counts"
counts(poa) <- as(counts(poa), "dgCMatrix")
sum(duplicated(rownames(poa))) # 65
rownames(poa)[duplicated(rownames(poa))]
summary(as.matrix(t(counts(poa[rownames(poa)[duplicated(rownames(poa))],]))))
poa <- poa[unique(rownames(poa)), ]
poa <- logNormCounts(poa)
dim(logcounts(poa))
assays(poa) <- assays(poa)[c("counts", "logcounts")]

readr::write_rds(poa, "/mnt/data/PRJNA453138/sce_moffitt_2018.rds")

poa <- readr::read_rds("/mnt/data/PRJNA453138/sce_moffitt_2018.rds")

poa_gal <- poa[,
               poa$cell_type %in% c(
                   "i8:Gal/Amigo2",
                   "i16:Gal/Th",
                   "i18:Gal/Tac2",
                   "i20:Gal/Moxd1",
                   "i22:Gal/Pmaip1",
                   "e22:Gal/Ucn3",
                   "e24:Gal/Rxfp1"
               )]
agg.n <- summarizeAssayByGroup(poa_gal,
                               statistics = "prop.detected",
                               ids = colData(poa_gal)[,c("cell_type")])
agg.n <- agg.n[rowSums(assay(agg.n)) > 1/7, ]
dim(agg.n)
head(assay(agg.n))
readr::write_csv(tibble::as_tibble(assay(agg.n),
                                   rownames = TRUE),
                 "/mnt/data/PRJNA453138/prop-detect_moffitt_2018_gal.csv")

library(Seurat)
library(SeuratDisk)
poa.gal <-
    as.Seurat(poa_gal,
              counts = "counts",
              data = "logcounts")
poa.gal %<>% RenameAssays(originalexp = "RNA")
poa.gal %<>%
    PercentageFeatureSet(pattern = "^mt-",
                                  col.name = "percent.mt") %>%
    SCTransform(vars.to.regress = "percent.mt",
                verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunTSNE(
        dims = 1:30,
        verbose = FALSE
    )

readr::write_rds(poa.gal, "/mnt/data/PRJNA453138/srt_moffitt_2018_gal.rds")

