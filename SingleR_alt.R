
suppressPackageStartupMessages(library(SummarizedExperiment))

pseudobulk_fil <- "adrian_spec/test_run/p2_outASTRID_clustering_level_11_pseudobulk_matrix.csv"
referens_fil <- "data/ASTRID_SingleR_Reference_20240701.Rds"
exkludera_fil <- "data/excludeSingleR.csv"
kluster_nummer <- 2
facit <- NA_character_ 

pseudobulk <- read.csv(pseudobulk_fil, row.names = 1, check.names = FALSE)
referens <- readRDS(referens_fil)

if (nzchar(exkludera_fil)) {
  exkludera <- read.csv(exkludera_fil, stringsAsFactors = FALSE)
  stopifnot("samples" %in% names(exkludera))
  referens <- referens[, !colnames(referens) %in% exkludera$samples]
}

referens_uttryck <- as.matrix(assay(referens, "counts"))
referens_typer <- as.character(referens$cellTypeGRINT)
stopifnot(ncol(referens_uttryck) > 0,
          length(referens_typer) == ncol(referens_uttryck),
          !anyNA(referens_typer), all(nzchar(referens_typer)),
          !anyDuplicated(colnames(pseudobulk)),
          !anyDuplicated(rownames(referens_uttryck)),
          kluster_nummer >= 1, kluster_nummer <= nrow(pseudobulk),
          kluster_nummer == as.integer(kluster_nummer))

kluster_namn <- rownames(pseudobulk)[kluster_nummer]

gener <- intersect(colnames(pseudobulk), rownames(referens_uttryck))
if (length(gener) < 3) stop("Farre an tre gemensamma gener: kontrollera input.")

kluster_uttryck <- as.numeric(unlist(pseudobulk[kluster_nummer, gener, drop = FALSE]))
referens_uttryck <- referens_uttryck[gener, , drop = FALSE]

if (any(!is.finite(kluster_uttryck)) || any(kluster_uttryck < 0) ||
    any(!is.finite(referens_uttryck)) || any(referens_uttryck < 0) ||
    sum(kluster_uttryck) <= 0 || any(colSums(referens_uttryck) <= 0)) {
  stop("Uttrycket maste vara numeriskt, ändligt, icke-negativt och ha positiv totalsumma.")
}

x <- log10(1 + 1e5 * kluster_uttryck / sum(kluster_uttryck))
y <- log10(1 + 1e5 * sweep(referens_uttryck, 2, colSums(referens_uttryck), "/"))

cat("\nKluster:", kluster_namn,
    "\nGemensamma gener:", length(gener),
    "\nReferensprofiler:", ncol(y), "\n")

resultat <- data.frame(
  referens = colnames(y),
  celltyp = referens_typer,
  RMSD = apply(y, 2, function(profil) sqrt(mean((x - profil)^2))),
  Pearson = apply(y, 2, function(profil) cor(x, profil, method = "pearson")),
  Spearman = apply(y, 2, function(profil) cor(x, profil, method = "spearman")),
  row.names = NULL
)
vinnare <- data.frame()
for (matt in c("RMSD", "Pearson", "Spearman")) {
  varden <- resultat[[matt]]
  giltiga <- which(is.finite(varden))

  cat("\nTre basta referenser enligt", matt, ":\n")
  ordning <- giltiga[order(varden[giltiga], decreasing = matt != "RMSD")]
  print(head(resultat[ordning, ], 3), row.names = FALSE)

  if (length(giltiga) == 0) {
    basta <- data.frame(matt = matt, referens = NA_character_,
                        celltyp = NA_character_, varde = NA_real_)
  } else {
    basta_varde <- if (matt == "RMSD") min(varden[giltiga]) else max(varden[giltiga])
    index <- giltiga[abs(varden[giltiga] - basta_varde) <= 1e-12]
    basta <- data.frame(matt = matt, referens = resultat$referens[index],
                        celltyp = resultat$celltyp[index], varde = varden[index])
  }
  vinnare <- rbind(vinnare, basta)
}

if (!is.na(facit)) {
  vinnare$facit <- facit
  vinnare$samma_celltyp <- vinnare$celltyp == facit
}
cat("\nVinnande referenser (flera rader per matt betyder delad forstaplats):\n")
print(vinnare, row.names = FALSE)