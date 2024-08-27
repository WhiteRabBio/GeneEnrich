suppressMessages({
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Gg.eg.db)
library(org.Dm.eg.db)
library(org.Dr.eg.db)
library(org.Bt.eg.db)
library(org.Cf.eg.db)
library(org.Ss.eg.db)
library(AnnotationDbi)
library(argparser)
})


argv <- arg_parser('extract go database from an existing environment')
argv <- add_argument(argv, "--species", help = "hsa or mmu")
argv <- add_argument(argv, "--ont", help = "BP MF or CC")
argv <- add_argument(argv, "--output", help="output directory")
argv <- parse_args(argv)


replace_species_with_orgdb <- function(species_input, database='/Personal/huangwanxiang/GeneEnrich/GeneEnrich/database/convert/') {
    file_path <- paste(database, '/species_orgdb.csv', sep='')
    df <- read.table(file_path, sep=',', header=TRUE)

    if (species_input %in% df$species) {
        orgdb_value <- df$orgdb[df$species == species_input]
    return(orgdb_value)

    } else {
        stop("The provided species is not supported.")
    }
}

OrgDb <- replace_species_with_orgdb(argv$species)
   
keyType <- "ENTREZID"

GO_DATA <- clusterProfiler:::get_GO_data(OrgDb, argv$ont, keyType)
qTermID  <- GO_DATA$EXTID2PATHID
qExtID2TermID.df <- data.frame(geneid=rep(paste(argv$species, ':', names(qTermID), sep=""),
                                    times=lapply(qTermID, length)), 
                               pathid=paste('path:', unlist(qTermID), sep=""))
PATHID2NAME <- get("PATHID2NAME", envir = GO_DATA)
Pathname.df <- data.frame(pathid=paste(names(PATHID2NAME), sep=""), 
                          pathname=PATHID2NAME)

write.table(qExtID2TermID.df, paste(argv$output, '/TermID_GO_', argv$ont, '_', argv$species, '_df.xls', sep=''), sep='\t', quote=F, row.names=F)
write.table(Pathname.df, paste(argv$output, '/PathwayName_GO_', argv$species, '_df.xls', sep=''), sep='\t', quote=F, row.names=F)

