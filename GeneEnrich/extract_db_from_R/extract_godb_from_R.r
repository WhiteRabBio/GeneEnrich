suppressMessages({
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(argparser)
})


argv <- arg_parser('extract go database from an existing environment')
argv <- add_argument(argv, "--species", help = "hsa or mmu")
argv <- add_argument(argv, "--ont", help = "BP MF or CC")
argv <- add_argument(argv, "--output", help="output directory")
argv <- parse_args(argv)

if (argv$species == 'hsa'){
   OrgDb <- "org.Hs.eg.db"
}else if (argv$species == 'mmu') {
   OrgDb <- "org.Mm.eg.db"
} else {
   stop("species is not supported...")
}
   
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

