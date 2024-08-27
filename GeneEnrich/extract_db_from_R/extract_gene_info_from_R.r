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
library(org.Mmu.eg.db)
library(argparser)
})


argv <- arg_parser('extract gene info')
argv <- add_argument(argv, "--species", help = "hsa, mmu, etc...")
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

for (from_type in c('UNIPROT', 'SYMBOL')) {
    all_genes <- keys(get(OrgDb), keytype = from_type)
    gene_mapping <- bitr(
        all_genes,
        fromType = from_type,
        toType = "ENTREZID",
        OrgDb = OrgDb
    )
    write.table(gene_mapping, paste(argv$output, '/', from_type, '_', argv$species, '.xls', sep=''), sep='\t', quote=F, row.names=F)
}
