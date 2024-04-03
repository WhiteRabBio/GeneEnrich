get_data_from_KEGG_db <- function(species) {
    PATHID2EXTID <- as.list(get_KEGG_db("KEGGPATHID2EXTID"))
    if (!any(grepl(species, names(PATHID2EXTID)))) {
        stop("input species is not supported by KEGG.db...")
    }
    idx <- grep(species, names(PATHID2EXTID))
    PATHID2EXTID <- PATHID2EXTID[idx]
    PATHID2EXTID.df <- stack(PATHID2EXTID)
    PATHID2EXTID.df <- PATHID2EXTID.df[, c(2,1)]
    PATHID2NAME <- as.list(get_KEGG_db("KEGGPATHID2NAME"))
    PATHID2NAME.df <- data.frame(path=names(PATHID2NAME),
                                 name=unlist(PATHID2NAME))
    build_Anno(PATHID2EXTID.df, PATHID2NAME.df)
}

get_KEGG_db <- function(kw) {
    annoDb <- "KEGG.db"
    suppressMessages(requireNamespace(annoDb))
    eval(parse(text=paste0(annoDb, "::", kw)))
}

build_Anno <- function(path2gene, path2name) {
    if (!exists(".Anno_clusterProfiler_Env", envir = .GlobalEnv)) {
        pos <- 1
        envir <- as.environment(pos) 
        assign(".Anno_clusterProfiler_Env", new.env(), envir = envir)
    }
    Anno_clusterProfiler_Env <- get(".Anno_clusterProfiler_Env", envir= .GlobalEnv)

    if(class(path2gene[[2]]) == 'list') {
        ## to compatible with tibble
        path2gene <- cbind(rep(path2gene[[1]],
                               times = vapply(path2gene[[2]], length, numeric(1))),
                           unlist(path2gene[[2]]))
    }

    path2gene <- as.data.frame(path2gene) 
    path2gene <- path2gene[!is.na(path2gene[,1]), ]
    path2gene <- path2gene[!is.na(path2gene[,2]), ]
    path2gene <- unique(path2gene)
    
    PATHID2EXTID <- split(as.character(path2gene[,2]), as.character(path2gene[,1]))
    EXTID2PATHID <- split(as.character(path2gene[,1]), as.character(path2gene[,2]))
    
    assign("PATHID2EXTID", PATHID2EXTID, envir = Anno_clusterProfiler_Env)
    assign("EXTID2PATHID", EXTID2PATHID, envir = Anno_clusterProfiler_Env)

    if ( missing(path2name) || is.null(path2name) || is.na(path2name)) {
        assign("PATHID2NAME", NULL, envir = Anno_clusterProfiler_Env)
    } else {
        path2name <- as.data.frame(path2name)
        path2name <- path2name[!is.na(path2name[,1]), ]
        path2name <- path2name[!is.na(path2name[,2]), ]
	path2name <- unique(path2name)
	PATH2NAME <- as.character(path2name[,2])
	names(PATH2NAME) <- as.character(path2name[,1]) 
        assign("PATHID2NAME", PATH2NAME, envir = Anno_clusterProfiler_Env)
    }
    return(Anno_clusterProfiler_Env)
}

library(argparser)
argv <- arg_parser('extract kegg database from an existing environment')
argv <- add_argument(argv, "--species", help = "hsa or mmu")
argv <- add_argument(argv, "--output", help="output directory")
argv <- parse_args(argv)

KEGG_DATA <- get_data_from_KEGG_db(argv$species)
qTermID  <- KEGG_DATA$EXTID2PATHID
qExtID2TermID.df <- data.frame(geneid=rep(paste(argv$species, ':', names(qTermID), sep=""),
                                    times=lapply(qTermID, length)), 
                               pathid=paste('path:', unlist(qTermID), sep=""))

PATHID2NAME <- get("PATHID2NAME", envir = KEGG_DATA)
Pathname.df <- data.frame(pathid=paste(argv$species, names(PATHID2NAME), sep=""), 
                          pathname=PATHID2NAME)

write.table(qExtID2TermID.df, paste(argv$output, '/TermID_KEGG_', argv$species, '_df.xls', sep=''), sep='\t', quote=F, row.names=F)
write.table(Pathname.df, paste(argv$output, '/PathwayName_KEGG_', argv$species, '_df.xls', sep=''), sep='\t', quote=F, row.names=F)
