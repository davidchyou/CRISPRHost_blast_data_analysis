#PPV statistics of CRISPRHost_blast stratified by kingdom (archaea and bacteria)
#Displayed as a ggplot bar chart (variable gg)

library("ggplot2")
library("dplyr")

lineage_tbl <- read.csv("NCBITax_lookup.csv")
genus_new <- gsub("Candidatus ", "", as.character(lineage_tbl$genus))
lineage_tbl <- cbind(lineage_tbl, data.frame(genus_new = genus_new))

lineage_tbl_agg <- lineage_tbl %>% group_by(genus_new) %>% slice(c(1)) %>% select(genus_new, family, order)
lineage_tbl_agg <- as.data.frame(lineage_tbl_agg)
lineage_tbl_agg$genus_new <- as.character(lineage_tbl_agg$genus_new)
lineage_tbl_agg$family <- gsub("Candidatus ", "", as.character(lineage_tbl_agg$family))
lineage_tbl_agg$order <- gsub("Candidatus ", "", as.character(lineage_tbl_agg$order))

all_genera <- unique(as.character(lineage_tbl_agg$genus_new))

genera_name <- function(str) {		
	ls <- strsplit(str, " ")[[1]]
	ls_common <- intersect(ls, all_genera)
	
	if (length(ls_common) > 0) {	
		return(ls_common[1])
	} else {
		return("NA")
	}
}

to_genera_list <- function(xlst) {
	str <- paste(xlst, collapse = ";")
	toks <- strsplit(str, ";")[[1]]
	toks <- unlist(lapply(toks, FUN=genera_name))
	toks <- unique(toks);
	str <- paste(toks, collapse = ";")		
	return(str)
}

collect_species <- function(x) {
	str <- paste(x, collapse = ";")
	toks <- strsplit(str, ";")[[1]]
	toks <- unique(toks);
	str <- paste(toks, collapse = ";")
	return(str)
}

prepare_match_table <- function(data_in, max_mismatch = 20) {
	expand_row <- function(x) {
		lst_pred_genera <- strsplit(as.character(x$PRED_GENERA)[1], ";")[[1]]
		
		nn <- length(lst_pred_genera)
		
		df <- data.frame(ASSEMBLY = rep(as.character(x$ASSEMBLY)[1], nn), 
					 	 KNOWN_GENERA = rep(as.character(x$KNOWN_GENERA)[1], nn),
		             	 PRED_GENERA = lst_pred_genera,
		             	 E_VALUE = rep(as.numeric(x$E_VALUE)[1], nn),
		             	 BITSCORE = rep(as.numeric(x$BITSCORE)[1], nn),
		             	 PROP_NT_MATCHED_IN_HYBRID = rep(as.numeric(x$PROP_NT_MATCHED_IN_HYBRID)[1], nn),
		             	 N_SPACER_NT_MISMATCHED = rep(as.numeric(x$N_SPACER_NT_MISMATCHED)[1], nn),
		             	 N_HOST_SPACER = rep(as.numeric(x$ N_HOST_SPACER)[1], nn))
		            
		return(df)
	}
	
	all_viral_genome_gcf <- gsub("(\\S+)\\|\\S+\\|\\S+", "\\1", data_in$VIRAL_GENOME)
	all_viral_genome_exp_genera <- gsub("\\S+\\|\\S+\\|(\\S+)", "\\1", data_in$VIRAL_GENOME)
	all_viral_genome_exp_genera[! all_viral_genome_exp_genera %in% all_genera] <- "NA"
	all_species_matched <- data_in$KNOWN_SPECIES
	nmis <- data_in$SPACER_LENGTH - round(data_in$N_NT_MATCHED_IN_HYBRID / 2)

	df_simple <- data.frame(ASSEMBLY = all_viral_genome_gcf, 
                        	KNOWN_GENERA = all_viral_genome_exp_genera, 
                        	PRED_SPECIES = all_species_matched,
                        	E_VALUE = data_in$E_VALUE,
                        	BITSCORE = data_in$BITSCORE,
                        	PROP_NT_MATCHED_IN_HYBRID = data_in$PROP_NT_MATCHED_IN_HYBRID,
                        	N_SPACER_NT_MISMATCHED = abs(nmis),
	                        N_HOST_SPACER = as.numeric(data_in$N_HOST_SPACER))
    
  	df_simple <- df_simple[df_simple$N_SPACER_NT_MISMATCHED <= max_mismatch,]
  	df_simple <- df_simple[df_simple$E_VALUE <= 1e-8,]
  	#df_simple <- df_simple[df_simple$N_HOST_SPACER > 1,]
  	
  	pred_genera <- unlist(lapply(as.character(df_simple$PRED_SPECIES), to_genera_list))
  	df_simple <- cbind(df_simple, data.frame(PRED_GENERA = pred_genera))
  	
  	nr <- dim(df_simple)[1];
  	df_simple_expanded <- do.call(rbind, lapply(1:nr, FUN = function(x){expand_row(df_simple[x,])}))
	
	df_simple_expanded <- df_simple_expanded[as.character(df_simple_expanded$PRED_GENERA) != "NA",]
	df_simple_expanded <- df_simple_expanded[as.character(df_simple_expanded$KNOWN_GENERA) != "NA",]
}

add_distance_data <- function(df, lineage_tbl) {
	all_genera <- unique(c(unique(as.character(df$PRED_GENERA)),unique(as.character(df$KNOWN_GENERA))))
	df_dist <- data.frame(KNOWN_GENERA = rep(all_genera, each = length(all_genera)), 
                      	  PRED_GENERA = rep(all_genera, length(all_genera)))
                          
    df_lineage_known <- data.frame(KNOWN_GENERA = unique(as.character(df_dist$KNOWN_GENERA)))
    df_lineage_known <- merge(df_lineage_known, lineage_tbl, by.x = "KNOWN_GENERA", by.y = "genus_new", all.x=TRUE)
    names(df_lineage_known) <- c("KNOWN_GENERA", "KNOWN_FAMILY", "KNOWN_ORDER")
    
    df_lineage_pred <- data.frame(PRED_GENERA = unique(as.character(df_dist$PRED_GENERA)))
    df_lineage_pred <- merge(df_lineage_pred, lineage_tbl, by.x = "PRED_GENERA", by.y = "genus_new", all.x=TRUE)
    names(df_lineage_pred) <- c("PRED_GENERA", "PRED_FAMILY", "PRED_ORDER")
    
    df_dist <- merge(df_dist, df_lineage_known, by = c("KNOWN_GENERA"), all.x=TRUE)
    df_dist <- merge(df_dist, df_lineage_pred, by = c("PRED_GENERA"), all.x=TRUE)
    
  	df <- merge(df, df_dist, by = c("KNOWN_GENERA", "PRED_GENERA"), all.x=TRUE)
  	
  	b1 <- as.numeric(as.character(df$KNOWN_GENERA) == as.character(df$PRED_GENERA))
  	b2 <- as.numeric(as.character(df$KNOWN_FAMILY) == as.character(df$PRED_FAMILY))
  	b3 <- as.numeric(as.character(df$KNOWN_ORDER) == as.character(df$PRED_ORDER))
  	b1[is.na(b1)] <- 0
  	b2[is.na(b2)] <- 0
  	b3[is.na(b3)] <- 0
  	b4 <- (1 - b1) * b2
  	b5 <- (1 - b1) * (1 - b2) * b3
  	b6 <- (1 - b1) * (1 - b2) * (1 - b3)
  	
  	sort_key <- rep(0, dim(df)[1])
  	sort_key[b4 == 1] <- 8 + df[b4 == 1,]$E_VALUE
  	
  	max_sort_key <- max(sort_key)
  	sort_key[b5 == 1] <- ceiling(max_sort_key) + df[b5 == 1,]$E_VALUE
  	
  	max_sort_key <- max(sort_key)
  	sort_key[b6 == 1] <- ceiling(max_sort_key) + df[b6 == 1,]$E_VALUE
  	
  	dist <- rep(0, dim(df)[1])
  	dist[b4 == 1] <- 8
  	dist[b5 == 1] <- 16
  	dist[b6 == 1] <- 32
  	
  	df <- cbind(df, data.frame(SORT_KEY = sort_key, DIST = dist))
	df <- df %>% group_by(ASSEMBLY) %>% slice(which.min(SORT_KEY))
	df <- as.data.frame(df)
	
	return(df)
}

ppv_stats <- function(data_in, lineage_tbl, max_mismatch, org) {
	df <- prepare_match_table(data_in, max_mismatch)
	df_dist <- add_distance_data(df, lineage_tbl)
	
	n_genome_asgn <- dim(df_dist)[1]
	n_exact <- sum(df_dist$DIST == 0)
	n_family <- sum(df_dist$DIST <= 8)
	n_order <- sum(df_dist$DIST <= 16)
	
	ppv <- n_order / n_genome_asgn
	ppv1 <- n_family / n_genome_asgn
	ppv2 <- n_exact / n_genome_asgn
	
	ppvv <- c(ppv, ppv1, ppv2)
	ppv_type <- c("PPV-Order", "PPV-Family", "PPV-Genus")
	group <- paste("Max mismatch = ", max_mismatch, sep = "")
	lbl <- paste(n_genome_asgn, "genomes assigned")
	
	df_out <- data.frame(ORGANISM = rep(org, 3),
	                     PPV = ppvv,
	                     PPV_TYPE = ppv_type,
	                     MISMATCH_ALLOWED = rep(max_mismatch, 3),
	                     N_GENOME_ASSIGNED = rep(n_genome_asgn,3),
	                     GROUP = rep(group,3),
	                     LBL = c(lbl, "", ""))
}

known_pred_freq_mtrx <- function(df, all_genera) {
  df_count <- df %>% group_by(KNOWN_GENERA, PRED_GENERA) %>% summarize(FREQ = n())
  df_count <- as.data.frame(df_count)
  
  df_count_base <- data.frame(KNOWN_GENERA = rep(all_genera, each = length(all_genera)), 
                              PRED_GENERA = rep(all_genera, length(all_genera)))
  
  df_count_base <- merge(df_count_base, df_count, by = c("KNOWN_GENERA", "PRED_GENERA"), all.x=TRUE)
  df_count_base$FREQ[is.na(df_count_base$FREQ)] <- 0 
  count_mtrx <- matrix(df_count_base$FREQ, length(all_genera), length(all_genera))
  rownames(count_mtrx) <- unique(df_count_base$KNOWN_GENERA)
  colnames(count_mtrx) <- unique(df_count_base$PRED_GENERA)
  count_mtrx <- count_mtrx[all_genera, all_genera]
  
  return(count_mtrx)
}

n_mismatches <- 0:2
orgs <- c("Archaea", "Bacteria")
pdata_paths <- c("validation_data_vhdb_archaea.csv", "validation_data_vhdb_bacteria.csv")
orgs_rep <- rep(orgs, each = length(n_mismatches))
pdata_paths_rep <- rep(pdata_paths, each = length(n_mismatches))
n_mismatches_rep <- rep(n_mismatches, length(orgs))
inds <- 1:(length(n_mismatches_rep))

df_ppv <- do.call(rbind, lapply(inds, FUN=function(x){print(c(n_mismatches_rep[x],orgs_rep[x]));
													  df_in <- read.csv(pdata_paths_rep[x]);  
                                                      ppvdf <- ppv_stats(df_in, lineage_tbl_agg, n_mismatches_rep[x],orgs_rep[x]);
                                                      return(ppvdf)}))
df_ppv$PPV_TYPE <- factor(df_ppv$PPV_TYPE, levels = c("PPV-Order", "PPV-Family", "PPV-Genus"))

gg <- ggplot(df_ppv, aes(x=PPV_TYPE,y=PPV)) + geom_bar(stat="identity",width=0.5) + geom_text(aes(x=2,y=1.1,label=LBL, size=3),show_guide=FALSE)
gg <- gg + facet_grid(GROUP~ORGANISM) + xlab("") + ylab("") + geom_text(aes(label=round(PPV,2)), vjust=1.6, color = "white", size=3)
gg <- gg + theme(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank(),
                 panel.border=element_rect(fill=NA,color="black"),
                 axis.text.x = element_text(angle=90,hjust=1.05,vjust=0.05))
gg <- gg + scale_y_continuous(breaks = seq(0, 1, 0.2), labels = seq(0, 1, 0.2))