# Function for excluding very frequently generated scars in lineage tracing data based on their occurrence in many bulk datasets.
# Written by Bastiaan Spanjaard

suppressPackageStartupMessages(library(ggplot2))
require(igraph)
suppressPackageStartupMessages(library(reshape2))
require(stringdist)
require(plyr)
require(data.table)
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(plyr))
#suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(pheatmap))


# Function for adding probabilities for scars in the data by comparing to the frequency the scar is found in in bulk sequencing data (from larvae).
# If probabilities of certain genes are not available, they will automatically be assigned a low probability to avoid being filtered out.

compare_scars <- function(input_dat_name, dat_name, gene_names, scar_probabilities_in, seur_dat, primer_set, min.presence, ...){

    all.scars.g1 <- input_dat_name
    all.scars.g1 <- all.scars.g1[all.scars.g1$Barcode %in% seur_dat$Barcode,]
    
    print(paste0('Starting with ', nrow(all.scars.g1),' entries in scar file.'))
    
    for(i in 1:length(gene_names)){
    
        gene_it <- gene_names[i]
        
        print(paste0('Comparing gene ', gene_it))
        
        Z1.scars <- f_scars_in_trnscrptm[f_scars_in_trnscrptm$Gene == gene_it,]
        
        
        if(nrow(Z1.scars) == 0){
        
            print(paste0(gene_it, " not found! Moving to next gene."))
            next
        }
        
        # Count presence of scars ####
        unique.scars <- data.frame(table(Z1.scars$Sequence))
        colnames(unique.scars) <- c("Sequence", "Freq.Z1")

        # This comes after merging scars from several datasets - maybe it is not needed here
        unique.scars[is.na(unique.scars)] <- 0
        unique.scars$Presence <- apply(as.data.frame(unique.scars[, -1]), 1,
                                       function(x) sum(x >= min.presence))

        all.CIGARs <- unique(Z1.scars[, c("Sequence", "CIGAR")])
        all.CIGARs <- all.CIGARs[!duplicated(all.CIGARs$Sequence), ]
        unique.scars <- merge(unique.scars, all.CIGARs)
        unique.scars$Sequence <- as.character(unique.scars$Sequence)
        
    
        # Compare presence with probabilities ####
        
         if(primer_set == "1"){                              
                                       
            if(gene_it == "actb1"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 5, 75)

            }else if(gene_it == "actb2"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 73)

            }else if(gene_it == "rpl39"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 14, 75)

            }else if(gene_it == "rpl18a"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 6, 75)
            }else if(gene_it == "cfl1"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 75)

                unique.scars$p <- 0.0001
                unique.scars$Embryos <- 0
                unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
                unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
                unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

                col_order <- c('Sequence.short','Sequence','Freq.Z1','Presence','CIGAR','p','Embryos','Scar')
                unique.scars_2 <- unique.scars[, col_order]
            } else if(gene_it == "dsRedRecCas"){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 116)

                unique.scars$p <- 0.0001
                unique.scars$Embryos <- 0
                unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
                unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
                unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

                col_order <- c('Sequence.short','Sequence','Freq.Z1','Presence','CIGAR','p','Embryos','Scar')
                unique.scars_2 <- unique.scars[, col_order]
            }


        }else if(primer_set == "2"){
            
            if(gene_it %in% c('actb1','actb2','rpl39')){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 75)

            }else if(gene_it %in% c('cfl1','cirbpb','ube2e1')){
                unique.scars$Sequence.short <- 
                substr(unique.scars$Sequence, 1, 75)

                unique.scars$p <- 0.0001
                unique.scars$Embryos <- 0
                unique.scars$p[unique.scars$Sequence.short == wildtype.clip] <- max(unique.scars$p, na.rm = T) + 0.1
                unique.scars <- unique.scars[order(-unique.scars$Presence, -unique.scars$p), ]
                unique.scars$Scar <- paste(0:(nrow(unique.scars)-1), unique.scars$CIGAR, sep = ":")

                col_order <- c('Sequence.short','Sequence','Freq.Z1','Presence','CIGAR','p','Embryos','Scar')
                unique.scars_2 <- unique.scars[, col_order]
        
            }
            
        }else{
            print("Please supply a valid name for the primer set used! At the moment options are Nora or Nina.")
            break   
        }
            
            
        if(!gene_it %in% c('cfl1','cirbpb','ube2e1','dsRedRecCas')){
            # Load scar probabilities. Because bulk scar sequencing is done differently than
            # single-cell scar sequencing, some single cell scars cannot be assigned a
            # probability because they cannot be observed in bulk sequencing. We filter 
            # these out.

            scar.probabilities <- scar_probabilities_in[scar_probabilities_in$Name == gene_it,]

            if(primer_set == "Nora"){      

                if(gene_it == "actb1"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,1,71)

                }else if(gene_it == "actb2"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,3,75)

                }else if(gene_it == "rpl39"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,1,62)

                }else if(gene_it == "rpl18a"){
                    scar.probabilities$Sequence <- substr(scar.probabilities$Sequence,1,75)

                }
            }


            duplicated_scars <- scar.probabilities[duplicated(scar.probabilities$Sequence) == T,]
            duplicated_scars <- scar.probabilities[scar.probabilities$Sequence %in% duplicated_scars$Sequence,]
            duplicated_scars <- duplicated_scars[order(duplicated_scars$Sequence),]
            length(unique(duplicated_scars$CIGAR))
            length(unique(duplicated_scars$Sequence))
            # Even before clipping the sequences in the scar-probabilities a little bit, there are duplicated sequences! This will lead to the creation of more than one row for a unique scar, when datasets are merged below! I'll therefore remove duplicated sequences, keeping the one that has the highest number of "Embryos".

            scar.probabilities <- scar.probabilities[order(scar.probabilities$Embryos, decreasing = T),]
            scar.probabilities <- scar.probabilities[order(scar.probabilities$Sequence, decreasing = T),]
            scar.probabilities <- scar.probabilities[duplicated(scar.probabilities$Sequence) == F,]

            # MERGE
            unique.scars_2 <- merge(unique.scars, scar.probabilities[, c("Sequence", "p", "Embryos")],
                                  by.x = "Sequence.short", by.y = "Sequence", all.x = T)
            unique.scars_2$Embryos[is.na(unique.scars_2$Embryos)] <- 0
            unique.scars_2$p[is.na(unique.scars_2$p)] <- min(unique.scars_2$p, na.rm = T)/100
            unique.scars_2$p[unique.scars_2$Sequence.short == wildtype.clip] <- max(unique.scars_2$p, na.rm = T) + 0.1
            unique.scars_2 <- unique.scars_2[order(-unique.scars_2$Presence, -unique.scars_2$p), ]
            unique.scars_2$Scar <- paste(0:(nrow(unique.scars_2)-1), unique.scars_2$CIGAR, sep = ":")
            
        }
        
        # Count total and unique scars per library ####
        unique.scars.nowt <- unique.scars_2[-1, ]
        scars.Z1 <- unique.scars.nowt$Scar[unique.scars.nowt$Freq.Z1 >= min.presence]

        sum(scars.Z1 %in% unique.scars.nowt$Scar[unique.scars.nowt$Presence == 1])

        # Write output ####
        Z1.scars.compared <- 
          merge(Z1.scars, unique.scars_2[, c("Sequence", "Presence", "p", "Embryos", "Scar")], by = "Sequence", all.x = T)
         
        write.csv(Z1.scars.compared, file = paste0("Z1_scars_compared_",dat_name,"_",gene_it,".csv"),
                  quote = F, row.names = F)
                
    }

}



