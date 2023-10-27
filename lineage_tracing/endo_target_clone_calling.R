# Functions for calling clones using sequencing data individual endogenous CRISPR-Cas9 lineage tracing targets:


# Written by Nora Fresmann 2022/2023

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))
require(plyr)
require(data.table)
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(pheatmap))


# 1. Call clones on a per-target basis using 

getclones <- function(dat_wd, dat_name, wt_seq_ids, seq_id_delim, scar_file_name, heatmap_after_initial_filt = F, heatmap_after_final_filt = F, cooc_fraction_cutoff, seur_dat, tums_list, genes_list){

# Loop through samples    
    for(z in 1:length(tums_list)){

        tum_it <- tums_list[z]

# Loop through genes
        for(y in 1:length(genes_list)){
    
            gene_it <- genes_list[y]        
            
# STEP 1: Data preparation           
            # Read in filtered scar data

            scar.input_full <- read.csv(paste0(dat_wd,"/",scar_file_name,"_",gene_it,".csv"), stringsAsFactors = F)

            scar.input_full <- merge(scar.input_full, seur_dat[,c('Barcode','sample_inf')], by = "Barcode")
            scar_freqs <- as.data.frame(table(scar.input_full$Sequence))
            scar.input_full$Presence <- 0
            scar.input_full$Presence[scar.input_full$Sequence %in% scar_freqs[scar_freqs$Freq > 1,]$Var1] <- 1
            
            print(paste0('entries in scar input object for gene ',gene_it,": ",nrow(scar.input_full)))
            
            scar.input_full$UMIs <- as.numeric(scar.input_full$UMIs)
            
            # Subset for cells from a single sample of interest, if per-sample analysis is desired
            if(tum_it == 'all'){
                scar.input <- scar.input_full
            }else{
                scar.input <- scar.input_full[scar.input_full$sample_inf == tum_it,]   
            }

            if(nrow(scar.input) == 0){
                next
            }



            # Get seq_ids, which are short stretches of sequence around the expected scar site. These will be used as identifiers for each sequence.
            # The length and positioning are defined by the input variable seq_id_delim and need to correspond to the supplied wildtype IDs in wt_seq_ids.
            scar.input$seq_id <- substr(scar.input$Sequence, seq_id_delim[[gene_it]][1],seq_id_delim[[gene_it]][2])
            
            # Get wildtype seq_id from wt_seq_ids.
            wt_seq_id <- wt_seq_ids[gene_it]                   

            # Include hard-clipped reads by giving them a different seq_id, unless they have only one UMI.
            # Then remove entries that still have an empty seq_id slot
            scar.input$seq_id[scar.input$seq_id == "" & scar.input$UMIs > 1] <- substr(scar.input$Sequence[scar.input$seq_id == "" & scar.input$UMIs > 1], 1,30)
            scar.input <- scar.input[scar.input$seq_id != "",]
            

            # Remove cell barcode seq_id duplicates
            scar.input <- scar.input[order(scar.input$UMIs, decreasing = T),]
            scar.input$BC_seq_id <- paste0(scar.input$Barcode, "_", scar.input$seq_id)
            scar.input <- scar.input[!duplicated(scar.input$BC_seq_id),]
            scar.input$BC_seq_id  <- NULL


# STEP 2: Selection of cells with two distinct alleles
            # First select cell barcodes that have been assigned two scars (wildtype sequences are allowed!)
            BCs_2scars <- scar.input$Barcode[duplicated(scar.input$Barcode) == T]
            scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]
            
            # Remove sequences that only appear in one cell. Since some cells that previously had two different alleles, might now only have one, remove cells that have been left with only a single scar.
            # Iterate through these two steps until there are no more scars that can only be found in a single cell. If this does not resolve within 20 rounds of iteration, send a notification and continue to next gene.

            Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)

            seqs_one_cell <- as.data.frame(t(as.matrix(colSums(Scar_clone_overview))))
            seqs_one_cell <- colnames(seqs_one_cell[,seqs_one_cell[1,] == 1])

            for(i in 1:20){
                scar.input.2scars <- scar.input.2scars[!scar.input.2scars$seq_id %in% seqs_one_cell,]

                BCs_2scars <- scar.input.2scars$Barcode[duplicated(scar.input.2scars$Barcode) == T]
                scar.input.2scars <- scar.input[scar.input$Barcode %in% BCs_2scars,]
                nrow(scar.input.2scars) == 2*length(BCs_2scars)

                Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$seq_id)
                seqs_one_cell <- colSums(Scar_clone_overview)
                seqs_one_cell <- names(which(seqs_one_cell == 1))

                if(identical(seqs_one_cell,character(0))  == T) {
                  break
                }
            }
            # MOVE TO NEXT gene, IF this converges to 20
            if(i == 20){
                print('Isolation of cells with two alleles that are also found in other cells did not converge. Skipping clone calling for gene ',gene_it,'.')
                next
            }

            # Hierarchically cluster cells and make a heatmap, if desired.
            if(heatmap_after_initial_filt == T){
                # Make a heatmap
                Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars$Barcode, scar.input.2scars$seq_id))
                scar_freqs <- sort(colSums(Scar_clone_overview), decreasing = T)

                ploto <-  pheatmap(Scar_clone_overview,
                  show_rownames = FALSE, show_colnames = FALSE,
                  breaks = seq(-1, +1, length = 101),
                  cluster_cols = F
                )

                ggsave(plot = ploto, filename = paste0(dat_wd,"/scarsintwocells_heatmap_",dat_name,"_tum_",tum_it,"_",gene_it,".png"), width = 7, height = 7, dpi = 300)
                rm(ploto)
            }


# STEP 3: Definition of scar pairs

            # Split the data into two dataframes. One contains the first allele of a cell and the second one contains the second allele of the cell.
            # We do this by sorting for seq_id first and then for the barcode to make sure that all cells with the same seq_id combination will have the two seq_ids appear in the same order in the dataframe.
            scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$seq_id),]
            scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$Barcode),]
            cell_scar_pairs <- scar.input.2scars[!duplicated(scar.input.2scars$Barcode),c("Barcode","seq_id","Sequence","CIGAR","UMIs","UMI_fraction")]
            colnames(cell_scar_pairs) <- c("Barcode","seq_id_1","Sequence_1","CIGAR_1","UMIs_1","UMI_fraction_1")

            cell_scar_pairs_2 <- scar.input.2scars[duplicated(scar.input.2scars$Barcode),c("Barcode","seq_id","Sequence","CIGAR","UMIs","UMI_fraction")]
            colnames(cell_scar_pairs_2) <- c("Barcode","seq_id_2","Sequence_2","CIGAR_2","UMIs_2","UMI_fraction_2")
            
            # Merge these two dataframes, so that each cell has a row with the two seq_ids in the columns seq_id_1 and seq_id_2
            cell_scar_pairs <- merge(cell_scar_pairs, cell_scar_pairs_2, by = "Barcode")
            cell_scar_pairs$UMI_sum <- cell_scar_pairs$UMIs_1 + cell_scar_pairs$UMIs_2

            # Remove combinations of seq_ids that only occur once or twice. If no seq_id combination is found more than once, send a note and move to next gene.
            cell_scar_pairs$comb_seq_id <- paste0(cell_scar_pairs$seq_id_1, cell_scar_pairs$seq_id_2) # Combine both seq_ids in a joint column
            comb_seq_ids <- as.data.frame(table(cell_scar_pairs$comb_seq_id))
            comb_seq_ids_over2 <- comb_seq_ids[comb_seq_ids$Freq > 2,]$Var1
            cell_scar_pairs_filt <- cell_scar_pairs[cell_scar_pairs$comb_seq_id %in% comb_seq_ids_over2,]

            if(nrow(cell_scar_pairs_filt) == 0){
                print(paste0("Sample ",tum_it,", gene ",gene_it,": There were no combinations that occur in more than 2 cells. Moving to the next gene"))
                next
            }
   
# STEP 4: Calculating scar co-occurrence for all scar combinations.

            # Go through all scar combinations that you see, using the cell_scar_pairs_filt dataframe as input.

            scar.input.2scars.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt$Barcode,]

            Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars.filt$Barcode, scar.input.2scars.filt$seq_id))

            for(i in 1:length(unique(cell_scar_pairs_filt$comb_seq_id))){

              comb_seq_id <- unique(cell_scar_pairs_filt$comb_seq_id)[i]
              scar_1 <- unique(cell_scar_pairs_filt$seq_id_1[cell_scar_pairs_filt$comb_seq_id == comb_seq_id])
              scar_2 <- unique(cell_scar_pairs_filt$seq_id_2[cell_scar_pairs_filt$comb_seq_id == comb_seq_id])

              subs <- Scar_clone_overview[, c(scar_1, scar_2)]
              subs <- subs[rowSums(subs) != 0,]

              scar_freq_1 <- length(subs[,1][subs[,1] == 1])
              scar_freq_2 <- length(subs[,2][subs[,2] == 1])

              subs$common <- rowSums(subs)
              shared_freq <- nrow(subs[subs$common == 2,])

              # calculate fraction of combinations of the two scars of the total combinations that the two scars contribute to.
              frac_scar_1 <- shared_freq/scar_freq_1
              frac_scar_2 <- shared_freq/scar_freq_2


              if(i == 1){
              cell_scar_pairs_filt$combo_fraction_1 <- 0
              cell_scar_pairs_filt$combo_fraction_2 <- 0
              }
              cell_scar_pairs_filt$combo_fraction_1[cell_scar_pairs_filt$comb_seq_id == comb_seq_id] <- frac_scar_1
              cell_scar_pairs_filt$combo_fraction_2[cell_scar_pairs_filt$comb_seq_id == comb_seq_id] <- frac_scar_2
            }

# STEP 5: Removing scar-pairs with scars that are ambiguous due to co-occurrence with multiple other scars. Cut-off defined by cooc_fraction_cutoff.

            # Consider all combinations, which make up over X % (defined by cooc_fraction_cutoff) of observations of one of the respective scars
            # If no scar pairs pass the filter, send a note and move to next gene.
            cell_scar_pairs_filt_filt <- cell_scar_pairs_filt[cell_scar_pairs_filt$combo_fraction_1 > cooc_fraction_cutoff | cell_scar_pairs_filt$combo_fraction_2 > cooc_fraction_cutoff ,]

            if(nrow(cell_scar_pairs_filt_filt) == 0){
                print(paste0("Sample ",tum_it,", gene ",gene_it,": There were no scar combinations that can confidently be called a true pair. Moving to the next gene"))
                next
            }        
            
            # Subset the filtered original scar data for the cells that passed the last filter
            scar.input.2scars.filt.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt_filt$Barcode,]


            if(heatmap_after_final_filt == TRUE){
                # Plot these scars in a heatmap
                Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars.filt.filt$Barcode, scar.input.2scars.filt.filt$seq_id))
                scar_freqs <- sort(colSums(Scar_clone_overview), decreasing = T)

                ploto <- pheatmap(Scar_clone_overview,
                                  show_rownames = FALSE, show_colnames = FALSE,
                                  breaks = seq(-1, +1, length = 101)
                                 )
                ggsave(plot = ploto, filename = paste0(dat_wd,"/scarsintwocells_seqcombs_over_",(cooc_fraction_cutoff*10),"_",dat_name,"_tum_",tum_it,"_",gene_it,".png"), width = 6, height = 6, dpi = 300)
                rm(ploto)
            }

# STEP 6: Classify each scar in a pair based on whether it was likely generated first ('parent scar') or last ('child scar')

            # Define the "child" scar that actually defines the clone -  the other might be a "parent" scar that also occurs in combinations with other scars.
            # For these, we'll define the clone-defining scar as the one that has a higher fraction of all its combinations with the paired scar.
            # For some scar combinations, this doesn't make much sense, as they are more or less exclusive to each other.
            # These get the label 'both' for the clone-defining scar.
            # Wildtype seq_ids are not allowed as child-sequences!

            # Define the "child", i.e. clone-defining scar
            cell_scar_pairs_filt_filt$clone_def_scar <- cell_scar_pairs_filt_filt$seq_id_1
            cell_scar_pairs_filt_filt$clone_def_scar[cell_scar_pairs_filt_filt$combo_fraction_2 > cell_scar_pairs_filt_filt$combo_fraction_1] <- cell_scar_pairs_filt_filt$seq_id_2[cell_scar_pairs_filt_filt$combo_fraction_2 > cell_scar_pairs_filt_filt$combo_fraction_1] 
            cell_scar_pairs_filt_filt$clone_def_scar[abs(cell_scar_pairs_filt_filt$combo_fraction_1 - cell_scar_pairs_filt_filt$combo_fraction_2) < 0.05] <- "both"

            cell_scar_pairs_filt_filt <- cell_scar_pairs_filt_filt[cell_scar_pairs_filt_filt$clone_def_scar != wt_seq_id,]

            if(nrow(cell_scar_pairs_filt_filt) == 0){
                print(paste0("Sample ",tum_it,", gene ",gene_it,": No scars retained after classifying child scars in scar combinations. Moving to next gene."))
                next
            }        
                            
            # Define the "parent" scar
            cell_scar_pairs_filt_filt$parent_scar <- "both"
            cell_scar_pairs_filt_filt$parent_scar[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_1] <- cell_scar_pairs_filt_filt$seq_id_2[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_1]
            cell_scar_pairs_filt_filt$parent_scar[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_2] <- cell_scar_pairs_filt_filt$seq_id_1[cell_scar_pairs_filt_filt$clone_def_scar == cell_scar_pairs_filt_filt$seq_id_2]

# STEP 7: Merge relevant data and write to file.
            
            # Add information on the comined seq_ids defining clones to the scar data
            clone_seqs <- as.data.frame(unique(cell_scar_pairs_filt_filt$comb_seq_id))
            colnames(clone_seqs)[1] <- "clone_seq"

            clone_seqs$cloneID <- c(1:length(unique(clone_seqs$clone_seq)))
            clone_seqs$cloneID <- paste0(gene_it,"_",clone_seqs$cloneID)

            cell_scar_pairs_filt_filt <- merge(cell_scar_pairs_filt_filt, clone_seqs, by.x = "comb_seq_id", by.y = "clone_seq", all.x = T)
            
            # Write extensive output to file
            write.table(cell_scar_pairs_filt_filt, paste0(dat_wd,"/final_clones_allinfo_",dat_name,"_tum_",tum_it,"_",gene_it,"_cutoff_",cooc_fraction_cutoff,".txt"), quote = F)

            # Add metadata information and reduce individual scar-cell information to most relevant.
            cell_scar_pairs_filt_filt_sub <- cell_scar_pairs_filt_filt[c("comb_seq_id","Barcode","seq_id_1","seq_id_2", "cloneID","clone_def_scar","parent_scar")]

            all_info <- merge(scar.input.2scars.filt.filt, cell_scar_pairs_filt_filt_sub, by = "Barcode", all.x = T)
            
            # Write analysis output to file
            write.table(all_info, paste0(dat_wd,"/final_clones_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",gene_it,"_cutoff_",cooc_fraction_cutoff,".txt"), quote = F, sep = ",")

            rm(all_info)
            
# STEP 8: Add cells that only have one scar in the data, when this scar was defined as a clone-defining scar in the previous analysis steps.

            # Load data back in (into variable name scar_dat)
            scar_dat <- read.delim(paste0(dat_wd,"/final_clones_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",gene_it,"_cutoff_",cooc_fraction_cutoff,".txt"),stringsAsFactors = F, sep = ',', row.names = 1, header = T)

            # Select cell barcodes that only have one scar
            BCs_1scar <- as.data.frame(table(scar.input$Barcode))
            BCs_1scar <- BCs_1scar$Var1[BCs_1scar$Freq == 1]
            scar.input.1scar <- scar.input[scar.input$Barcode %in% BCs_1scar,]

            # Select cell barcodes that have a clone-defining scar
            scar.input.1scar <- scar.input.1scar[scar.input.1scar$seq_id %in% scar_dat$clone_def_scar,]

            # Only keep cell barcodes that have more than one UMI for that scar
            scar.input.1scar <- scar.input.1scar[scar.input.1scar$UMIs > 1,]

            
            # If no cells passed the filters above, move to next gene. Otherwise, fill in information for cells to merge with the previous analysis output.
            if(nrow(scar.input.1scar) != 0){

                # Infer or fill in additional columns that are needed to merge this data with the scar-data for cells with two scars.
                # Add cell barcode with corresponding parent scar!

                scar.input.1scar$comb_seq_id <- NA
                scar.input.1scar$seq_id_1 <- scar.input.1scar$seq_id
                scar.input.1scar$seq_id_2 <- NA

                # Get a table of clone names and corresponding defining scar

                clone_def <- scar_dat[duplicated(scar_dat$cloneID) == F,]
                clone_def$seq_id <- clone_def$clone_def_scar
                clone_def <- clone_def[,c("seq_id","cloneID", "clone_def_scar", "parent_scar")]
                
                # Merge the two
                scar.input.1scar <- left_join(scar.input.1scar, clone_def, by = "seq_id")
                    dim(scar.input.1scar)

                # PRINT OUT WARNINGS, IF SOMETHING IS WRONG HERE
                if(TRUE %in% names(table(scar.input.1scar$Barcode %in% scar_dat$Barcode))){
                    print("Warning: WRONG table merger 1. Please check your input and the output from the two-scar-clone analysis!")  
                }else if(FALSE %in% table(colnames(scar.input.1scar) == colnames(scar_dat))){
                    print("Warning: WRONG table merger 2. Please check your input and the output from the two-scar-clone analysis!")  
                }
                
                # Merge data from previous step with the additional cells recovered here.
                scar_dat_ext <- rbind(scar_dat, scar.input.1scar)
                print(paste0("Sample ",tum_it,", gene ",gene_it,": ","adding ",nrow(scar.input.1scar)," cells with a single child scar"))
                
                # Write data to file
                write.table(scar_dat_ext, paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",gene_it,"_cutoff_",cooc_fraction_cutoff,".txt"), quote = F, sep = ",")
                rm(scar_dat_ext)

            }else{
                
                print(paste0("Sample ",tum_it,", gene ",gene_it,": No additional cells with only one scar could be classified! Moving to next gene."))
                
            }
            
            
            
        }

    }
}

                



## Merge on target pairs
merge_target_pairs <- function(dat_wd, dat_name, cooc_fraction_cutoff, tums_list, genes_list){

    gene_combs <- as.data.frame(t(combn(genes_list, 2)))

    for(z in 1:length(tums_list)){

        tum_it <- tums_list[z]
        
        for(y in 1:nrow(gene_combs)){
            
            targ_1 <- gene_combs[y,1]
            targ_2 <- gene_combs[y,2]

# STEP 1: Load in data for two target genes.
        
            ## Load the first two scar datasets

            if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_withUncert.txt"))){
                scars_1 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_withUncert.txt"), stringsAsFactors = F, sep = ',')
                
            }else if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,"_SingleScarExt.txt"))){
                scars_1 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,"_SingleScarExt.txt"), stringsAsFactors = F, sep = ',')

            }else if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,"_with_wt.txt"))){
                scars_1 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,"_with_wt.txt"), stringsAsFactors = F, sep = ',')

            }else if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,".txt"))){
                scars_1 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,".txt"), stringsAsFactors = F, sep = ',')

            }else if(file.exists(paste0(dat_wd,"/final_clones_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,".txt"))){
                scars_1 <- read.delim(paste0(dat_wd,"/final_clones_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_1,"_cutoff_",cooc_fraction_cutoff,".txt"), stringsAsFactors = F, sep = ',')

            }else{
                print(paste0(tum_it," ,target ",targ_1," not found"))
                next
            }


            if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_withUncert.txt"))){
                scars_2 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_withUncert.txt"), stringsAsFactors = F, sep = ',')
                
            }else if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,"_SingleScarExt.txt"))){
                scars_2 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,"_SingleScarExt.txt"), stringsAsFactors = F, sep = ',')

            }else if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,"_with_wt.txt"))){
                scars_2 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,"_with_wt.txt"), stringsAsFactors = F, sep = ',')

            }else if(file.exists(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,".txt"))){
                scars_2 <- read.delim(paste0(dat_wd,"/final_clones_with_singlescarcells_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,".txt"), stringsAsFactors = F, sep = ',')

            }else if(file.exists(paste0(dat_wd,"/final_clones_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,".txt"))){
                scars_2 <- read.delim(paste0(dat_wd,"/final_clones_allinfo_per_cell_scar_",dat_name,"_tum_",tum_it,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,".txt"), stringsAsFactors = F, sep = ',')

            }else{
                print(paste0(tum_it," ,target ",targ_2," not found"))
                next
            }
                
            
            # Keep one row per barcode
            scars_1 <- scars_1[duplicated(scars_1$Barcode) == F,]
            scars_2 <- scars_2[duplicated(scars_2$Barcode) == F,]

            # Merge the two datasets
            all_scars <- rbind(scars_1,scars_2)

# STEP 2: Selection of cells that have clone information for both targets

            BCs_2scars <- all_scars$Barcode[duplicated(all_scars$Barcode) == T]
            scar.input.2scars <- all_scars[all_scars$Barcode %in% BCs_2scars,]

            # Remove sequences that only appear in one cell. Since some cells that previously had two different alleles, might now only have one, remove cells that have been left with only a single scar.
            # Iterate through these two steps until there are no more scars that can only be found in a single cell. If this does not resolve within 20 rounds of iteration, send a notification and continue to next gene.
            Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$cloneID)

            seqs_one_cell <- colSums(Scar_clone_overview)
            seqs_one_cell <- names(which(seqs_one_cell == 1))

            for(i in 1:20){
                scar.input.2scars <- scar.input.2scars[!scar.input.2scars$cloneID %in% seqs_one_cell,]

                BCs_2scars <- scar.input.2scars$Barcode[duplicated(scar.input.2scars$Barcode) == T]
                scar.input.2scars <- all_scars[all_scars$Barcode %in% BCs_2scars,]
                
                Scar_clone_overview <- table(scar.input.2scars$Barcode, scar.input.2scars$cloneID)
                seqs_one_cell <- colSums(Scar_clone_overview)
                seqs_one_cell <- names(which(seqs_one_cell == 1))

                if(identical(seqs_one_cell,character(0))  == T) {
                  break
                }
            }

            # MOVE TO NEXT tumour, IF this converges to 20
            if(i == 20){
                print(paste0('Isolation of cells with clone information on both targets did not converge. Skipping target-pair-based clone calling for genes ',targ_1, '_', targ_2,'.'))
                next
            }

            Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars$Barcode, scar.input.2scars$cloneID))
            scar_freqs <- sort(colSums(Scar_clone_overview), decreasing = T)


# STEP 3: Definition of target-clone pairs

            # Split the data into two dataframes. One contains the first allele of a cell and the second one contains the second allele of the cell.
            # We do this by sorting for cloneID first and then for the barcode to make sure that all cells with the same cloneID combination will have the two cloneIDs appear in the same order in the dataframe.
            scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$cloneID),]
            scar.input.2scars <- scar.input.2scars[order(scar.input.2scars$Barcode),]
            cell_scar_pairs <- scar.input.2scars[!duplicated(scar.input.2scars$Barcode),c("Barcode","comb_seq_id","seq_id_1","seq_id_2","cloneID","clone_def_scar","CIGAR")]
            colnames(cell_scar_pairs) <- c("Barcode", paste0("targ_1_",c("comb_seq_id","seq_id_1","seq_id_2","cloneID","clone_def_scar","CIGAR")))

            cell_scar_pairs_2 <- scar.input.2scars[duplicated(scar.input.2scars$Barcode),c("Barcode","comb_seq_id","seq_id_1","seq_id_2","cloneID","clone_def_scar","CIGAR")]
            colnames(cell_scar_pairs_2) <- c("Barcode", paste0("targ_2_",c("comb_seq_id","seq_id_1","seq_id_2","cloneID","clone_def_scar","CIGAR")))

            cell_scar_pairs <- merge(cell_scar_pairs, cell_scar_pairs_2, by = "Barcode")

            # Remove cloneID-combinations that only occur once:
            cell_scar_pairs$targ_1_targ_2_clone_id <- paste0(cell_scar_pairs$targ_1_cloneID,"_", cell_scar_pairs$targ_2_cloneID)
            targ_1_targ_2_clone_ids <- as.data.frame(table(cell_scar_pairs$targ_1_targ_2_clone_id))

            targ_1_targ_2_clone_ids_over2 <- targ_1_targ_2_clone_ids[targ_1_targ_2_clone_ids$Freq > 1,]$Var1
            cell_scar_pairs_filt <- cell_scar_pairs[cell_scar_pairs$targ_1_targ_2_clone_id %in% targ_1_targ_2_clone_ids_over2,]
            dim(cell_scar_pairs_filt)


            scar.input.2scars.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt$Barcode,]
            
            if(nrow(scar.input.2scars.filt) == 0){
                next
            }
            
# STEP 4: Calculating cloneID co-occurrence for all cloneID combinations.
            
            Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars.filt$Barcode, scar.input.2scars.filt$cloneID))

            for(i in 1:length(unique(cell_scar_pairs_filt$targ_1_targ_2_clone_id))){

              targ_1_targ_2_clone_id <- unique(cell_scar_pairs_filt$targ_1_targ_2_clone_id)[i]
              scar_1 <- unique(cell_scar_pairs_filt$targ_1_cloneID[cell_scar_pairs_filt$targ_1_targ_2_clone_id == targ_1_targ_2_clone_id])
              scar_2 <- unique(cell_scar_pairs_filt$targ_2_cloneID[cell_scar_pairs_filt$targ_1_targ_2_clone_id == targ_1_targ_2_clone_id])

              subs <- Scar_clone_overview[, c(scar_1, scar_2)] # error!
              subs <- subs[rowSums(subs) != 0,]

              scar_freq_1 <- length(subs[,1][subs[,1] == 1])
              scar_freq_2 <- length(subs[,2][subs[,2] == 1])

              subs$common <- rowSums(subs)
              shared_freq <- nrow(subs[subs$common == 2,])

              # calculate fraction of combinations of the two scars of the total combinations that the two cloneIDs contribute to.
              frac_scar_1 <- shared_freq/scar_freq_1
              frac_scar_2 <- shared_freq/scar_freq_2


              if(i == 1){
              cell_scar_pairs_filt$combo_fraction_targ_1 <- 0
              cell_scar_pairs_filt$combo_fraction_targ_2 <- 0
              }

              cell_scar_pairs_filt$combo_fraction_targ_1[cell_scar_pairs_filt$targ_1_targ_2_clone_id == targ_1_targ_2_clone_id] <- frac_scar_1
              cell_scar_pairs_filt$combo_fraction_targ_2[cell_scar_pairs_filt$targ_1_targ_2_clone_id == targ_1_targ_2_clone_id] <- frac_scar_2
            }


# STEP 5: Removing scar-pairs with scars that are ambiguous due to co-occurrence with multiple other scars. Cut-off defined by cooc_fraction_cutoff.

            # Consider all combinations, which make up over X % of observations of one of the respective cloneIDs
            cell_scar_pairs_filt_filt <- cell_scar_pairs_filt[cell_scar_pairs_filt$combo_fraction_targ_1 > 0.8 | cell_scar_pairs_filt$combo_fraction_targ_2 > 0.8 ,]

            if(nrow(cell_scar_pairs_filt_filt) == 0){
                print(paste0("Sample ",tum_it,", genes ",targ_1,"_",targ_2,": There were no scar combinations that can confidently be called a true pair. Moving to the next gene"))
                next
            }

            scar.input.2scars.filt.filt <- scar.input.2scars[scar.input.2scars$Barcode %in% cell_scar_pairs_filt_filt$Barcode,]

            Scar_clone_overview <- as.data.frame.matrix(table(scar.input.2scars.filt.filt$Barcode, scar.input.2scars.filt.filt$cloneID))
            scar_freqs <- sort(colSums(Scar_clone_overview), decreasing = T)


# STEP 6: Classify each cloneID in a pair based on whether it was likely generated first ('parent cloneID') or last ('child cloneID')

            # Define the "child" cloneID that actually defines the clone -  the other might be a "parent" cloneID that also occurs in combinations with other cloneIDs
            # For these, we'll define the clone-defining cloneID as the one that has a higher fraction of all its combinations with the paired cloneID
            # For some cloneID combinations, this doesn't make much sense, as they are more or less exclusive to each other.
            # These get the label 'both' for the clone-defining cloneID
            # Wildtype seq_ids are not allowed as child-sequences!


            # 1. Define the "child", i.e. clone-defining scar

            cell_scar_pairs_filt_filt$clone_def_cloneID <- cell_scar_pairs_filt_filt$targ_1_cloneID
            cell_scar_pairs_filt_filt$clone_def_cloneID[cell_scar_pairs_filt_filt$combo_fraction_targ_2 > cell_scar_pairs_filt_filt$combo_fraction_targ_1] <- cell_scar_pairs_filt_filt$targ_2_cloneID[cell_scar_pairs_filt_filt$combo_fraction_targ_2 > cell_scar_pairs_filt_filt$combo_fraction_targ_1] 
            cell_scar_pairs_filt_filt$clone_def_cloneID[abs(cell_scar_pairs_filt_filt$combo_fraction_targ_1 - cell_scar_pairs_filt_filt$combo_fraction_targ_2) < 0.05] <- "both"

            # Define the "parent" clone

            cell_scar_pairs_filt_filt$parent_cloneID <- "both"
            cell_scar_pairs_filt_filt$parent_cloneID[cell_scar_pairs_filt_filt$clone_def_cloneID == cell_scar_pairs_filt_filt$targ_1_cloneID] <- cell_scar_pairs_filt_filt$targ_2_cloneID[cell_scar_pairs_filt_filt$clone_def_cloneID == cell_scar_pairs_filt_filt$targ_1_cloneID]
            cell_scar_pairs_filt_filt$parent_cloneID[cell_scar_pairs_filt_filt$clone_def_cloneID == cell_scar_pairs_filt_filt$targ_2_cloneID] <- cell_scar_pairs_filt_filt$targ_1_cloneID[cell_scar_pairs_filt_filt$clone_def_cloneID == cell_scar_pairs_filt_filt$targ_2_cloneID]



            cell_scar_pairs_filt_filt_sub <- cell_scar_pairs_filt_filt[c("targ_1_targ_2_clone_id","Barcode","targ_1_cloneID","targ_2_cloneID","clone_def_cloneID", "parent_cloneID")]

            colnames(cell_scar_pairs_filt_filt_sub) <- gsub("targ_1",paste0(targ_1),colnames(cell_scar_pairs_filt_filt_sub))
            colnames(cell_scar_pairs_filt_filt_sub) <- gsub("targ_2",paste0(targ_2),colnames(cell_scar_pairs_filt_filt_sub))


# STEP 7: Merge relevant data and write to file.

            all_info <- merge(scar.input.2scars.filt.filt, cell_scar_pairs_filt_filt_sub, by = "Barcode", all.x = T)

            write.table(all_info, paste0(dat_wd,"/final_combclones_",dat_name,"_tum_",tum_it,"_",targ_1,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,".txt"), quote = F, sep = ",")

            
# STEP 8: Add cells that only have one cloneID (within the two genes assessed here) in the data, when this scar was defined as a clone-defining cloneID in the previous analysis steps.

            # Write the output from the previous steps into the variable scar_dat
            scar_dat <- all_info
            rm(all_info)

            # Define cloneID pairs in order to be able to impute cloneIDs later on
            scarpairs <- as.data.frame(table(scar_dat[,paste0(targ_1,"_cloneID")], scar_dat[,paste0(targ_2,"_cloneID")]))
            scarpairs_1 <- scarpairs[scarpairs$Freq > 0,]
            scarpairs_2 <- scarpairs_1

            rm(scarpairs)
            scarpairs_2$Var1 <- scarpairs_1$Var2
            scarpairs_2$Var2 <- scarpairs_1$Var1

            scarpairs <- rbind(scarpairs_1,scarpairs_2)
            scarpairs[,3] <- NULL

            colnames(scarpairs) <- c("cloneID_1","cloneID")

            scar_dat_min <- scar_dat
            scar_dat_min$Barcode <- NULL
            all_info_pairs <- left_join(scarpairs, scar_dat_min, by = "cloneID")
            all_info_pairs <- all_info_pairs[!duplicated(all_info_pairs$cloneID_1),]


            # Select cell barcodes that only have one cloneID
            BCs_1scar <- as.data.frame(table(all_scars$Barcode))
            BCs_1scar <- BCs_1scar$Var1[BCs_1scar$Freq == 1]
            scar.input.1scar <- all_scars[all_scars$Barcode %in% BCs_1scar,]


            defclones <- c(scar_dat$clone_def_cloneID, scar_dat$targ_1_cloneID[scar_dat$clone_def_cloneID == "both"], scar_dat$targ_2_cloneID[scar_dat$clone_def_cloneID == "both"])
            unique(defclones)

            # Select cell barcodes that have a clone-defining cloneID
            scar.input.1scar <- scar.input.1scar[scar.input.1scar$cloneID %in% defclones,]

            # IF THE DATAFRAME IS EMPTY, MOVE ON TO NEXT GENE
            if(nrow(scar.input.1scar) == 0){
                print(paste0("Sample ",tum_it,", gene combination ",targ_1, ", ",targ_2,": No additional cells with only one scar could be classified! Moving to next gene."))
                next
            }

            # Infer or fill in additional columns that are needed to merge this data with the scar-data for cells with two scars

            # Get a table of clone names and corresponding defining scar
            clone_def <- scar_dat[duplicated(scar_dat$cloneID) == F,]

            #clone_def$cloneID <- clone_def$clone_def_cloneID
            clone_def$targ_1_cloneID <- NA
            clone_def$targ_2_cloneID <- NA
            clone_def$targ_1_cloneID[clone_def$cloneID %like any% "targ_1%"] <- clone_def$cloneID[clone_def$cloneID %like any% "targ_1%"]
            clone_def$targ_2_cloneID[clone_def$cloneID %like any% "targ_2%"] <- clone_def$cloneID[clone_def$cloneID %like any% "targ_2%"]
            clone_def <- clone_def[,c("cloneID",paste0(targ_1,"_",targ_2,"_clone_id"),paste0(targ_1, "_cloneID"),paste0(targ_2, "_cloneID"),"clone_def_cloneID", "parent_cloneID")]

            # Merge the two
            scar.input.1scar <- left_join(scar.input.1scar, clone_def, by = "cloneID")

            colnames(scar.input.1scar) <- gsub("targ_1",paste0(targ_1),colnames(scar.input.1scar))
            colnames(scar.input.1scar) <- gsub("targ_2",paste0(targ_2),colnames(scar.input.1scar))


            # Impute the second clone row needed for the combined clone
            scar.input.1scar_2 <- scar.input.1scar
            cols_pre <- colnames(scar.input.1scar_2)
            scar.input.1scar_2 <- scar.input.1scar_2[,c("Barcode","cloneID")]
            colnames(scar.input.1scar_2) <- c("Barcode","cloneID_1")



            scar.input.1scar_2 <- left_join(scar.input.1scar_2, all_info_pairs, by = "cloneID_1")
            scar.input.1scar_2$cloneID_1 <- NULL
            scar.input.1scar_2 <- scar.input.1scar_2[,cols_pre]

            scar.input.1scar_2$UMIs <- "imputed"
            scar.input.1scar_2$UMI_fraction <- "imputed"
            scar.input.1scar_2$UMI_cfraction <- "imputed"
            scar.input.1scar_2$Fraction <- "imputed"
            scar.input.1scar_2$UMI_cf_min <- "imputed"
            scar.input.1scar_2$Presence <- "imputed"
            scar.input.1scar_2$p <- "imputed"
            scar.input.1scar_2$Embryos <- "imputed"

            table(colnames(scar.input.1scar_2) == cols_pre)
            table(duplicated(scar.input.1scar_2$Barcode))


            scar.input.1scar <- rbind(scar.input.1scar,scar.input.1scar_2)


            # PRINT OUT WARNINGS, IF SOMETHING IS WRONG HERE
            if(TRUE %in% table(scar.input.1scar$Barcode %in% scar_dat$Barcode)){
                print("WRONG table merger. Check input, output and code!")  
            }else if(FALSE %in% table(colnames(scar.input.1scar) == colnames(scar_dat))){
                print("WRONG table merger. Check input, output and code!")  
            }

            scar_dat_ext <- rbind(scar_dat, scar.input.1scar)

            if(FALSE %in% (nrow(scar_dat_ext) == (nrow(scar_dat) + nrow(scar.input.1scar)))){
                print("WRONG table merger.  Check input, output and code!")  
            }

            write.table(scar_dat_ext, paste0(dat_wd,"/final_combclones_",dat_name,"_tum_",tum_it,"_",targ_1,"_",targ_2,"_cutoff_",cooc_fraction_cutoff,"_with_singletargcells.txt"), quote = F, sep = ",")


            rm(all_scars)
            rm(scars_1)
            rm(scars_2)
            

        }

    }

}



            

            

            

    