library("dplyr")
library("ggplot2")
library("tidyr")
library("betareg")
library("effectsize")
library("ggrepel")
library("scales")
library("rstatix")
library("metap")

setwd("/Users/shannon/uni_onedrive/Honours/1KG/analysis/")

x_chromosome_loci <- c("BCLAF3", "ARX_2", "ARX_1", "DMD", "AR", "ZIC3", "SOX3", "FMR1", "AFF2", "TMEM185A")
read_len <- 150
fig_dir <- "../figures"

get_genotype <- function(gt) {
  if (is.na(gt)) {
    return(NA)
  }
  if (as.character(gt) == "./." | as.character(gt) == ".") {
    return(NA)
  }
  return(as.numeric(unlist(strsplit(as.character(gt), "/"))))
}

get_all_genotypes <- function(gts) {
  lapply(gts, get_genotype)
}

get_other_parent <- function(df, parent, child) {
  return(df[df$ParentID != parent & df$ChildID == child,])
}

create_pileup_table <- function(df, filename) {
  df %>%
    mutate(Gene = replace(Gene, Gene == "NUTM2B.AS1", "NUTM2B-AS1")) %>%
    mutate(ParentPath = paste(paste(ParentID, ParentID, sep="/"), Gene, "svg", sep = ".")) %>%
    mutate(OtherParentPath = paste(paste(OtherParentID, OtherParentID, sep="/"), Gene, "svg", sep = ".")) %>%
    mutate(ChildPath = paste(paste(ChildID, ChildID, sep="/"), Gene, "svg", sep = ".")) %>%
    mutate(Directory = paste(ParentID, ParentAllele, ".", ChildID, ChildAllele, Gene, sep="_")) %>%
    dplyr::select(c(ParentPath, OtherParentPath, ChildPath, Directory)) %>%
    write.table(.,file = paste("../pileups/", filename, sep=""), sep=",", row.names = F, col.names = F, quote = F)
}

ggsave_wrap <- function(filename, width, height) {
  ggsave(paste(fig_dir, filename, sep="/"), width = width, height = height)
}

#--------------------------------- PREP ------------------------------------------------------

# Load data
df_og <- read.csv("str_disease_trios.csv")
info_df <- read.csv("../../catalog/disease_loci_info_df.csv")
ped <- read.csv("scripts/1kg_v3.ped", sep = "\t")
trb_df <- read.csv("1KGTRB.csv")

# Correct POLG variant names 
names(df_og)[names(df_og) == 'POLG_chr15.89333589.89333595'] <- 'POLG'
# names(df_og)[names(df_og) == 'POLG_chr15.89333589.89333595_CI'] <- 'POLG_CI'

# Correct NUTM2B-AS1
names(df_og)[names(df_og) == 'NUTM2B.AS1'] <- 'NUTM2B-AS1'
info_df[info_df$gene == "NUTM2B.AS1",]$gene <- "NUTM2B-AS1"

# Rename CIs
# names(df_og) <- lapply(names(df_og), gsub, pattern="_CI", replacement="-CI")

# Change HG02300 from female to male (incorrect information in pedigree)
df_og[df_og$SampleID == "HG02300",]$Sex <- 1

# Add family IDs and population
names(ped)[names(ped) == 'Individual.ID'] <- "SampleID"
names(ped)[names(ped) == 'Family.ID'] <- "FamilyID"
df_og <- merge(df_og, ped[, c("SampleID", "FamilyID", "Population")], by="SampleID")
df_og <- df_og %>%
  dplyr::select(FamilyID, Population, everything())

# Remove confidence intervals
# df <- select(df_og, !ends_with("-CI"))
df <- df_og

# Change cell values to [3, 3] instead of '3/3' etc
df <- df %>%
  dplyr::mutate( across(-c(FamilyID, Population, SampleID, CpgID, FatherID, MotherID, Sex), get_all_genotypes)
  )

# ------------------------------ Find missing values ---------------------------

na_counts <- colSums(is.na(df))
na_columns <- na_counts[na_counts > 0]
na_columns

only_genotypes_df <- df %>%
  dplyr::select(-c(FamilyID, Population, SampleID, CpgID, FatherID, MotherID, Sex))
sum(is.na(only_genotypes_df))

#--------------------------------- Summary Stats  ------------------------------

summary_stats <- function(locus, binwidth=1, limit = max(unlist(df[[locus]])), sampleid_subset) {
  values <- df %>% filter(SampleID %in% complete) %>% pull(locus) %>% unlist()
  print(length(values))
  plot <- ggplot(data.frame(values), aes(x=values, y=after_stat(count / sum(count)))) +
    geom_histogram(binwidth=binwidth, colour="black", fill="grey") +
    xlab("Repeat Length (RU)") +
    ylab("Frequency (%)") +
    theme_classic() +
    theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
    scale_y_continuous(labels = percent_format()) +
    xlim(c(0,limit))
  print(paste(locus, " summary stats", sep=""))
  print(summary(unlist(df[[locus]])))
  print(paste("SD: ", sd(unlist(df[[locus]])), sep=""))
  plot
}

# Distribution of XYLT1
summary_stats("XYLT1", limit=100)
summary_stats("RFC1", limit=700)
max(unlist(df[["RFC1"]]), na.rm = T)
length(unlist(df[["RFC1"]]))
summary_stats("FMR1")

#-------------------------------- Find quads -----------------------------------
# quads <- df %>%
#   filter(MotherID != 0 & FatherID != 0) %>%
#   group_by(MotherID, FatherID) %>%
#   filter(n() == 2) 
# nrow(quads) / 2
# quints <- df %>%
#   filter(MotherID != 0 & FatherID != 0) %>%
#   group_by(MotherID, FatherID) %>%
#   filter(n() == 3) 
# nrow(quints) / 3

#------------------ Find multi-generational transmissions ----------------------
df %>%
  filter(SampleID %in% MotherID | SampleID %in% FatherID) %>%
  filter(!(FatherID == 0 & MotherID == 0))

#------------------------ Find germline mutations ------------------------------

disease_loci <- colnames(df)[8:length(colnames(df))]
related <- df[df$SampleID %in% c(df$FatherID, df$MotherID) | 
                     (df$FatherID != 0 & df$MotherID != 0 & 
                        df$FatherID %in% df$SampleID & df$MotherID %in% df$SampleID),]
children <- related[related$FatherID != 0 & related$MotherID != 0,]
complete <- union(children$SampleID, union(children$FatherID, children$MotherID))

# Return list(c(allele from p1, corresponding allele in c), c(allele from p2, corresponding allele in c))
# Assume input has no NAs or 0s
# Mother is p1 and father is p2
# NA means edge case and cannot tell
# One return result means X chromosome gene and transmission from mother -> son
return_transmissions <- function(p1, p2, ch) {
  order <- c(NA, NA) # order[1]: parent who contributed c[1] and order[2]: parent who contributed c[2]
  
  # Returns the optimal assignment of child allele to parent while prioritising unmutated alleles
  get_optimal_order <- function(p1_p2_mins, p2_p1_mins, child_alleles) {
    if (child_alleles[1] == child_alleles[2]) {
      return(list(p1, p2)) # Order doesn't matter if child is homozygous
    }
    if (0 %in% p1_p2_mins && 0 %in% p2_p1_mins) { # Could have passed down an unchanged allele with either order
      if (max(p1_p2_mins) == max(p2_p1_mins)) { # Same length distances either way! Cannot tell -> edge case
        return(list(NA, NA))
      }
      if (max(p1_p2_mins) < max(p2_p1_mins)) {
        return(list(p1, p2))
      } else {
        return(list(p2, p1))
      }
    } else if (0 %in% p1_p2_mins) { # c[1] from mother and c[2] from father is optimal
      return(list(p1, p2))
    } else if (0 %in% p2_p1_mins) { # c[1] from father and c[2] from mother is optimal
      return(list(p2, p1))
    } else { # Double mutation! Take averages of minimum distances to decide order
        if (mean(p1_p2_mins) == mean(p2_p1_mins) ) { # Cannot decide order if both arrangements are optimal, edge case
        return(list(NA, NA)) 
      } else if (mean(p1_p2_mins) < mean(p2_p1_mins)) {
        return(list(p1, p2))
      } else {
        return(list(p2, p1))
      }
    }
  }
  
  # Return most likely parent allele using child allele
  # NA if there are two options (one contraction one expansion but same distance)
  get_original_parental_allele <- function(parent, child_allele) {
    if (abs(parent[1] - child_allele) == 0 & abs(parent[2] - child_allele) == 0) { # Stable either way
      return(parent[1])
    }
    if (abs(parent[1] - child_allele) == abs(parent[2] - child_allele)) { # Both alleles same distance away
      if (parent[1] == parent[2]) { # Parental alleles are the same, can select this as the original allele
        return(parent[1])
      } else {
        return(NA) # Parental alleles are different but same length distance, cannot tell
      }
    }
    else if (abs(parent[1] - child_allele) < abs(parent[2] - child_allele)) { # Choose allele with smallest distance
      return(parent[1])
    } else {
      return(parent[2])
    }
  }
  
  # X chromosome
  if (length(p2) == 1) {
    if (length(ch) == 1) { # Son
      if (ch %in% p1) return(list(c(ch, ch)))
      if ((abs(p1[1] - ch) == abs(p1[2] - ch)) & (p1[1] != p1[2])) return(list(NA))
      if (abs(p1[1] - ch) <= abs(p1[2] - ch)) return(c(p1[1], ch)) else return(list(c(p1[2], ch)))
    } else { # Daughter
      if (ch[1] %in% p1 & ch[2] %in% p2) return(list(c(ch[1], ch[1]), c(ch[2], ch[2])))
      else if (ch[2] %in% p1 & ch[1] %in% p2) return(list(c(ch[1], ch[1]), c(ch[2], ch[2])))
      
      # Determine optimal arrangement of child alleles to parent
      p1_p2_mins <- c(min(abs(p1[1] - ch[1]), abs(p1[2] - ch[1])), abs(p2[1] - ch[2]))
      p2_p1_mins <- c(abs(p2[1] - ch[1]), min(abs(p1[1] - ch[2]), abs(p1[2] - ch[2])))
      order <- get_optimal_order(p1_p2_mins, p2_p1_mins, ch)
      if (all(is.na(order))) return(list(NA, NA))
      
      # Determine optimal allele from mother
      maternal_transmission <- NA
      if (length(order[[1]]) == 2) {
        maternal_transmission <- c(get_original_parental_allele(order[[1]], ch[1]), ch[1])
        if (is.na(maternal_transmission[1])) return(list(NA, NA))
        return(list(maternal_transmission, c(p2[1], ch[2])))
      } else {
        maternal_transmission <- c(get_original_parental_allele(order[[2]], ch[2]), ch[2])
        if (is.na(maternal_transmission[1])) return(list(NA, NA))
        return(list(maternal_transmission, c(p2[1], ch[1])))
      }
    }
  }
  
  # No mutations
  if ((ch[1] %in% p1 & ch[2] %in% p2)) {
    return(list(c(ch[1], ch[1]), c(ch[2], ch[2])))
  } else if (ch[1] %in% p2 & ch[2] %in% p1) {
    return(list(c(ch[2], ch[2]), c(ch[1], ch[1])))
  }
  
  # Calculate all absolute differences between child alleles and all parental alleles
  diffs <- matrix(NA, nrow = 2, ncol = 4)  # Initialize a matrix to store differences
  
  # Differences for c[1] and c[2] compared to all parental alleles
  diffs[1, 1] <- abs(ch[1] - p1[1])
  diffs[1, 2] <- abs(ch[1] - p1[2])
  diffs[1, 3] <- abs(ch[1] - p2[1])
  diffs[1, 4] <- abs(ch[1] - p2[2])
  diffs[2, 1] <- abs(ch[2] - p1[1])
  diffs[2, 2] <- abs(ch[2] - p1[2])
  diffs[2, 3] <- abs(ch[2] - p2[1])
  diffs[2, 4] <- abs(ch[2] - p2[2])
  
  # Assign child's alleles to parents based on minimum distance to parental alleles
  p1_p2_mins <- c(min(diffs[1,1], diffs[1,2]), min(diffs[2,3], diffs[2,4])) # c1 from p1 and c2 from p2
  p2_p1_mins <- c(min(diffs[1,3], diffs[1,4]), min(diffs[2,1], diffs[2,2])) # c1 from p2 and c2 from p1
  order <- get_optimal_order(p1_p2_mins, p2_p1_mins, ch)
  if (all(is.na(order))) return(list(NA, NA))
  
  # Find original parent alleles in each individual (smallest distance)
  # Edge case if can't tell (same distance in opposite directions)
  c1_transmission <- c(get_original_parental_allele(order[[1]], ch[1]), ch[1])
  c2_transmission <- c(get_original_parental_allele(order[[2]], ch[2]), ch[2])
  
  # Reject both if either could not be determined
  if (any(is.na(c1_transmission[1]), is.na(c2_transmission[1]))) return(list(NA, NA))
  
  if (all(order[[1]] == p1)) {
    return(list(c1_transmission, c2_transmission))
  } else {
    return(list(c2_transmission, c1_transmission))
  }
}

## Find number of trios with germline transmissions
edge_cases <- list()
transmission_rows <- apply(children, 1, function (row) {
  father <- df_filt[df_filt$SampleID == row[["FatherID"]],]
  mother <- df_filt[df_filt$SampleID == row[["MotherID"]],]
  
  temp_tr_rows <- list()
  for (locus in disease_loci) {
    child_a <- row[[locus]]
    father_a <- father[[locus]][[1]]
    mother_a <- mother[[locus]][[1]]
    ru_length <- info_df[info_df$gene == locus,]$repeatunitlen
    ref_copies <- info_df[info_df$gene == locus,]$reference_length
    min_intermediate <- info_df[info_df$gene == locus,]$intermediate_min
    min_pathogenic <- info_df[info_df$gene == locus,]$pathogenic_min
  
    # Find which parental allele became which child allele
    if (0 %in% mother_a | 0 %in% father_a | 0 %in% child_a | 
        any(is.na(mother_a)) | any(is.na(father_a)) | any(is.na(child_a))) {
      next
    }
    transmissions <- return_transmissions(mother_a, father_a, child_a)
    if (any(is.na(transmissions))) {
      edge_cases[[length(edge_cases) + 1]] <<- data.frame(ChildID = row[["SampleID"]], gene = locus)
    }
    
    # Add each transmission as a new row
    if (length(transmissions[[1]]) == 2) {
      
      p <- NA
      i <- NA
      x <- NA
      y <- NA
      if (!is.na(min_pathogenic)) {
        if (transmissions[[1]][1] < min_pathogenic & transmissions[[1]][2] >= min_pathogenic) {
          p <- T
        }
        if (transmissions[[1]][1] >= min_pathogenic & transmissions[[1]][2] < min_pathogenic) {
          x <- T
        }
      }
      if (!is.na(min_intermediate)) {
        if (transmissions[[1]][1] < min_intermediate & transmissions[[1]][2] >= min_intermediate) {
          i <- T
        }
        if (transmissions[[1]][1] >= min_intermediate & transmissions[[1]][2] < min_intermediate) {
          y <- T
        }
      }
      
      temp_tr_rows[[length(temp_tr_rows) + 1]] <- data.frame(Population = row[["Population"]], 
                                                             FamilyID = row[["FamilyID"]], 
                                                             ParentID = row[["MotherID"]],
                                                             OtherParentID = row[["FatherID"]],
                                                             ChildID = row[["SampleID"]], 
                                                             Gene = locus, ParentSex = 1,
                                                             MotifLength = ru_length,
                                                             RefCopies = ref_copies,
                                                             ParentAllele = transmissions[[1]][1], 
                                                             ChildAllele = transmissions[[1]][2],
                                                             Magnitude = transmissions[[1]][2] - transmissions[[1]][1], 
                                                             ParentBpLength = transmissions[[1]][1] * ru_length,
                                                             ChildBpLength = transmissions[[1]][2] * ru_length,
                                                             BecamePathogenic = p,
                                                             BecameIntermediate = i,
                                                             path_to_ben = x,
                                                             int_to_ben = y)
    }
    if (length(transmissions) == 2) {
      if (length(transmissions[[2]]) == 2) {
        p <- NA
        i <- NA
        x <- NA
        y <- NA
        if (!is.na(min_pathogenic)) {
          if (transmissions[[2]][1] < min_pathogenic & transmissions[[2]][2] >= min_pathogenic) {
            p <- T
          }
          if (transmissions[[2]][1] >= min_pathogenic & transmissions[[2]][2] < min_pathogenic) {
            print(locus)
            print(transmissions[[2]])
            print(min_pathogenic)
            x <- T
          }
        }
        if (!is.na(min_intermediate)) {
          if (transmissions[[2]][1] < min_intermediate & transmissions[[2]][2] >= min_intermediate) {
            i <- T
          }
          if (transmissions[[2]][1] >= min_intermediate & transmissions[[2]][2] < min_intermediate) {
            y <- T
          }
        }
        
        temp_tr_rows[[length(temp_tr_rows) + 1]] <- data.frame(Population = row[["Population"]], 
                                                               FamilyID = row[["FamilyID"]], 
                                                               ParentID = row[["FatherID"]],
                                                               OtherParentID = row[["MotherID"]],
                                                               ChildID = row[["SampleID"]], 
                                                               Gene = locus, ParentSex = 0,
                                                               MotifLength = ru_length,
                                                               RefCopies = ref_copies,
                                                               ParentAllele = transmissions[[2]][1], 
                                                               ChildAllele = transmissions[[2]][2],
                                                               Magnitude = transmissions[[2]][2] - transmissions[[2]][1],
                                                               ParentBpLength = transmissions[[2]][1] * ru_length,
                                                               ChildBpLength = transmissions[[2]][2] * ru_length,
                                                               BecamePathogenic = p,
                                                               BecameIntermediate = i,
                                                               path_to_ben = x,
                                                               int_to_ben = y)
      }
    }
  }
  return(temp_tr_rows)
})
transmissions_df <- do.call(rbind, unlist(transmission_rows, recursive = F))
rm(transmission_rows)
edge_cases_df <- do.call(rbind, edge_cases)

#------------------------------- Pre QC analysis -------------------------------

# Mendelian concordance (97.4%) (20 100% concordant across all trios)
nrow(transmissions_df[transmissions_df$Magnitude == 0,]) / nrow(transmissions_df)
concordant_100 <- transmissions_df %>%
  group_by(Gene) %>%
  filter(all(Magnitude == 0)) %>%
  summarise(n = n_distinct(Gene)) %>%
  select(Gene) %>%
  unlist(use.names = F)

# Reference length vs discordance rate
reflength_discordance_df <- transmissions_df %>%
  group_by(Gene) %>%
  summarise(
    RefCopies = mean(RefCopies),
    RefLength = mean(RefCopies) * mean(MotifLength),
    MeanParentRU = mean(ParentAllele),
    MeanParentLength = mean(ParentBpLength),
    DiscordanceRate = sum(Magnitude != 0) / n()) %>%
  mutate(Inclusion = ifelse(Gene %in% c("XYLT1", "STARD7", "FGF14", "YEATS2", "RFC1"), "Excluded", "Included"))
ggplot(reflength_discordance_df, aes(x = MeanParentLength, y = DiscordanceRate, label = Gene, colour = Inclusion)) +
  geom_point(aes(fill = Inclusion)) +
  labs(x = "Mean Parental Allele Length (bp)", y = "Discordance Rate (%)") +
  theme_classic() +
  scale_y_continuous(labels = percent_format()) +
  geom_text_repel() +
  scale_color_manual(values = c("Excluded" = "red", "Included" = "black"))
ggsave_wrap("reflengthbp_scatter.png", 6, 5)

#------------------------------- QC filtering ----------------------------------

# # Download pileups for all unstable transmissions (not XYLT1)
# transmissions_df %>%
#   filter(Magnitude != 0, Gene != "XYLT1" ) %>%
#   create_pileup_table(filename="pileups_preQC_noXYLT1.csv")
# 
# Download pileups for all stable TCF4 over 25RU
transmissions_df %>%
  filter(Magnitude == 0, Gene == "TCF4", ParentAllele > 25) %>%
  create_pileup_table(filename="pileups_preQC_TCF4_stable.csv")
# 
# # Download specific FMR1 trio pileup for miscall example
# transmissions_df %>%
#   filter(ParentID == "HG00577", ChildID == "HG00579", Gene == "FMR1") %>%
#   create_pileup_table(filename="pileups_FMR1_miscall_eg.csv")

# Remove transmissions where parent allele is longer than read length
tr_df_filt <- transmissions_df[transmissions_df$ParentBpLength <= 150,]

# Remove XYLT1
tr_df_filt <- tr_df_filt[!(tr_df_filt$Gene %in% c('XYLT1')),]

# Download pileups for unstable transmissions involving >100bp alleles
tr_df_filt %>%
  filter((Magnitude != 0) & ((ChildBpLength >= 100) | (ParentBpLength >=100))) %>%
  create_pileup_table(filename="pileups_for_download.csv")

# Remove transmissions involving individuals that didn't pass the visual QC
inclusion_list <- read.csv("../pileups/transmission_qc_results_passed.csv", header = F, col.names = c("ParentID", "ChildID", "Gene"))
over_100_passed <- semi_join(tr_df_filt, inclusion_list, by = c("ParentID", "ChildID", "Gene"))
tr_df_filt <- rbind(tr_df_filt %>% filter(!((Magnitude != 0) & ((ChildBpLength >= 100) | (ParentBpLength >=100)))), over_100_passed)

# Remove bad quality loci
tr_df_filt <- tr_df_filt[!(tr_df_filt$Gene %in% c('STARD7', "FGF14", "YEATS2", "RFC1")),]

# NOTCH2NLC QC
create_pileup_table(tr_df_filt %>% filter(Gene == "NOTCH2NLC" & Magnitude != 0), "pileups_for_download_NOTCH2NLC.csv")
NOTCH2NLC_failed <- read.csv("../pileups/NOTCH2NLC_failed.csv", header = F, col.names = c("ParentID", "ChildID", "Note")) %>% 
  dplyr::select(c("ParentID", "ChildID")) %>% 
  mutate(Gene = "NOTCH2NLC")

tr_final <- anti_join(tr_df_filt, NOTCH2NLC_failed, by = join_by(ParentID, ChildID, Gene))

# # Download pileups for all unstable transmissions after QC
# tr_final %>%
#   filter(Magnitude != 0) %>%
#   create_pileup_table(filename="pileups_postQC.csv")

# ------------------------------ Summaries -------------------------------------

# Mendelian concordance
nrow(tr_final[tr_final$Magnitude == 0,]) / nrow(tr_final)
tr_final %>%
  group_by(Gene) %>%
  filter(all(Magnitude == 0)) %>%
  summarise(n = n_distinct(Gene))

# Average number of de novo STR mutations per sample
tr_final %>%
  group_by(ChildID) %>%                        
  summarise(non_zero_count = sum(Magnitude != 0)) %>%
  summarise(avg_non_zero = mean(non_zero_count)) 

# Parental inheritance (PATERNAL)
nrow(tr_final %>% filter(Magnitude != 0 & ParentSex == 0)) / nrow(tr_final %>% filter(Magnitude != 0))
nrow(tr_final %>% filter(ParentSex == 0))

# Estimated mutation rate across all loci
nrow(tr_final %>% filter(Magnitude != 0)) * 100 / nrow(tr_final)

nrow(tr_final %>% filter(!(Gene %in% x_chromosome_loci), ParentSex == 0))
nrow(tr_final %>% filter(!(Gene %in% x_chromosome_loci), ParentSex == 1))

# # Magnitude changes across motif lengths
# tr_final %>%
#   mutate(Magnitude = case_when(
#     Magnitude == 0 ~ "Consistent",
#     abs(Magnitude) == 1 ~ "Off by one",
#     abs(Magnitude) > 1  ~ "Larger error"
#   )) %>%
#   group_by(MotifLength, Magnitude) %>%
#   summarise(count = n()) %>%
#   mutate(percentage = count / sum(count) * 100) %>%
#   ungroup() %>%
#   ggplot(aes(x = factor(MotifLength), y = percentage, fill = Magnitude)) +
#   geom_bar(stat = "identity") + 
#   labs(x = "Motif Length", y = "Frequency", fill = "Error Type") +
#   theme_classic() + 
#   scale_fill_manual(values = c("Consistent" = "green4", "Off by one" = "orange", "Larger error" = "red3"))

# --------------------Benign -> Pathogenic analysis-----------------------------
tr_final %>% filter(!is.na(BecamePathogenic))
tr_final %>% filter(!is.na(BecameIntermediate))
tr_final %>% filter(!is.na(path_to_ben))
tr_final %>% filter(!is.na(int_to_ben))

create_pileup_table(pathogenic_new, "pileups_for_download_newpathogenic.csv")

#------------------------------- Plotting ---------------------------------------------------------#

# ---------------------- INITIAL SUMMARY-- ----------------------------------- #
tr_final %>%
  group_by(MotifLength) %>%
  summarise(
    non_zero_count = sum(Magnitude != 0),
    total_count = n(),
    perc_non_zero = (non_zero_count / total_count)
  ) %>%
  ggplot(aes(x=MotifLength, y=perc_non_zero)) +
  geom_col() +
  theme_classic() +
  labs(x = "Motif Length (bp)", y = "Discordance Rate (%)") +
  scale_y_continuous(labels = percent_format()) +
  geom_text(aes(label=paste(round(perc_non_zero*100, 3), "%", sep=""), y=perc_non_zero), vjust=-0.8)
ggsave_wrap("discordance_motif_barplot.png", 6, 5)

# # # Paternal vs maternal transmissions
# patmatexp <- tr_final %>%
#   filter(Magnitude != 0) %>%
#   group_by(ParentSex) %>%
#   summarise(count = n(), .groups = 'drop') %>%
#   mutate(percentage = count / sum(count) * 100)
# ggplot(patmatexp, aes(x = "", y = percentage, fill = factor(ParentSex, labels = c("Paternal", "Maternal")))) +
#   geom_bar(stat = "identity", position = "stack") +
#   labs(x = NULL, y = "Percentage of STR mutations", fill = "Parent of origin") +
#   scale_fill_manual(values = c("Paternal" = "lightblue", "Maternal" = "pink")) +
#   scale_y_continuous(labels = scales::percent_format(scale = 1)) +
#   theme_minimal() +
#   theme(legend.position = "right") +
#   geom_text(aes(label = round(percentage, 2)), position = position_stack(vjust = 0.5), size = 4)


# ---------------------- OVERALL (SUMMARY) ----------------------------------- #
# Magnitude of jumps (all)
ggplot(tr_final[tr_final$Magnitude != 0,], aes(x=Magnitude)) +
  geom_histogram(binwidth = 1, aes(y=after_stat(count / sum(count))), col="white") +
  theme_classic() +
  ylab("Percentage of Unstable Transmissions (%)") +
  xlab("Magnitude of Jump (RU)") +
  scale_y_continuous(labels = percent_format()) +
  geom_vline(xintercept = 0, colour = "red", linewidth = 1, linetype = 2)
ggsave_wrap("jump_histogram.png", 7, 5)

# Examining outliers of above
magnitude_outliers <- tr_final[tr_final$Magnitude < -10,]
create_pileup_table(magnitude_outliers, "pileups_for_download_magnitude_outliers.csv") # Download pileups for visual examination

# Magnitude of jumps (faceted by motif length)
ggplot(tr_final[tr_final$Magnitude != 0,] %>% 
         mutate(MotifLength = ifelse(MotifLength == 3, "3", ">3")) %>%
         mutate(MotifLength = factor(MotifLength, levels = c("3", ">3"))), aes(x=Magnitude)) +
  geom_histogram(binwidth = 1, aes(y=after_stat(count / tapply(..count.., ..PANEL.., sum)[..PANEL..])), col="white") +
  theme_classic() +
  ylab("Percentage of Unstable Transmissions (%)") +
  xlab("Magnitude of Jump (RU)") +
  geom_vline(xintercept = 0, colour = "red", linewidth = 1, linetype = 2) +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~MotifLength, nrow = 1)
ggsave_wrap("jump_histogram_faceted.png", 11, 5)

# Magnitude of jumps (faceted by motif length) (BEFORE QC)
ggplot(, aes(x=Magnitude)) +
  geom_histogram(data = transmissions_df[transmissions_df$Magnitude != 0,] %>%
                   filter(Gene != "XYLT1") %>%
                   mutate(MotifLength = ifelse(MotifLength == 3, "Trinucleotide Motifs", "Motifs > 3bp")) %>%
                   mutate(MotifLength = factor(MotifLength, levels = c("Trinucleotide Motifs", "Motifs > 3bp"))), binwidth = 1, alpha=0.5, position = "identity", aes(y=after_stat(count)), col="white", fill = "#615cf7") +
  geom_histogram(data = tr_final[tr_final$Magnitude != 0,] %>% 
                   mutate(MotifLength = ifelse(MotifLength == 3, "Trinucleotide Motifs", "Motifs > 3bp")) %>%
                   mutate(MotifLength = factor(MotifLength, levels = c("Trinucleotide Motifs", "Motifs > 3bp"))), binwidth = 1, alpha=0.5, position = "identity", aes(y=after_stat(count)), col="white", fill = "green") +
  theme_classic() +
  ylab("Count") +
  xlab("Magnitude of Jump (RU)") +
  geom_vline(xintercept = 0, colour = "red", linewidth = 0.5, linetype = 2) +
  facet_wrap(~MotifLength, nrow = 1) 
ggsave_wrap("jump_histogram_preQC_XYLT1_faceted.png", 9, 4)

# Magnitude numbers (all)
tr_final %>%
  filter(Magnitude != 0) %>%
  group_by(Magnitude) %>%
  summarise(MagnitudeN = n()) %>%
  mutate(MagnitudeChange = cut(
    abs(Magnitude),
    breaks = c(-Inf, 1, 2, 3, 4, Inf),
    labels = c("±1", "±2", "±3", "±4", ">±4"),
    right = TRUE )) %>%
  aggregate(MagnitudeN ~ MagnitudeChange, sum) %>%
  mutate(percentage = round(MagnitudeN / nrow(tr_final[tr_final$Magnitude != 0,]), 2) * 100) %>%
  rename("n" = "MagnitudeN")

tr_final %>%
  filter(Magnitude != 0) %>%
  mutate(AbsMagnitude = abs(Magnitude),
         AbsMagnitudeGroup = cut(
           AbsMagnitude,
           breaks = c(-Inf, 1, 2, 3, 4, Inf),
           labels = c("±1", "±2", "±3", "±4", ">±4"),
           right = TRUE)) %>%
  group_by(AbsMagnitudeGroup) %>%
  summarise(
    expansion_perc = sum(Magnitude > 0) / n() * 100
  )

print("stable n:")
print(nrow(tr_final %>% filter(Magnitude == 0)))

# ---------------------- INSTABILITY (STABLE OR MUTATED)---------------------- #

### STAT TESTS
logreg_instability <- tr_final %>%
  mutate(MagnitudeType = ifelse(Magnitude == 0, "Stable", "Mutated")) %>%
  mutate(MagnitudeType = factor(MagnitudeType))
logreg_instability$MagnitudeType = relevel(logreg_instability$MagnitudeType, ref = "Stable")
logreg_instability$Population <- relevel(as.factor(logreg_instability$Population), ref = "CEU")


## Logistic regression as a function of motif length/ parent allele
# RU
logistic_model <- glm(MagnitudeType ~ MotifLength + ParentAllele,
                      family = binomial, data = logreg_instability)
summary(logistic_model)
exp(coef(logistic_model))["ParentAllele"]
exp(coef(logistic_model))["MotifLength"]
# BP
logistic_model <- glm(MagnitudeType ~ MotifLength + ParentBpLength,
                      family = binomial, data = logreg_instability)
summary(logistic_model)
exp(coef(logistic_model))["ParentBpLength"]
# RU (>3bp motif)
logistic_model <- glm(MagnitudeType ~ MotifLength + ParentAllele,
                      family = binomial, data = logreg_instability %>% filter(MotifLength > 3))
summary(logistic_model)
exp(coef(logistic_model))["ParentAllele"]
# BP (>3bp motif)
logistic_model <- glm(MagnitudeType ~ ParentBpLength + MotifLength,
                      family = binomial, data = logreg_instability %>% filter(MotifLength > 3))
summary(logistic_model)
exp(coef(logistic_model))["ParentBpLength"]

## Chisqr for parental sex and instability
tab <- table(logreg_instability$ParentSex, logreg_instability$MagnitudeType)
print(tab)
chisq.test(tab)
cramers_v(tab)
#male : 0.004017184623
#female : 0.002926153135

## Logistic regression for population
logistic_model <- glm(MagnitudeType ~ Population + ParentAllele,
                      family = binomial, data = logreg_instability)
summary(logistic_model)
exp(coef(logistic_model))["PopulationPEL"]
exp(coef(logistic_model))["PopulationBEB"]

# Repeat length vs discordance rate (beta regression)
parentA_discordance_df <- tr_final %>%
  group_by(ParentAllele) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4) + 0.00001) %>%
  filter(total >= 100)
model <- betareg(discordance_rate ~ ParentAllele, data = parentA_discordance_df, link = "log")
summary(model)
coef(model)[1]
coef(model)[2]
shapiro.test(residuals(model)) # NORMAL!
exp(coef(model)["ParentAllele"])
parentA_discordance_df$fitted_values <- predict(model, type = "response")
ggplot(parentA_discordance_df, aes(x=ParentAllele, y = discordance_rate)) +
  geom_point() +
  geom_line(aes(y = fitted_values), color = "purple") +
  theme_classic() +
  ylab("Discordance Rate (%)") +
  xlab("Parental Repeat Length (RU)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.03))
ggsave_wrap("discordance_parent_ru_scatter.png", 5.5, 5)

# Total BP length vs discordance rate (beta regression)
parentbp_discordance_df <- tr_final %>%
  group_by(ParentBpLength) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4) + 0.00001) %>%
  filter(total >= 100)
model <- betareg(discordance_rate ~ ParentBpLength, data = parentbp_discordance_df, link = "log")
summary(model)
coef(model)[1]
coef(model)[2]
shapiro.test(residuals(model))
exp(coef(model)["ParentBpLength"])
parentbp_discordance_df$fitted_values <- predict(model, type = "response")
ggplot(parentbp_discordance_df, aes(x=ParentBpLength, y = discordance_rate)) +
  geom_point() +
  geom_line(aes(y = fitted_values), color = "orange3") +
  theme_classic() +
  ylab("Discordance Rate %)") +
  xlab("Parental Allele Length (bp)") +
  scale_y_continuous(labels = percent_format())
ggsave_wrap("discordance_parent_bp_scatter.png", 5.5, 5)

# Repeat length vs discordance rate (motif length >3 only)
parentA_no3_discordance_df <- tr_final %>%
  filter(MotifLength > 3) %>%
  group_by(ParentAllele) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4) + 0.00001) %>%
  filter(total >= 50)
model <- betareg(discordance_rate ~ ParentAllele, data = parentA_no3_discordance_df, link = "log")
summary(model)
coef(model)[1]
coef(model)[2]
shapiro.test(residuals(model)) # NOT NORMAL!
exp(coef(model)["ParentAllele"])
parentA_no3_discordance_df$fitted_values <- predict(model, type = "response")
ggplot(parentA_no3_discordance_df, aes(x=ParentAllele, y = discordance_rate)) +
  geom_point() +
  geom_line(aes(y = fitted_values), color = "darkgreen") +
  theme_classic() +
  ylab("Discordance Rate (%)") +
  xlab("Parental Repeat Length (RU)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.03))
ggsave_wrap("discordance_no3_parent_ru_scatter.png", 5.5, 5)

# bp length vs discordance rate (motif length >3 only)
parentbp_no3_discordance_df <- tr_final %>%
  filter(MotifLength > 3) %>%
  group_by(ParentBpLength) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4) + 0.00001) %>%
  filter(total >= 50)
model <- betareg(discordance_rate ~ ParentBpLength, data = parentbp_no3_discordance_df, link = "log")
summary(model)
coef(model)[1]
coef(model)[2]
shapiro.test(residuals(model)) # not NORMAL!
exp(coef(model)["ParentBpLength"])
parentbp_no3_discordance_df$fitted_values <- predict(model, type = "response")
ggplot(parentbp_no3_discordance_df, aes(x=ParentBpLength, y = discordance_rate)) +
  geom_point() +
  geom_line(aes(y = fitted_values), color = "blue") +
  theme_classic() +
  ylab("Discordance Rate (%)") +
  xlab("Parental Allele Length (bp)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 0.03))
ggsave_wrap("discordance_no3_parent_bp_scatter.png", 5.5, 5)

# # Number of expansions and contractions
# totals <- tr_final %>%
#   filter(Magnitude != 0) %>%
#   mutate(Sign = ifelse(Magnitude < 0, "Contraction", "Expansion")) %>%
#   group_by(Sign) %>%
#   summarise(total = n(), .groups = 'drop')
# tr_final %>%
#   filter(Magnitude != 0) %>%
#   mutate(Sign = ifelse(Magnitude < 0, "Contraction", "Expansion"),
#          MagnitudeCat = factor(case_when(
#            abs(Magnitude) == 1 ~ "1",
#            abs(Magnitude) == 2 ~ "2",
#            abs(Magnitude) == 3 ~ "3",
#            abs(Magnitude) == 4 ~ "4",
#            abs(Magnitude) > 4 ~ ">4"
#          ), levels = c("1", "2", "3", "4", ">4"))) %>%
#   group_by(Sign, MagnitudeCat) %>%
#   summarise(SignN = n(), .groups = 'drop') %>%
#   ggplot(aes(x = Sign, y = SignN, fill = MagnitudeCat)) +
#   geom_col() +
#   theme_classic() +
#   ylab("Count") +
#   xlab("Type") +
#   geom_text(aes(label = SignN), 
#             position = position_stack(vjust = 0.5), size = 4) +
#   geom_text(data = totals, aes(x = Sign, y = total + 1, label = paste("total = ", total, sep="")), 
#             inherit.aes = F, vjust = -0.5, size = 4) +
#   scale_fill_brewer(palette = "Pastel1") +
#   guides(fill=guide_legend(title="Magnitude"))
# ggsave_wrap("exp_contr_stacked.png", 5.5, 6)

# Parental sex and expansion vs contraction
total_vec <- c(tr_final %>% filter(ParentSex == 0, Magnitude != 0) %>% nrow(), 0,
               tr_final %>% filter(ParentSex == 1, Magnitude != 0) %>% nrow(), 0)
ratio_vec <- c(tr_final %>% filter(ParentSex == 0, Magnitude > 0) %>% nrow() / tr_final %>% filter(ParentSex == 0, Magnitude < 0) %>% nrow(), 0,
               tr_final %>% filter(ParentSex == 1, Magnitude > 0) %>% nrow() / tr_final %>% filter(ParentSex == 1, Magnitude < 0) %>% nrow(), 0)
ratio_vec <- round(ratio_vec, 2)
tr_final %>%
  group_by(ParentSex) %>%
  mutate(ParentSex = dplyr::recode(ParentSex, '0' = 'Paternal', '1' = 'Maternal')) %>%
  summarise(Expansions = sum(Magnitude > 0), Contractions = sum(Magnitude < 0)) %>%
  mutate(Ratio = round(Expansions / Contractions, 2))
patmatexp <- tr_final %>%
  filter(Magnitude != 0) %>%
  mutate(MutationType = ifelse(Magnitude > 0, "Expansion", "Contraction")) %>%
  group_by(ParentSex, MutationType) %>%
  summarise(count = n(), .groups = 'drop') %>%
  mutate(total = total_vec, ratio = ratio_vec)
ggplot(patmatexp, aes(x = factor(ParentSex, labels = c("Paternal", "Maternal")), 
                      y = count, 
                      fill = MutationType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Parent of Origin", y = "Number of Unstable Transmissions", fill = "Direction") +
  scale_fill_manual(values = c("Expansion" = "#fab5ae", "Contraction" = "#b3cce2")) +
  theme_classic() +
  theme(legend.position = "right") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 4) +
  geom_text(aes(label = ifelse(total > 0, paste(paste("total = ", total, sep=""), "\n", paste("ratio = ", ratio, sep="")), "")), 
            position = position_stack(vjust = 1.2), size = 4)
ggsave_wrap("parent_exp_contr_barplot.png", 5.25, 5.25)


# ---------------------- EXPANSION OR CONTRACTION ---------------------------- #

# Logistic regression for probability of expansion vs contraction as function of motiflength/parentallele
logreg_direction <- tr_final %>%
  filter(Magnitude != 0,) %>%
  mutate(MagnitudeType = ifelse(Magnitude > 0, "Expansion", "Contraction")) %>%
  mutate(MagnitudeType = factor(MagnitudeType))
logreg_direction$MagnitudeType = relevel(logreg_direction$MagnitudeType, ref = "Contraction")
logistic_model <- glm(MagnitudeType ~ MotifLength + ParentAllele, 
                      family = binomial, data = logreg_direction)
summary(logistic_model)
exp(coef(logistic_model)["MotifLength"])
exp(coef(logistic_model)["ParentAllele"])

# Chisqr for expansion vs contraction diff
mutation_table <- table(logreg_direction$MagnitudeType)
chisq_test <- chisq.test(mutation_table)
chisq_test

## Chisqr for parental sex and direction
tab <- table(logreg_direction$ParentSex, logreg_direction$MagnitudeType)
chisq.test(tab)
cramers_v(tab)

## Chisqr for population and direction
tab <- table(logreg_direction$Population, logreg_direction$MagnitudeType)
chisq.test(tab, simulate.p.value = T)
cramers_v(tab)

# # Expansion percentage with parental repeat length
# parentRU_perc_df <- tr_final %>%
#   filter(Magnitude != 0) %>%
#   group_by(ParentAllele) %>%
#   summarise(
#     total = sum(Magnitude != 0),
#     total_expansion = sum(Magnitude > 0),
#     perc_expansion = total_expansion / total,
#     se = sqrt((perc_expansion * (1 - perc_expansion)) / total),
#     ci_lower = perc_expansion - 1.96 * se,  # Lower bound of 95% CI
#     ci_upper = perc_expansion + 1.96 * se   # Upper bound of 95% CI
#   ) %>%
#   filter(total >= 5)
# ggplot(parentRU_perc_df, aes(x = ParentAllele, y = perc_expansion)) +
#   geom_point() +
#   # geom_point(size = 3, color = "hotpink") + 
#   # geom_line(color = "hotpink", size = 1) +
#   geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "black", size = 0.5, alpha = 0.5) +
#   theme_classic() +
#   geom_smooth(method = "loess", 
#               se = TRUE,              
#               color = "hotpink",  
#               fill = "pink",
#               span = 2) + 
#   ylab("Percentage of expansions") +
#   xlab("Parental repeat length (RU)") +
#   scale_y_continuous(labels = percent_format()) +
#   coord_cartesian(ylim = c(0, 1))
# ggsave_wrap("expansion_parentRU_line.png", 8, 5)

parentA_perc_df <- tr_final %>%
  filter(Magnitude != 0, ParentAllele > 1, MotifLength == 3) %>%
  group_by(ParentAllele) %>%
  summarise(
    total = sum(Magnitude != 0),
    total_expansion = sum(Magnitude > 0),
    total_contraction = sum(Magnitude < 0),
    perc_expansion = total_expansion / total,
    perc_contraction = total_contraction / total
  ) %>%
  filter(total >= 1) %>%
  pivot_longer(cols = c(perc_expansion, perc_contraction),
               names_to = "MutationType",
               values_to = "Percentage") %>%
  mutate(MutationType = factor(MutationType, levels = c("perc_contraction", "perc_expansion")))
ggplot(parentA_perc_df, aes(x = ParentAllele, y = Percentage, fill = MutationType)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("perc_contraction" = "#b3cce2", "perc_expansion" = "#fab5ae"),
                    labels = c("Contraction", "Expansion")) +
  ylab("Percentage of Unstable Transmissions (%)") +
  xlab("Parental Repeat Length (RU)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(fill = "Direction") +
  xlim(c(7,38))
ggsave_wrap("expansion_parentA_stacked_trinucleotide.png", 8, 5)

# --------------------- DEGREE OF MUTATION ----------------------------------- #

# Unstable transmissions only with absolute value magnitude
degree_df <- tr_final %>%
  filter(Magnitude != 0) %>%
  mutate(AbsMagnitude = abs(Magnitude),
         AbsBpMagnitude = abs(Magnitude) * MotifLength,
         MutationType = ifelse(Magnitude > 0, "Expansion", "Contraction"))
parental_allele_groups <- unique(degree_df$ParentAllele)

shapiro.test(degree_df$AbsMagnitude) # Not normal

# Parental STR length
model <- lm(AbsMagnitude ~ ParentAllele + MotifLength, data = degree_df)
summary(model)
model <- lm(AbsMagnitude ~ ParentBpLength, data = degree_df)
summary(model)
model <- lm(AbsMagnitude ~ ParentBpLength + MotifLength, data = degree_df) # Motif length
summary(model)

# Table by motif length
degree_df %>%
  group_by(MotifLength) %>%
  summarise(median_magnitude = median(AbsMagnitude),
         mean_magnitude = mean(AbsMagnitude),
         sd_magnitude = sd(AbsMagnitude))

# Mean magnitude by parental repeat length
mean_data <- degree_df %>%
  filter(MotifLength == 3) %>%
  group_by(ParentAllele) %>%
  summarise(MeanMagnitude = mean(AbsMagnitude, na.rm = TRUE),
            total = n()) %>%
  filter(total > 6)
ggplot(mean_data, aes(x = ParentAllele, y = MeanMagnitude)) +
  geom_point(size = 2) +
  geom_line(size = 0.8, alpha = 0.5) +
  labs(x = "Parental Repeat Length (RU)", 
       y = "Mean Magnitude in Trinucleotides (RU)", 
       color = "Motif Length") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave_wrap("absolute_magnitude_by_length_tri.png", 5, 5)
mean_data <- degree_df %>%
  filter(MotifLength == 3) %>%
  group_by(ParentAllele) %>%
  summarise(MeanMagnitude = mean(Magnitude, na.rm = TRUE),
            total = n()) %>%
  filter(total > 6)
ggplot(mean_data, aes(x = ParentAllele, y = MeanMagnitude)) +
  geom_point(size = 2) +
  geom_line(size = 0.8, alpha = 0.5) + 
  labs(x = "Parental Repeat Length (RU)", 
       y = "Mean Change of Length in Trinucleotides (RU)", 
       color = "Motif Length") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  geom_hline(yintercept = 0)
ggsave_wrap("magnitude_by_length_tri.png", 5, 5)

# Parent sex stratified kruskal's test
p_values <- c()
eta_squared_values <- c()
results <- lapply(parental_allele_groups, function(allele) {
  subset_data <- degree_df[degree_df$ParentAllele == allele, ]
  if (length(unique(subset_data$ParentSex)) > 1) {
    kruskal_result <- kruskal.test(AbsMagnitude ~ ParentSex, data = subset_data)
    
    H <- kruskal_result$statistic 
    k <- length(unique(subset_data$ParentSex)) 
    n <- nrow(subset_data) 
    
    eta_squared <- max((H - (k - 1)) / (n - k), 0)
    
    p_values <<- c(p_values, kruskal_result$p.value)
    eta_squared_values <<- c(eta_squared_values, eta_squared)
    return(list(p_value = kruskal_result$p.value, eta_squared = eta_squared))
  }
})
summary_p_value <- sumlog(na.omit(p_values))$p
summary_eta_squared <- mean(na.omit(eta_squared_values[!is.infinite(eta_squared_values)]), na.rm = TRUE)
cat("Summary p-value parent sex (Fisher's method):", summary_p_value, "\n")
cat("Summary eta-squared parent sex (mean):", summary_eta_squared, "\n")

# Parental biases for magnitude
ggplot(degree_df, aes(x = factor(ParentSex), y = AbsMagnitude)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("0" = "Paternal", "1" = "Maternal")) +
  labs(x = "Parent of Origin", y = "Magnitude (RU)") +
  theme_classic()
ggsave_wrap("sex_boxplot_RU.png", 5, 5)
ggplot(degree_df, aes(x = factor(ParentSex), y = AbsBpMagnitude)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("0" = "Paternal", "1" = "Maternal")) +
  labs(x = "Parent of Origin", y = "Magnitude (bp)") +
  theme_classic()
ggsave_wrap("sex_boxplot_bp.png", 5, 5)

# ---------------------- LOCI SPECIFIC --------------------------------------- #

# Locus stratified kruskal's test for magnitude
p_values <- c()
eta_squared_values <- c()
results <- lapply(parental_allele_groups, function(allele) {
  subset_data <- degree_df[degree_df$ParentAllele == allele, ]
  if (length(unique(subset_data$Gene)) > 1) {
    kruskal_result <- kruskal.test(AbsMagnitude ~ Gene, data = subset_data)
    
    H <- kruskal_result$statistic
    k <- length(unique(subset_data$Gene))
    n <- nrow(subset_data)
    
    eta_squared <- (H - (k - 1)) / (n - k)
    
    p_values <<- c(p_values, kruskal_result$p.value)
    eta_squared_values <<- c(eta_squared_values, eta_squared)
    return(list(p_value = kruskal_result$p.value, eta_squared = eta_squared))
  }
})
summary_p_value <- sumlog(na.omit(p_values))$p
summary_eta_squared <- mean(na.omit(eta_squared_values[!is.infinite(eta_squared_values)]), na.rm = TRUE)
cat("Summary p-value gene (Fisher's method):", summary_p_value, "\n")
cat("Summary eta-squared gene (mean):", summary_eta_squared, "\n")

# GLM for gene vs expansion or contraction
logistic_model <- glm(MagnitudeType ~ Gene + ParentAllele, 
                      family = binomial, data = logreg_direction)
summary(logistic_model)

# GLM for gene vs stable or unstable
logistic_model <- glm(MagnitudeType ~ Gene + ParentAllele, 
                      family = binomial, data = logreg_instability)
summary(logistic_model)

# Instability summary in all loci
mutation_data <- tr_final %>%
  group_by(Gene) %>%
  summarise(total = n(), unstable_n = sum(Magnitude != 0), 
            Contraction = sum(Magnitude < 0), Expansion = sum(Magnitude > 0)) %>%
  filter(unstable_n > 0) %>%
  pivot_longer(
    cols = c("Expansion", "Contraction"),
    names_to = "Type",
    values_to = "values"
  ) %>%
  mutate(percentage = (unstable_n / total) * 100)
mutation_data %>%
  ggplot(aes(x = Gene, y = values, fill = Type)) +
  scale_fill_manual(values = c("Expansion" = "#fab5ae", "Contraction" = "#b3cce2")) +
  geom_col() +
  theme_classic() +
  geom_text(aes(label = ifelse(values == 0, "", values)), 
            position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
  scale_y_continuous(
    name = "Total Unstable Transmissions",
    sec.axis = sec_axis(~ . / max(mutation_data$total) * 100, name = "Percentage of Total Transmissions at Locus (%)", labels = scales::label_percent(scale = 1))
  ) +
  theme(
    axis.text.x = element_text(size = 12, angle = 70, vjust = 1, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.position = c(0.1, 0.8)
  ) +
  labs(x = "Locus", fill = "Direction")
ggsave_wrap("instability_by_gene_exp_contr_barplot.png", 10, 7)
  
# Magnitude by locus
magnitude_data <- tr_final %>%
  group_by(Gene) %>%
  summarise(total = n(),
            unstable_n = sum(Magnitude != 0), 
            TotalMagnitude = sum(Magnitude, na.rm = TRUE),
            Contraction = sum(Magnitude[Magnitude < 0], na.rm = TRUE),
            Expansion = sum(Magnitude[Magnitude > 0], na.rm = TRUE)) %>% 
  filter(unstable_n > 0) %>%
  pivot_longer(
    cols = c("Expansion", "Contraction"),
    names_to = "Type",
    values_to = "MagnitudeValues"
  )
magnitude_data %>%
  ggplot(aes(x = Gene, y = MagnitudeValues, fill = Type)) +
  scale_fill_manual(values = c("Expansion" = "#fab5ae", "Contraction" = "#b3cce2")) +
  geom_col() +
  theme_classic() +
  scale_y_continuous(
    name = "Total Magnitude Change (RU)",
  ) +
  theme(
    axis.text.x = element_text(size = 12, angle = 70, vjust = 1, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.position = c(0.1, 0.2)
  ) +
  labs(x = "Locus", fill = "Direction")
ggsave_wrap("instability_by_gene_magnitude_barplot.png", 10, 7)

# Mean parental allele by locus
allele_data <- tr_final %>%
  filter(Magnitude != 0) %>%
  group_by(Gene) %>%
  summarise(AverageParentAllele = mean(ParentBpLength, na.rm = TRUE))
ggplot(allele_data, aes(x = reorder(Gene, AverageParentAllele), y = AverageParentAllele)) +
  geom_col() +
  theme_classic() +
  labs(x = "Locus", y = "Average Unstable Parental Allele Length (bp)") +
  theme(
    axis.text.x = element_text(size = 12, angle = 70, vjust = 1, hjust = 1, face = "italic"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15)
  )
ggsave_wrap("parentlength_by_gene_barplot.png", 10, 7)

# definitely check for sex biases for each gene
logistic_model <- glm(ParentSex ~ Gene + ParentSex, 
                      family = binomial, data = tr_final)
summary(logistic_model)
exp(coef(logistic_model)["MotifLength"])
exp(coef(logistic_model)["ParentAllele"])

# Expansion vs contraction in trinucleotides
parentA_perc_df <- tr_final %>%
  filter(Magnitude != 0, MotifLength == 3, ParentAllele > 1) %>%
  group_by(ParentAllele) %>%
  summarise(
    total = sum(Magnitude != 0),
    total_expansion = sum(Magnitude > 0),
    total_contraction = sum(Magnitude < 0),
    perc_expansion = total_expansion / total,
    perc_contraction = total_contraction / total
  ) %>%
  filter(total >= 6) %>%
  pivot_longer(cols = c(perc_expansion, perc_contraction),
               names_to = "MutationType",
               values_to = "Percentage") %>%
  mutate(MutationType = factor(MutationType, levels = c("perc_contraction", "perc_expansion")))
ggplot(parentA_perc_df, aes(x = ParentAllele, y = Percentage, fill = MutationType)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("perc_contraction" = "#b3cce2", "perc_expansion" = "#fab5ae"),
                    labels = c("Contraction", "Expansion")) +
  ylab("Percentage of Unstable Transmissions (%)") +
  xlab("Parental Repeat Length (RU)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(fill = "Direction")

#---------------- Alleles by gene (not transmissions) ----------------

# Df of alleles for each individual
child_ids <- unique(tr_final$ChildID)
filtered_df <- tr_final[!tr_final$ParentID %in% child_ids, ]
filtered_df_new <- filtered_df %>%
  mutate(SampleID = ParentID, AlleleLength = ParentAllele) %>%
  dplyr::select(SampleID, Gene, AlleleLength)
original_df_new <- tr_final %>%
  mutate(SampleID = ChildID, AlleleLength = ChildAllele) %>%
  dplyr::select(SampleID, Gene, AlleleLength)
alleles_df <- bind_rows(filtered_df_new, original_df_new)
alleles_df$AlleleLength <- as.numeric(alleles_df$AlleleLength)

# Summaries for genes
summarised_allele_df <- alleles_df %>%
  group_by(Gene) %>%
  summarise(AverageAlleleLength = mean(AlleleLength, na.rm = TRUE),
            MedianAlleleLength = median(AlleleLength, na.rm = TRUE),
            StdDevAlleleLength = sd(AlleleLength, na.rm = TRUE),
            Range = max(AlleleLength, na.rm = T) - min(AlleleLength, na.rm = T))

# 100% concordant analysis
summarised_allele_df %>%
  filter(Gene %in% concordant_100)
  

# Histogram comparison TRB
summary_stats("XYLT1", binwidth = 1, limit=100, sampleid_subset = complete)
ggsave_wrap("bench_xylt1_test.png", 5, 4)
trb_df %>% filter(Gene == "XYLT1") %>%
  ggplot(aes(x=AlleleRU, y=after_stat(count / sum(count)))) +
  geom_histogram(binwidth=1, colour="black", fill="#97bbdb") +
  xlab("Repeat Length (RU)") +
  ylab("Frequency (%)") +
  scale_y_continuous(labels = percent_format()) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  xlim(c(0,100))
ggsave_wrap("bench_xylt1_ref.png", 5, 4)
summary_stats("RFC1", limit=700, sampleid_subset = complete)
ggsave_wrap("bench_rfc1_test.png", 5, 4)
trb_df %>% filter(Gene == "RFC1") %>%
  ggplot(aes(x=AlleleRU, y=after_stat(count / sum(count)))) +
  geom_histogram(binwidth=1, colour="black", fill="#97bbdb") +
  xlab("Repeat Length (RU)") +
  ylab("Frequency (%)") +
  theme_classic() +
  scale_y_continuous(labels = percent_format()) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  xlim(c(0,700))
ggsave_wrap("bench_rfc1_ref.png", 5, 4)
summary_stats("FMR1", limit=60, sampleid_subset = complete)
ggsave_wrap("bench_fmr_test.png", 5, 4)
trb_df %>% filter(Gene == "FMR1") %>%
  ggplot(aes(x=AlleleRU, y=after_stat(count / sum(count)))) +
  geom_histogram(binwidth=1, colour="black", fill="#97bbdb") +
  xlab("Repeat Length (RU)") +
  ylab("Frequency (%)") +
  scale_y_continuous(labels = percent_format()) +
  theme_classic() +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12)) +
  xlim(c(0,60))
ggsave_wrap("bench_fmr_ref.png", 5, 4)

#---------------- Pentanucleotide analysis ----------------
tr_final %>%
  filter(MotifLength == 5) %>%           
  group_by(Gene) %>%                
  summarise(
    non_zero_count = sum(Magnitude != 0),
    total_count = n(),
    perc_non_zero = (non_zero_count / total_count)
  ) %>%
  ggplot(aes(x = Gene, y = perc_non_zero)) + 
  geom_col() +
  theme_classic() +
  labs(x = "Locus", y = "Discordance Rate (%)") +
  scale_y_continuous(labels = percent_format()) +
  geom_text(aes(label = paste(round(perc_non_zero * 100, 3), "%", sep = ""), y = perc_non_zero), 
            vjust = -0.8) +
  theme(
    axis.text.x = element_text(face = "italic")
  )
ggsave_wrap("discordance_penta.png", 6, 5)
