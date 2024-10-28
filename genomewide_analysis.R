library("dplyr")
library("ggplot2")
library("tidyr")
library("readr")
library("scales")
library("effectsize")
library("rstatix")
library("metap")
library("RColorBrewer")

read_len <- 150
fig_dir <- "NEW_FIGS"
info_df <- read.csv("genomewide_info.csv")
pathogenic_list <- readLines("../../../../catalogs/disease_loci.txt")

ggsave_wrap <- function(filename, width, height) {
  ggsave(paste(fig_dir, "/gwide.", filename, sep=""), width = width, height = height)
}

homopolymers <- info_df %>%
  filter(REPEATLENGTH == 1) %>%
  select(VARID) %>%
  unlist()

df <- read_csv("genome_wide_transmissions.csv")
#------------------------------- pre-QC summaries ------------------------------

print("----------- PRE QC SUMMARIES -----------")

# Still remove homopolymers first
df <- df[!(df$Gene %in% homopolymers),]

# Total loci
print("Total loci")
print(length(unique(df[["Gene"]])))

print("Total transmissions")
print(nrow(df))

print("Total unstable transmissions")
print(nrow(df %>% filter(Magnitude != 0)))

# Mendelian concordance
print("mendelian concordance (number of 100% concordant loci)")
nrow(df[df$Magnitude == 0,]) / nrow(df)
df %>%
  group_by(Gene) %>%
  filter(all(Magnitude == 0)) %>%
  summarise(n = n_distinct(Gene)) %>%
  nrow()

# Average number of de novo STR mutations per sample
print("Average DNM per sample")
df %>%
  group_by(ChildID) %>%
  summarise(non_zero_count = sum(Magnitude != 0)) %>%
  summarise(avg_non_zero = mean(non_zero_count))

# Parental inheritance (PATERNAL)
print("Paternal inheritance (%)")
nrow(df %>% filter(Magnitude != 0 & ParentSex == 0)) / nrow(df %>% filter(Magnitude != 0))

# Estimated mutation rate across all loci
print("Average mutation rate across all loci")
nrow(df %>% filter(Magnitude != 0)) * 100 / nrow(df)
df %>%
  group_by(MotifLength) %>%
  summarise(
    non_zero_count = sum(Magnitude != 0),
    total_count = n(),
    percentage_non_zero = (non_zero_count / total_count) * 100)

#------------------------------- QC filtering ----------------------------------

# Remove bad loci
df_filt <- df[!(df$Gene %in% c('XYLT1', 'FGF14', 'RFC1', 'YEATS2', 'STARD7')),]

# Remove transmissions where parent allele is longer than read length
df_filt <- df_filt[df_filt$ParentBpLength <= 150,]

# Remove alleles from visual inspection
inclusion_list <- read.csv("../../loci_passed_visual_qc/transmission_qc_results_passed.csv", header = F, col.names = c("ParentID", "ChildID", "Gene"))
over_100_passed_pathogenic <- semi_join(df_filt, inclusion_list, by = c("ParentID", "ChildID", "Gene"))
df_filt <- rbind(df_filt %>%
                        filter(!((Magnitude != 0) &
                              ((ChildBpLength >= 100) | (ParentBpLength >= 100)) &
                              (Gene %in% pathogenic_list))), over_100_passed_pathogenic)

# Remove NOTCH2NLC from visual inspection
NOTCH2NLC_failed <- read.csv("../../loci_passed_visual_qc/NOTCH2NLC_failed.csv", header = F, col.names = c("ParentID", "ChildID", "Note")) %>%
  select(c("ParentID", "ChildID")) %>%
  mutate(Gene = "NOTCH2NLC")
tr_final <- anti_join(df_filt, NOTCH2NLC_failed, by = join_by(ParentID, ChildID, Gene))
rm(df_filt, df)
gc()

# NUMBER OF UNSTABLE TRANSMISSIONS AFTER 100BP
print("NUMBER OF UNSTABLE TRANSMISSIONS AFTER 100BP")
tr_final %>% filter(Magnitude != 0, ParentBpLength > 100 | ChildBpLength > 100) %>% nrow()
# saveRDS(tr_final, "all_trans_qced.rds")
# tr_final <- readRDS("all_trans_qced.rds")

#--------------------------- missing values? -----------------------------------

only_genotypes_df <- df %>%
  select(-c(FamilyID, Population, SampleID, CpgID, FatherID, MotherID, Sex))
sum(is.na(only_genotypes_df))
rm(only_genotypes_df)

# ------------------------------ Summaries -------------------------------------
print("----------- QCED FROM HERE -----------")

# Total loci
print("Total loci")
print(length(unique(tr_final[["Gene"]])))

print("Total transmissions")
print(nrow(tr_final))

print("Total unstable transmissions")
print(nrow(tr_final %>% filter(Magnitude != 0)))

# Mendelian concordance
print("mendelian concordance + number 100% concordant")
nrow(tr_final[tr_final$Magnitude == 0,]) / nrow(tr_final)
tr_final %>%
  group_by(Gene) %>%
  filter(all(Magnitude == 0)) %>%
  summarise(n = n_distinct(Gene))

# Average number of de novo STR mutations per sample
print("average dnms per sample")
tr_final %>%
  group_by(ChildID) %>%
  summarise(non_zero_count = sum(Magnitude != 0)) %>%
  summarise(avg_non_zero = mean(non_zero_count))

# Parental inheritance (PATERNAL)
print("paternal inheritance (%)")
nrow(tr_final %>% filter(Magnitude != 0 & ParentSex == 0)) / nrow(tr_final %>% filter(Magnitude != 0))
nrow(tr_final %>% filter(ParentSex == 0))

print("mutation rate in male then female")
nrow(tr_final %>% filter(Magnitude != 0, ParentSex == 0)) / nrow(tr_final %>% filter(ParentSex == 0))
nrow(tr_final %>% filter(Magnitude != 0, ParentSex == 1)) / nrow(tr_final %>% filter(ParentSex == 1))
# Estimated mutation rate across all loci
print("mutation rate across all loci")
nrow(tr_final %>% filter(Magnitude != 0)) * 100 / nrow(tr_final)

tr_final %>%
  group_by(MotifLength) %>%
  summarise(maxlength = max(ParentAllele))

print("expansion n")
tr_final %>% filter(Magnitude > 0) %>% nrow()
print("contraction n")
tr_final %>% filter(Magnitude < 0) %>% nrow()
#------------------------------- Plotting ---------------------------------------------------------#

# ---------------------- INITIAL SUMMARY-- ----------------------------------- #

info_df %>%
  filter(VARID %in% tr_final$Gene) %>%
  group_by(REPEATLENGTH) %>%
  summarise(n = n()) %>%
  mutate(percentage = n / sum(n)) %>%
  ggplot(aes(x=REPEATLENGTH, y=percentage)) +
  geom_col() +
  theme_classic() +
  xlab("Motif Length (bp)") +
  ylab("Percentage of Total (%)") +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(aes(label = n),
            position = position_dodge(width=0.9), size = 4, vjust=-0.5)
ggsave_wrap("motif_lengths_gwide.png", width = 5, height = 5)

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

# Absolute numbers by bp and RU
# BP
print("absolute numbers by bp")
max_value <- 150
filtered_data <- tr_final %>%
  mutate(LengthBin = cut(ParentBpLength,
                         breaks = seq(1, max_value + 10, by = 10),
                         include.lowest = TRUE,
                         right = TRUE,
                         labels = paste0("[", seq(1, max_value, by = 10), ", ", seq(10, max_value + 9, by = 10), "]"))) %>%
  group_by(LengthBin) %>%
  summarise(TotalCount = n(),
            UnstableCount = sum(Magnitude != 0),
            Percentage = UnstableCount / TotalCount,
            .groups = 'drop')
ggplot(filtered_data, aes(x = LengthBin, y = UnstableCount)) +
  geom_bar(stat = "identity") +
  labs(x = "Parental Allele Length (bp)",
       y = "Number of Unstable Transmissions") +
  theme_classic()
ggsave_wrap("absolute_numbers_bp.png", 10, 5)
print(filtered_data, n=50)

# RU
print("RU")
filtered_data <- tr_final %>%
  mutate(LengthBin = cut(ParentAllele,
                         breaks = seq(1, ceiling(max(ParentAllele) / 5) * 5 + 5, by = 5),
                         include.lowest = TRUE,
                         right = TRUE,
                         labels = paste0("[", seq(1, ceiling(max(ParentAllele) / 5) * 5, by = 5), ", ", seq(5, ceiling(max(ParentAllele) / 5) * 5 + 4, by = 5), "]"))) %>%
  group_by(LengthBin) %>%
  summarise(TotalCount = n(),
            UnstableCount = sum(Magnitude != 0),
            Percentage = UnstableCount / TotalCount,
            .groups = 'drop')
ggplot(filtered_data, aes(x = LengthBin, y = UnstableCount)) +
  geom_bar(stat = "identity") +
  labs(x = "Parental Repeat Length (RU)",
       y = "Number of Unstable Transmissions") +
  theme_classic()
ggsave_wrap("absolute_numbers_RU.png", 10, 5)
print(filtered_data, n=50)
rm(filtered_data)

# ---------------------- OVERALL (SUMMARY) ----------------------------------- #
# Magnitude of jumps (all)
ggplot(tr_final[tr_final$Magnitude != 0,], aes(x=Magnitude)) +
  geom_histogram(binwidth = 1, aes(y = after_stat(count / sum(count))), col="white") +
  theme_classic() +
  ylab("Percentage of Unstable Transmissions (%)") +
  xlab("Magnitude of Jump (RU)") +
  scale_y_continuous(labels = percent_format()) +
  geom_vline(xintercept = 0, colour = "red", linewidth = 1, linetype = 2) +
  xlim(c(-20, 20))
print("number over 20")
nrow(tr_final[tr_final$Magnitude > 20,])
print("number under 20")
nrow(tr_final[tr_final$Magnitude < -20,])
ggsave_wrap("jump_histogram.png", 7, 5)

# Magnitude of jumps (faceted by motif length)
ggplot(tr_final[tr_final$Magnitude != 0,], aes(x = Magnitude)) +
  geom_histogram(binwidth = 1, aes(y = after_stat(count / tapply(..count.., ..PANEL.., sum)[..PANEL..])), col = "white") +
  theme_classic() +
  ylab("Percentage of Unstable Transmissions (%)") +
  scale_y_continuous(labels = percent_format()) +
  xlab("Magnitude of Jump (RU)") +
  geom_vline(xintercept = 0, colour = "red", linewidth = 0.5, linetype = 2) +
  xlim(c(-20, 20)) +
  facet_wrap(~MotifLength, nrow = 1) +
  theme(strip.text.x = element_text(size = 13, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 13))
ggsave_wrap("jump_histogram_faceted.png", 17, 4)

# # Magnitude of jumps (faceted by motif length) (BEFORE QC)
# ggplot(, aes(x=Magnitude)) +
#   geom_histogram(data = df[df$Magnitude != 0,] %>% filter(!(Gene %in% homopolymers) & Gene != "XYLT1"),
#                  binwidth = 1, alpha=0.5, position = "identity", aes(y=after_stat(count)), fill = "#615cf7") +
# geom_histogram(data = tr_final[tr_final$Magnitude != 0,], binwidth = 1, alpha=0.5, position = "identity", aes(y=after_stat(count)), fill = "green") +
#   theme_classic() +
#   ylab("Count") +
#   xlab("Magnitude of Jump (RU)") +
#   geom_vline(xintercept = 0, colour = "red", linewidth = 0.5, linetype = 2) +
#   facet_wrap(~MotifLength, nrow = 1) +
#   theme(strip.text.x = element_text(size = 13, face = "bold"),
#         axis.title = element_text(size = 17),
#         axis.text = element_text(size = 13))
# ggsave_wrap("jump_histogram_preQC_XYLT1_faceted.png", 20, 5)

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
print("stable n:")
print(nrow(tr_final %>% filter(Magnitude == 0)))
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


#---------------------- INSTABILITY (STABLE OR MUTATED)---------------------- #

print("(RU) (stable vs mutated)")
# Repeat length vs discordance rate
parentA_discordance_df <- tr_final %>%
  group_by(ParentAllele) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4) + 0.00001) %>%
  mutate(log_discordance_rate = log(discordance_rate)) %>%
  filter(total >= 1000)
model <- lm(log(discordance_rate) ~ ParentAllele, data = parentA_discordance_df)
summary(model)$coefficients[,4]
print(summary(model)$r.squared)
exp(coef(model)["ParentAllele"])
beta0 <- coef(model)[1]
beta1 <- coef(model)[2]
A <- exp(beta0)  # base discordance rate
B <- exp(beta1)  # growth factor
print(A)
print(B)
parentA_discordance_df$fitted_values <- A * (B ^ parentA_discordance_df$ParentAllele)
ggplot(parentA_discordance_df, aes(x = ParentAllele, y = discordance_rate)) +
  geom_point() +
  # geom_smooth(method = "loess", se = TRUE, color = "#b342f5", fill = "#e2c4f5", span = 0.7) +
  geom_line(aes(y = fitted_values), color = "#b342f5", size = 1, alpha = 0.5) +
  theme_classic() +
  ylab("Discordance Rate (%)") +
  xlab("Parental Repeat Length (RU)") +
  scale_y_continuous(labels = percent_format())
ggsave_wrap("discordance_parent_ru_scatter.png", 5.5, 5)
rm(parentA_discordance_df)

print(" (bp) (stable vs mutated)")
# Total BP length vs discordance rate
parentbp_discordance_df <- tr_final %>%
  group_by(ParentBpLength) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4)) %>%
  filter(total >= 1000)
ggplot(parentbp_discordance_df, aes(x=ParentBpLength, y = discordance_rate)) +
  geom_point() +
  theme_classic() +
  ylab("Discordance Rate (%)") +
  xlab("Parental Allele Length (bp)") +
  scale_y_continuous(labels = percent_format())
ggsave_wrap("discordance_parent_bp_scatter.png", 5.5, 5)
rm(parentbp_discordance_df)

print("scatter (RU, overlay) (stable vs mutated)")
parentA_discordance_df <- tr_final %>%
  group_by(ParentAllele, MotifLength) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4) + 0.00001) %>%
  filter(total >= 200)
parentA_discordance_df$MotifLength <- as.factor(parentA_discordance_df$MotifLength)
fit_model_and_predict <- function(df) {
  model <- lm(log(discordance_rate) ~ ParentAllele, data = df)
  print(summary(model)$coefficients[,4])
  print(summary(model)$r.squared)
  print(exp(coef(model)["ParentAllele"]))
  beta0 <- coef(model)[1]
  beta1 <- coef(model)[2]
  A <- beta0
  B <- beta1
  print(paste("MotifLength:", unique(df$MotifLength), "A =", A, "B =", B))
  df$fitted_values <- A * (B ^ df$ParentAllele)
  return(df)
}
fitted_df <- parentA_discordance_df %>%
  group_by(MotifLength) %>%
  group_modify(~ fit_model_and_predict(.x))
pal <- brewer.pal(5, "Set1")
names(pal) <- levels(fitted_df$MotifLength)
colScale <- scale_colour_manual(name = "Motif Length", values = pal)
ggplot(fitted_df, aes(x = ParentAllele, y = discordance_rate, colour = MotifLength)) +
  geom_point() +
  geom_line(aes(y = fitted_values, alpha = ifelse(MotifLength == 2, 0.1, 0.4)), size = 1) +
  theme_classic() +
  ylab("Discordance Rate (%)") +
  xlab("Parental Repeat Length (RU)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  colScale +
  guides(alpha = "none")
ggsave_wrap("discordance_overlay_all_motif_parent_ru_scatter_200.png", 6.5, 5)
rm(parentA_discordance_df, fitted_df)

print("scatter (bp, overlay) (stable vs mutated)")
parentbp_discordance_df <- tr_final %>%
  group_by(ParentBpLength, MotifLength) %>%
  summarise(total = n(), mutated = sum(Magnitude != 0)) %>%
  mutate(discordance_rate = round(mutated / total, 4) + 0.00001) %>%
  filter(total >= 200)
parentbp_discordance_df$MotifLength <- as.factor(parentbp_discordance_df$MotifLength)
fit_model_and_predict <- function(df) {
  model <- lm(log(discordance_rate) ~ ParentBpLength, data = df)
  print(summary(model)$coefficients[,4])
  print(summary(model)$r.squared)
  print(exp(coef(model)["ParentBpLength"]))
  beta0 <- coef(model)[1]
  beta1 <- coef(model)[2]
  A <- beta0
  B <- beta1
  print(paste("MotifLength:", unique(df$MotifLength), "A =", A, "B =", B))
  df$fitted_values <- A * (B ^ df$ParentBpLength)
  return(df)
}
fitted_df <- parentbp_discordance_df %>%
  group_by(MotifLength) %>%
  group_modify(~ fit_model_and_predict(.x))
pal <- brewer.pal(5, "Set1")
names(pal) <- levels(fitted_df$MotifLength)
colScale <- scale_colour_manual(name = "Motif Length", values = pal)
ggplot(fitted_df, aes(x = ParentBpLength, y = discordance_rate, colour = MotifLength)) +
  geom_point() +
  geom_line(aes(y = fitted_values, alpha = ifelse(MotifLength == 2, 0.1, 0.4)), size = 1) +
  theme_classic() +
  ylab("Discordance Rate (%)") +
  xlab("Parental Allele Length (bp)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  colScale +
  guides(alpha = "none")
ggsave_wrap("discordance_overlay_all_motif_parent_bp_scatter_200.png", 6.5, 5)
rm(parentbp_discordance_df, fitted_df)


# Number of expansions and contractions
#totals <- tr_final %>%
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
# ggsave_wrap("exp_contr_stacked.png", 6, 8)
# rm(totals)

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
rm(patmatexp)

# ---------------------- EXPANSION OR CONTRACTION ---------------------------- #

# Expansion percentage with parental repeat length
parentRU_perc_df <- tr_final %>%
  filter(Magnitude != 0, ParentAllele > 1) %>%
  group_by(ParentAllele) %>%
  summarise(
    total = sum(Magnitude != 0),
    total_expansion = sum(Magnitude > 0),
    perc_expansion = total_expansion / total,
    se = sqrt((perc_expansion * (1 - perc_expansion)) / total),
    ci_lower = perc_expansion - 1.96 * se,
    ci_upper = perc_expansion + 1.96 * se
  ) %>%
  filter(total >= 50)
model <- lm(perc_expansion ~ ParentAllele, data = parentRU_perc_df)
summary(model)

parentbp_perc_df <- tr_final %>%
  filter(Magnitude != 0, ParentBpLength > 6) %>%
  group_by(ParentBpLength) %>%
  summarise(
    total = sum(Magnitude != 0),
    total_expansion = sum(Magnitude > 0),
    perc_expansion = total_expansion / total,
  ) %>%
  filter(total >= 50)
model <- lm(perc_expansion ~ ParentBpLength, data = parentbp_perc_df)
summary(model)

#Stacked barplot of expansion/contraction by parental STR length
#BP
parentbp_perc_df <- tr_final %>%
  filter(Magnitude != 0, ParentBpLength > 6) %>%
  mutate(LengthBin = cut(ParentBpLength,
                         breaks = seq(7, ceiling(max(ParentBpLength) / 2) * 2, by = 2),
                         include.lowest = TRUE,
                         right = TRUE)) %>%
  group_by(LengthBin) %>%
  summarise(
    total = sum(Magnitude != 0),
    total_expansion = sum(Magnitude > 0),
    total_contraction = sum(Magnitude < 0),
    perc_expansion = total_expansion / total,
    perc_contraction = total_contraction / total
  ) %>%
  filter(total >= 50) %>%
  pivot_longer(cols = c(perc_expansion, perc_contraction),
               names_to = "MutationType",
               values_to = "Percentage") %>%
  mutate(MutationType = factor(MutationType, levels = c("perc_contraction", "perc_expansion")))
ggplot(parentbp_perc_df, aes(x = LengthBin, y = Percentage, fill = MutationType)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("perc_contraction" = "#b3cce2", "perc_expansion" = "#fab5ae"),
                    labels = c("Contraction", "Expansion")) +
  ylab("Percentage of Unstable Transmissions (%)") +
  xlab("Parental Allele Length (bp)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(fill = "Direction") +
  scale_x_discrete(breaks = c("(25,27]", "(49,51]", "(75,77]", "(99,101]", "(125,127]", "(149,151]"),
                   labels = c("25", "50", "75", "100", "125", "150"))
ggsave_wrap("expansion_parentbp_stacked.png", 8, 5)
parentA_perc_df <- tr_final %>%
  filter(Magnitude != 0, ParentAllele > 1) %>%
  group_by(ParentAllele) %>%
  summarise(
    total = sum(Magnitude != 0),
    total_expansion = sum(Magnitude > 0),
    total_contraction = sum(Magnitude < 0),
    perc_expansion = total_expansion / total,
    perc_contraction = total_contraction / total
  ) %>%
  filter(total >= 50) %>%
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
ggsave_wrap("expansion_parentA_stacked.png", 8, 5)

#stacked expansion contraction in RU faceted
parentA_perc_motif_df <- tr_final %>%
  filter(Magnitude != 0, ParentAllele > 1) %>%
  group_by(ParentAllele, MotifLength) %>%
  summarise(
    total = sum(Magnitude != 0),
    total_expansion = sum(Magnitude > 0),
    total_contraction = sum(Magnitude < 0),
    perc_expansion = total_expansion / total,
    perc_contraction = total_contraction / total
  ) %>%
  pivot_longer(cols = c(perc_expansion, perc_contraction),
               names_to = "MutationType",
               values_to = "Percentage") %>%
  mutate(MutationType = factor(MutationType, levels = c("perc_contraction", "perc_expansion")))
ggplot(parentA_perc_motif_df, aes(x = ParentAllele, y = Percentage, fill = MutationType)) +
  geom_bar(stat = "identity", colour = "black", linewidth = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("perc_contraction" = "#b3cce2", "perc_expansion" = "#fab5ae"),
                    labels = c("Contraction", "Expansion")) +
  ylab("Percentage of Unstable Transmissions (%)") +
  xlab("Parental Repeat Length (RU)") +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(fill = "Direction") +
  facet_wrap(~ MotifLength, ncol=1)
# ggsave_wrap("expansion_parentA_bymotif_stacked.png", 6, 14)
ggsave_wrap("expansion_parentA_bymotif_stacked.png", 6, 7)
ggsave_wrap("expansion_parentA_bymotif_stacked75.png", 6, 7.5)

# stacked expansion contraction in bp faceted
# # Create bins for ParentBpLength based on MotifLength
# parentbp_perc_motif_df <- tr_final %>%
#   filter(Magnitude != 0) %>%
#   mutate(LengthBin = case_when(
#     MotifLength == 2 ~ cut(ParentBpLength, breaks = seq(7, ceiling(max(ParentBpLength) / 2) * 2, by = 2), include.lowest = TRUE, right = TRUE),
#     MotifLength == 3 ~ cut(ParentBpLength, breaks = seq(7, ceiling(max(ParentBpLength) / 3) * 3, by = 3), include.lowest = TRUE, right = TRUE),
#     MotifLength == 4 ~ cut(ParentBpLength, breaks = seq(7, ceiling(max(ParentBpLength) / 4) * 4, by = 4), include.lowest = TRUE, right = TRUE),
#     TRUE ~ as.factor(ParentBpLength)
#   )) %>%
#   group_by(LengthBin, MotifLength) %>%
#   summarise(
#     total = sum(Magnitude != 0),
#     total_expansion = sum(Magnitude > 0),
#     total_contraction = sum(Magnitude < 0),
#     perc_expansion = total_expansion / total,
#     perc_contraction = total_contraction / total
#   ) %>%
#   filter(total >= 1) %>%
#   pivot_longer(cols = c(perc_expansion, perc_contraction),
#                names_to = "MutationType",
#                values_to = "Percentage") %>%
#   mutate(MutationType = factor(MutationType, levels = c("perc_contraction", "perc_expansion")))
# 
# # Plot with stacked bars and custom labels
# ggplot(parentbp_perc_motif_df, aes(x = LengthBin, y = Percentage, fill = MutationType)) +
#   geom_bar(stat = "identity", colour = "black", linewidth = 0.1) +
#   theme_classic() +
#   scale_fill_manual(values = c("perc_contraction" = "#b3cce2", "perc_expansion" = "#fab5ae"),
#                     labels = c("Contraction", "Expansion")) +
#   ylab("Percentage of Unstable Transmissions (%)") +
#   xlab("Parental Allele Length (bp)") +
#   scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
#   labs(fill = "Direction") +
#   # scale_x_discrete(breaks = c("(25,27]", "(49,51]", "(75,77]", "(99,101]", "(125,127]", "(149,151]"),
#   #                  labels = c("25", "50", "75", "100", "125", "150")) +
#   facet_wrap(~ MotifLength, ncol = 1, scales = "free_x")
# ggsave_wrap("expansion_parentbp_bymotif_stacked.png", 6, 14)
# rm(parentbp_perc_df, parentA_perc_df)
# gc(verbose = F)
# --------------------- DEGREE OF MUTATION ----------------------------------- #

# Unstable transmissions only with absolute value magnitude
degree_df <- tr_final %>%
  filter(Magnitude != 0) %>%
  mutate(AbsMagnitude = abs(Magnitude),
         AbsBpMagnitude = abs(Magnitude) * MotifLength,
         MutationType = ifelse(Magnitude > 0, "Expansion", "Contraction"))

# Table by motif length
degree_df %>%
  group_by(MotifLength) %>%
  summarise(median_magnitude = median(AbsMagnitude),
            mean_magnitude = mean(AbsMagnitude),
            sd_magnitude = sd(AbsMagnitude))

# #----------------Inaccurate analysis (inc >150bp transmissions) ----------------
# print(">150 inc but without xylt1")
#
# # Magnitude of jumps (all) (>150bp included)
# ggplot(df[df$Magnitude != 0,] %>% filter(!(Gene %in% homopolymers) & Gene != "XYLT1"), aes(x=Magnitude)) +
#   geom_histogram(binwidth = 1, colour = "white") +
#   theme_classic() +
#   ylab("Percentage of Unstable Transmissions (%)") +
#   xlab("Magnitude of Jump (RU)") +
#   geom_vline(xintercept = 0, colour = "red", linewidth = 0.5, linetype = 2) +
#   xlim(c(-50, 50 ))
# ggsave_wrap("magnitude_over150_noxylt1.png", 10, 5)
#
# # Magnitude numbers (all)
# df %>%
#   filter(Magnitude != 0) %>%
#   group_by(Magnitude) %>%
#   summarise(MagnitudeN = n()) %>%
#   mutate(bin = cut(
#     abs(Magnitude),
#     breaks = c(-Inf, 1, 2, Inf),
#     labels = c("±1", "±2", ">±2"),
#     right = TRUE )) %>%
#   aggregate(MagnitudeN ~ bin, sum) %>%
#   mutate(percentage = round(MagnitudeN / nrow(df[df$Magnitude != 0,]), 2) * 100)
#
# gc(verbose = F)

# ---------------------- Statistical tests ------------------------------------#

# STABLE VS NOT STABLE
print("STABLE VS NOT STABLE")
logreg_instability <- tr_final %>%
  mutate(MagnitudeType = ifelse(Magnitude == 0, "Stable", "Mutated")) %>%
  mutate(MagnitudeType = factor(MagnitudeType))
logreg_instability$MagnitudeType = relevel(logreg_instability$MagnitudeType, ref = "Stable")
logreg_instability$Population <- relevel(as.factor(logreg_instability$Population), ref = "CEU")

## Logistic regression as a function of motif length/ parent allele
# RU
print("RU")
logistic_model <- glm(MagnitudeType ~ ParentAllele + MotifLength,
                      family = binomial, data = logreg_instability)
summary(logistic_model)
exp(coef(logistic_model))["ParentAllele"]
exp(coef(logistic_model))["MotifLength"]
# BP
print("BP")
logistic_model <- glm(MagnitudeType ~ ParentBpLength + MotifLength,
                      family = binomial, data = logreg_instability)
summary(logistic_model)
exp(coef(logistic_model))["ParentBpLength"]
exp(coef(logistic_model))["MotifLength"]
# RU (>2bp motif)
print("RU (>2bp motif)")
logistic_model <- glm(MagnitudeType ~ ParentAllele + MotifLength,
                      family = binomial, data = logreg_instability %>% filter(MotifLength > 2))
summary(logistic_model)
exp(coef(logistic_model))["ParentAllele"]
exp(coef(logistic_model))["MotifLength"]

# Logistic regression for population
print("pop")
logistic_model <- glm(MagnitudeType ~ Population + ParentAllele,
                      family = binomial, data = logreg_instability)
summary(logistic_model)
coef_summary <- summary(logistic_model)$coefficients
tidy_results <- data.frame(
  term = rownames(coef_summary),
  estimate = coef_summary[, "Estimate"],
  std.error = coef_summary[, "Std. Error"]
) %>%
  mutate(
    conf.low = estimate - qnorm(0.975) * std.error,
    conf.high = estimate + qnorm(0.975) * std.error,
    OddsRatio = exp(estimate),
    LowerCI = exp(conf.low),
    UpperCI = exp(conf.high)
  )
print(tidy_results)
rm(logistic_model, tidy_results)

# parental sex and instability
print("parent sex instability chisq")
tab <- table(logreg_instability$ParentSex, logreg_instability$MagnitudeType)
print(tab)
chisq.test(tab)
cramers_v(tab)
rm(logreg_instability, tab)
gc(verbose = F)

## EXPANSION VS CONTRACTION
print("EXPANSION VS CONTRACTION")
logreg_direction <- tr_final %>%
  filter(Magnitude != 0,) %>%
  mutate(MagnitudeType = ifelse(Magnitude > 0, "Expansion", "Contraction")) %>%
  mutate(MagnitudeType = factor(MagnitudeType))
logreg_direction$MagnitudeType = relevel(logreg_direction$MagnitudeType, ref = "Contraction")
logreg_direction$Population <- relevel(as.factor(logreg_direction$Population), ref = "CEU")

# Logistic regression motiflength/parentallele
print("parental allele / motif length")
logistic_model <- glm(MagnitudeType ~ MotifLength + ParentAllele,
                      family = binomial, data = logreg_direction)
summary(logistic_model)
exp(coef(logistic_model)["MotifLength"])
exp(coef(logistic_model)["ParentAllele"])
rm(logistic_model)

print("parental bp length / motif length")
logistic_model <- glm(MagnitudeType ~ MotifLength + ParentBpLength,
                      family = binomial, data = logreg_direction)
summary(logistic_model)
exp(coef(logistic_model)["MotifLength"])
exp(coef(logistic_model)["ParentBpLength"])
rm(logistic_model)

## Chisqr for parental sex and direction
print("parent sex")
tab <- table(logreg_direction$ParentSex, logreg_direction$MagnitudeType)
print(tab)
chisq.test(tab)
cramers_v(tab)

# Chisqr for population and direction
print("population")
logistic_model <- glm(MagnitudeType ~ Population + ParentAllele,
                      family = binomial, data = logreg_direction)
summary(logistic_model)
coef_summary <- summary(logistic_model)$coefficients
tidy_results <- data.frame(
  term = rownames(coef_summary),
  estimate = coef_summary[, "Estimate"],
  std.error = coef_summary[, "Std. Error"]
) %>%
  mutate(
    conf.low = estimate - qnorm(0.975) * std.error,
    conf.high = estimate + qnorm(0.975) * std.error,
    OddsRatio = exp(estimate),
    LowerCI = exp(conf.low),
    UpperCI = exp(conf.high)
  )
print(tidy_results)
rm(logistic_model, tidy_results, logreg_direction, tab)
gc(verbose = F)

# DEGREE OF MUTATION
#Unstable transmissions only with absolute value magnitude
print("DEGREE OF MUTATION")
degree_df <- tr_final %>%
  filter(Magnitude != 0) %>%
  mutate(AbsMagnitude = abs(Magnitude),
         AbsBpMagnitude = abs(Magnitude) * MotifLength,
         MutationType = ifelse(Magnitude > 0, "Expansion", "Contraction"))
parental_allele_groups <- unique(degree_df$ParentAllele)

# Parental STR length
model <- lm(AbsMagnitude ~ ParentAllele + MotifLength, data = degree_df)
summary(model)
model <- lm(AbsMagnitude ~ ParentBpLength + MotifLength, data = degree_df) # Motif length
summary(model)


# Mean magnitude by parental repeat length
mean_data <- degree_df %>%
  group_by(ParentAllele, MotifLength) %>%
  summarise(MeanMagnitude = mean(AbsMagnitude, na.rm = TRUE),
            total = n()) %>%
  filter(total >= 10)
ggplot(mean_data, aes(x = ParentAllele, y = MeanMagnitude, color = factor(MotifLength))) +
  geom_point(size = 2) +
  geom_line(aes(group = factor(MotifLength)), size = 0.8, alpha = 0.5) +  
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +
  labs(x = "Parental Repeat Length (RU)",
       y = "Mean Magnitude (RU)",
       color = "Motif Length") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
ggsave_wrap("absolute_magnitude_by_length_all_10.png", 7, 5)
mean_data <- degree_df %>%
  group_by(ParentAllele, MotifLength) %>%
  summarise(MeanMagnitude = mean(Magnitude, na.rm = TRUE),
            total = n()) %>%
  filter(total >= 10)
ggplot(mean_data, aes(x = ParentAllele, y = MeanMagnitude, color = factor(MotifLength))) +
  geom_point(size = 2) +
  geom_line(aes(group = factor(MotifLength)), size = 0.8, alpha = 0.5) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) + 
  labs(x = "Parental Repeat Length (RU)",
       y = "Mean Change of Length (RU)",
       color = "Motif Length") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  ) +
  geom_hline(yintercept=0)
ggsave_wrap("magnitude_by_length_all_10.png", 7, 5)

# parental biases for magnitude
ggplot(degree_df, aes(x = factor(ParentSex), y = AbsMagnitude)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("0" = "Paternal", "1" = "Maternal")) +
  labs(x = "Parent of Origin", y = "Magnitude (RU)") +
  ylim(c(0, 25)) +
  theme_classic()
ggsave_wrap("sex_boxplot_limited_RU.png", 5, 5)
ggplot(degree_df, aes(x = factor(ParentSex), y = AbsBpMagnitude)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("0" = "Paternal", "1" = "Maternal")) +
  labs(x = "Parent of Origin", y = "Magnitude (bp)") +
  ylim(c(0, 50)) +
  theme_classic()
ggsave_wrap("sex_boxplot_limited_bp.png", 5, 5)


# Parent sex stratified kruskal's test
print("parent sex stratified kruskals for magnitude")
p_values <- c()
eta_squared_values <- c()
results <- lapply(parental_allele_groups, function(allele) {
  subset_data <- degree_df[degree_df$ParentAllele == allele, ]
  if (length(unique(subset_data$ParentSex)) > 1) {
    kruskal_result <- kruskal.test(AbsMagnitude ~ ParentSex, data = subset_data)

    H <- kruskal_result$statistic
    k <- length(unique(subset_data$ParentSex)) 
    n <- nrow(subset_data)

    eta_squared <- (H - (k - 1)) / (n - k)

    p_values <<- c(p_values, kruskal_result$p.value)
    eta_squared_values <<- c(eta_squared_values, eta_squared)
    return(list(p_value = kruskal_result$p.value, eta_squared = eta_squared))
  }
})
summary_p_value <- sumlog(na.omit(p_values))$p
summary_eta_squared <- mean(na.omit(eta_squared_values[!is.infinite(eta_squared_values)]), na.rm = TRUE)
cat("Summary p-value parent sex (Fisher's method):", summary_p_value, "\n")
cat("Summary eta-squared parent sex (mean):", summary_eta_squared, "\n")

#Ancestry stratified kruskal's test
print("population")
p_values <- c()
eta_squared_values <- c()
results <- lapply(parental_allele_groups, function(allele) {
  subset_data <- degree_df[degree_df$ParentAllele == allele, ]
  if (length(unique(subset_data$Population)) > 1) {
    kruskal_result <- kruskal.test(AbsMagnitude ~ Population, data = subset_data)

    H <- kruskal_result$statistic 
    k <- length(unique(subset_data$Population)) 
    n <- nrow(subset_data)

    eta_squared <- (H - (k - 1)) / (n - k)

    p_values <<- c(p_values, kruskal_result$p.value)
    eta_squared_values <<- c(eta_squared_values, eta_squared)
    return(list(p_value = kruskal_result$p.value, eta_squared = eta_squared))
  }
})
summary_p_value <- sumlog(na.omit(p_values))$p
summary_eta_squared <- mean(na.omit(eta_squared_values[!is.infinite(eta_squared_values)]), na.rm = TRUE)
cat("Summary p-value population (Fisher's method):", summary_p_value, "\n")
cat("Summary eta-squared population (mean):", summary_eta_squared, "\n")

## LOCI SPECIFIC stratified
print("LOCI SPECIFIC")
p_values <- c()
eta_squared_values <- c()
results <- lapply(parental_allele_groups, function(allele) {
  subset_data <- degree_df[degree_df$ParentAllele == allele, ]
  if (length(unique(subset_data$Gene)) > 1) {
    kruskal_result <- kruskal.test(AbsMagnitude ~ Gene, data = subset_data)

    H <- kruskal_result$statistic # Kruskal-Wallis test statistic
    k <- length(unique(subset_data$Gene)) # Number of groups
    n <- nrow(subset_data) # Total sample size

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

#Allele distributions for africans vs other populations
Df of alleles for each individual
parent_ids <- unique(tr_final$ParentID)
child_ids <- unique(tr_final$ChildID)
filtered_df <- tr_final[tr_final$ParentID %in% parent_ids & !(tr_final$ParentID %in% child_ids),]
filtered_df_new <- filtered_df %>%
  mutate(SampleID = ParentID,
         AlleleLength = ParentAllele,
         Ancestry = ifelse(Population %in% c("ASW", "ACB", "ESN", "MSL", "GWD", "YRI"), "African Ancestry", "Other"),
         IsPeruvian = ifelse(Population == "PEL", "Peruvian", "Other"),
         IsSriLankan = ifelse(Population == "STU", "Sri Lankan", "Other")) %>%
  dplyr::select(SampleID, Gene, AlleleLength, Ancestry, IsPeruvian, IsSriLankan)
# alleles_df <- bind_rows(filtered_df_new, original_df_new)
alleles_df <- filtered_df_new
alleles_df$AlleleLength <- as.numeric(alleles_df$AlleleLength)
alleles_df$Ancestry <- relevel(as.factor(alleles_df$Ancestry), ref = "Other")
alleles_df$IsPeruvian <- relevel(as.factor(alleles_df$IsPeruvian), ref = "Other")
alleles_df$IsSriLankan <- relevel(as.factor(alleles_df$IsSriLankan), ref = "Other")

print("AlleleLength ~ Ancestry")
wilcox.test(AlleleLength ~ Ancestry, alleles_df, conf.int = T)
print("AlleleLength ~ IsPeruvian")
wilcox.test(AlleleLength ~ IsPeruvian, alleles_df, conf.int = T)
print("AlleleLength ~ IsSriLankan")
wilcox.test(AlleleLength ~ IsSriLankan, alleles_df, conf.int = T)

means <- alleles_df %>%
  group_by(Ancestry) %>%
  summarise(mean = mean(AlleleLength),
            median = median(AlleleLength))
alleles_df %>%
  ggplot(aes(x = AlleleLength, y = after_stat(count / tapply(..count.., ..PANEL.., sum)[..PANEL..]))) +
  geom_histogram(colour = "white", binwidth = 1) +
  theme_classic() +
  ylab("Percentage of Total (%)") +
  xlab("Allele Length (RU)") +
  geom_text(data = means, aes(x = mean, y = Inf),
            label = "\U2022", vjust = 1, colour = "red") +
  geom_text(data = means, aes(x = median, y = Inf),
            label = "\U2022", vjust = 1, colour = "blue") +
  scale_y_continuous(labels = percent_format()) +
  facet_wrap(~Ancestry, ncol=1)
ggsave_wrap("african_mean.png", 5, 6)


