library(tidyverse)
library(magrittr)

subset_junctions <- function(population_fam, splice_junctions_fraction) {
    splice_junctions_subset <- splice_junctions_fraction[,population_fam$X1]
    splice_junctions_subset_prop10 <- splice_junctions_subset >= 0.1
    splice_junctions_subset_keep <- splice_junctions_subset[rowSums(splice_junctions_subset_prop10) >= (0.25 * ncol(splice_junctions_subset)),]
    splice_junctions_subset_keep
}

estimate_hellinger_distance <- function(grouped_df, grouping_df) {
    cluster_mat <- as.matrix(grouped_df)
    centroid <- rowMeans(cluster_mat)
    hellinger_dist <- mean(colSums((sqrt(cluster_mat) - sqrt(centroid))) ^ 2)
    tibble(cluster = grouping_df$cluster, hellinger_dist = hellinger_dist)
}

population_hellinger_distance <- function(population_fam, splice_junctions_fraction) {
    splice_junctions_subset <- select(splice_junctions_fraction, chr:full_name, one_of(population_fam$X1))
    hellinger_dists <- select(splice_junctions_subset, -chr, -start, -end, -full_name) |>
        group_by(cluster) |>
        group_map(estimate_hellinger_distance) |>
    hellinger_dists
    bind_rows()
}

esn_fam <- read_delim("../data/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x.ESN.maf05.geno01.fam", col_names = FALSE, delim = " ")
gwd_fam <- read_delim("../data/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x.GWD.maf05.geno01.fam", col_names = FALSE, delim = " ")
lwk_fam <- read_delim("../data/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x.LWK.maf05.geno01.fam", col_names = FALSE, delim = " ")
mkk_fam <- read_delim("../data/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x.MKK.maf05.geno01.fam", col_names = FALSE, delim = " ")
msl_fam <- read_delim("../data/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x.MSL.maf05.geno01.fam", col_names = FALSE, delim = " ")
yri_fam <- read_delim("../data/collapsed.hwepass01.withgnomad.var.merged.processedMKK.split1kG30x.YRI.maf05.geno01.fam", col_names = FALSE, delim = " ")

dir.create(file.path('./data', 'pop_splicing'), showWarnings = F)

splice_junctions_numerator <- vroom("./data/AFGRv2_perind_numerator.tsv")
splice_junctions_denominator <- vroom("./data/AFGRv2_perind_denominator.tsv")
splice_junctions_numerator_mat <- as.data.frame(splice_junctions_numerator) |> column_to_rownames("chrom") |> as.matrix
splice_junctions_denominator_mat <- as.data.frame(splice_junctions_denominator) |> column_to_rownames("chrom") |> as.matrix

splice_junctions_fraction <- splice_junctions_numerator_mat / (splice_junctions_numerator_mat + splice_junctions_denominator_mat)
splice_junctions_fraction[is.nan(splice_junctions_fraction)] <- 0

# Fix sample names
colnames(splice_junctions_fraction) %<>% str_replace_all("GM", "NA")

splice_junctions_df <- str_split_fixed(splice_junctions_numerator$chrom, ":", 4) |> set_colnames(c("chr", "start", "end", "cluster")) %>% as_tibble
splice_junctions_df$start %<>% as.integer
splice_junctions_df$end %<>% as.integer
splice_junctions_df$full_name <- splice_junctions_numerator$chrom

cluster_counts <- group_by(splice_junctions_df, cluster) %>% tally
cluster_counts_remove <- filter(cluster_counts, n > 10)
splice_junctions_df_keep <- filter(splice_junctions_df, !is_in(cluster, cluster_counts_remove$cluster))

splice_junctions_fraction_keep <- splice_junctions_fraction[splice_junctions_df_keep$full_name,]

subset_junctions_list <- list(esn_fam, gwd_fam, lwk_fam, mkk_fam, msl_fam, yri_fam) |>
    map(subset_junctions, splice_junctions_fraction_keep)

subset_junctions_intersect <- map(subset_junctions_list, rownames) %>% reduce(intersect)
splice_junctions_subset_df <- str_split_fixed(subset_junctions_intersect, ":", 4) |> set_colnames(c("chr", "start", "end", "cluster")) %>% as_tibble
splice_junctions_subset_df$start %<>% as.integer
splice_junctions_subset_df$end %<>% as.integer
splice_junctions_subset_df$full_name <- subset_junctions_intersect
splice_junctions_subset_tally <- group_by(splice_junctions_subset_df, cluster) %>% tally
splice_junctions_min2 <- filter(splice_junctions_subset_tally, n >= 2)
splice_junctions_min2_df <- filter(splice_junctions_subset_df, is_in(cluster, splice_junctions_min2$cluster))

splice_junctions_fraction_min2 <- splice_junctions_fraction_keep[splice_junctions_min2_df$full_name,]
splice_junctions_fraction_min2_df <- cbind(splice_junctions_min2_df, as_tibble(splice_junctions_fraction_min2)) %>% as_tibble

splice_junctions_annot <- select(splice_junctions_fraction_min2_df, chr:cluster)
write_tsv(splice_junctions_annot, "./data/splicing/splice_junctions.tsv")

hellinger_dists <- select(splice_junctions_fraction_min2_df, -chr, -start, -end, -full_name) |>
    group_by(cluster) |>
    group_map(estimate_hellinger_distance) |>
    bind_rows()

hellinger_dist_0.01 <- filter(all_dists, hellinger_dist > 0.01)

splice_junctions_fraction_min_01 <- filter(splice_junctions_fraction_min2_df, is_in(cluster, hellinger_dist_0.01$cluster))

esn_splice_junctions <- select(splice_junctions_fraction_min_01, full_name, one_of(esn_fam$X1))
colnames(esn_splice_junctions)[1] <- "cluster"
write_tsv(esn_splice_junctions, "./data/pop_splicing/ESN_splice_junctions.tsv")
gwd_splice_junctions <- select(splice_junctions_fraction_min_01, full_name, one_of(gwd_fam$X1))
colnames(gwd_splice_junctions)[1] <- "cluster"
write_tsv(gwd_splice_junctions, "./data/pop_splicing/GWD_splice_junctions.tsv")
lwk_splice_junctions <- select(splice_junctions_fraction_min_01, full_name, one_of(lwk_fam$X1))
colnames(lwk_splice_junctions)[1] <- "cluster"
write_tsv(lwk_splice_junctions, "./data/pop_splicing/LWK_splice_junctions.tsv")
mkk_splice_junctions <- select(splice_junctions_fraction_min_01, full_name, one_of(mkk_fam$X1))
colnames(mkk_splice_junctions)[1] <- "cluster"
write_tsv(mkk_splice_junctions, "./data/pop_splicing/MKK_splice_junctions.tsv")
msl_splice_junctions <- select(splice_junctions_fraction_min_01, full_name, one_of(msl_fam$X1))
colnames(msl_splice_junctions)[1] <- "cluster"
write_tsv(msl_splice_junctions, "./data/pop_splicing/MSL_splice_junctions.tsv")
yri_splice_junctions <- select(splice_junctions_fraction_min_01, full_name, one_of(yri_fam$X1))
colnames(yri_splice_junctions)[1] <- "cluster"
write_tsv(yri_splice_junctions, "./data/pop_splicing/YRI_splice_junctions.tsv")

population_dists <- list(esn_fam, gwd_fam, lwk_fam, mkk_fam, msl_fam, yri_fam) |>
    map(population_hellinger_distance, splice_junctions_fraction_min2_df)
colnames(population_dists[[1]])[2] <- "hellinger_dist_esn"
colnames(population_dists[[2]])[2] <- "hellinger_dist_gwd"
colnames(population_dists[[3]])[2] <- "hellinger_dist_lwk"
colnames(population_dists[[4]])[2] <- "hellinger_dist_mkk"
colnames(population_dists[[5]])[2] <- "hellinger_dist_msl"
colnames(population_dists[[6]])[2] <- "hellinger_dist_yri"

all_pop_dists <- reduce(population_dists, left_join)
all_dists <- left_join(hellinger_dists, all_pop_dists)
write_tsv(all_dists, "./data/pop_splicing/splice_junctions_hellinger_dists.tsv")
