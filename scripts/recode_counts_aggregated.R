# setwd('/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/scripts')
# wd <- here::here()

setwd( '/oak/stanford/groups/smontgom/jaljuli/afgr_eqtls/pipe_repo/scripts')

setwd('../data')



counts_aggregated <- readr::read_tsv('counts_aggregated.txt')
reco_cols <- colnames(counts_aggregated)[startsWith( colnames(counts_aggregated), 'GM')]
reco_cols <- stringr::str_replace(reco_cols, 'GM', 'NA')
reco_cols -> colnames(counts_aggregated)[startsWith( colnames(counts_aggregated), 'GM')]

readr::write_tsv(counts_aggregated, 'counts_aggregated_GM2NA.txt')
