library(tidyverse)
fam_name <- '../data/bed_bim_fam/MKK original file names and ids/MKK.filtered.biallelic.maf05.autosomal.chr.id.fam'
fam_file <- read_delim(fam_name, col_names = FALSE)
fam_ids <- fam_file$X2
id_key_path <- "../data/MKK.dictionary.txt"
id_key <- read_tsv(id_key_path, col_names = FALSE) %>% drop_na()
txt1 <- paste(id_key$X1, id_key$X2, sep = "' ~ '") # creating coding in txt1 
txt1 <- paste0("'",txt1)
txt1 <- paste0(txt1,"'")
txt2 <- paste0(c('case_match(fam_ids', txt1, ' .default = fam_ids)'), collapse = ', ')
fam_ids <- eval(parse(text=txt2))
fam_ids -> fam_file$X2

fam_name_new <- '../data/bed_bim_fam/all.MKK.30x.ID.fam'
fam_file <- write_delim(fam_file, 
                        fam_name_new, col_names = FALSE)
