
bigtable = c()
for (path in list.files(".","features.*.tbl.rds")) {
    cat(path, "\n")
    bigtable = rbind(bigtable, readRDS(path))
}

saveRDS(bigtable, "feature_mapping_empirical.tbl")
