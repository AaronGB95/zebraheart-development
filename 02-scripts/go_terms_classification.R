require(GOfuncR)

heart_growth <- "GO:0060419"
heart_development <- "GO:0007507"
cardiac_contraction <- "GO:0003015"

child_nodes <- get_child_nodes(heart_growth)

write.table(child_nodes$child_go_id,
            file = paste0(dir_docs, "heart_growth.txt"),
            sep = "\n",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

child_nodes <- get_child_nodes(heart_development)

write.table(child_nodes$child_go_id,
            file = paste0(dir_docs, "heart_development.txt"),
            sep = "\n",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

child_nodes <- get_child_nodes(cardiac_contraction)

write.table(child_nodes$child_go_id,
            file = paste0(dir_docs, "cardiac_contraction.txt"),
            sep = "\n",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)