require(GOfuncR)

heart_development <- "GO:0007507"
cardiac_contraction <- "GO:0003015"
calcium_transport <- "GO:0006816"

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

child_nodes <- get_child_nodes(calcium_transport)

write.table(child_nodes$child_go_id,
            file = paste0(dir_docs, "calcium_transport.txt"),
            sep = "\n",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
