# Sample data
analysis1 <- c(1, 2, 3, 4, 5)
analysis2 <- c(3, 4, 5, 6, 7)
analysis3 <- c(1, 2, 6, 7, 8)

category1 <- c(1, 2, 3)
category2 <- c(4, 5, 6)
category3 <- c(7, 8, 9)

# Create a function to count matches for each analysis and category, including 'Other'
count_matches <- function(analysis, categories) {
  matches <- sapply(categories, function(category) sum(analysis %in% category))
  other_count <- length(analysis) - sum(matches)
  result_table <- data.frame(Category = c(categories, "Other"), Count = c(matches, other_count))
  return(result_table)
}

# Count matches for each analysis and category, including 'Other'
result_analysis1 <- count_matches(analysis1, list(category1, category2, category3))
result_analysis2 <- count_matches(analysis2, list(category1, category2, category3))
result_analysis3 <- count_matches(analysis3, list(category1, category2, category3))

# Display the results
print("Analysis 1:")
print(result_analysis1)

print("\nAnalysis 2:")
print(result_analysis2)

print("\nAnalysis 3:")
print(result_analysis3)





