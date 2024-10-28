# Return list(c(allele from p1, corresponding allele in c), c(allele from p2, corresponding allele in c))
# p1: Mother p2: Father
# Assume input has no NAs or 0s
# NA in return means edge case and cannot tell optimal
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
