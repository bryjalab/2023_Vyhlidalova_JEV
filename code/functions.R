################################################################################
#############################     FUNCTIONS     ################################
################################################################################

# SCRIPT: 05_filter-out-B-samples
# Filter out patient-specifically B samples
filter_B <- function (data, B, U, S) {
  require(dplyr)
  require(lazyeval)
  B <- enquo(B)
  U <- enquo(U)
  S <- enquo(S)
  
  U.B <- data %>%
    filter(!!U != 0 & !!B == 0)
  
  S.B <- data %>%
    filter(!!S != 0 & !!B == 0)
  
  tmp <- tmp %>%
    left_join(., U.B %>% select(Suggested.Symbol, !!U), by = "Suggested.Symbol") %>%
    left_join(., S.B %>% select(Suggested.Symbol, !!S), by = "Suggested.Symbol") %>%
    select(-Suggested.Symbol)
  
  data.filtered <- cbind(data.filtered, tmp)
}

# SCRIPT: 06_methods-intersection
# Function to create intersections for each patient between UC and SEC methods
make_intersections <- function (data, U, S, number) {
  require(dplyr)
  require(lazyeval)
  
  U <- enquo(U)
  S <- enquo(S)
  
  pU <- data %>%
    filter(!is.na(!!U)) %>%
    select(Suggested.Symbol)
  
  pS <- data %>%
    filter(!is.na(!!S)) %>%
    select(Suggested.Symbol)
  
  patient_tmp <- list( pU, pS)
  names_list <- c(paste("U", number, sep = ""),
                  paste("S", number, sep = ""))
  names(patient_tmp) <- names_list
  
  patients <- append(patients, list(patient_tmp))
}