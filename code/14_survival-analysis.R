################################################################################
############################## SURVIVAL ANALYSIS ###############################
################################################################################

# Libraries required:
library(here)
library(tidyverse)
library(survival)
library(gtsummary)
library(survminer)


# Load the survival data:
os <- read.csv2(here("data", "survival-data.csv")) %>% 
  mutate(Patient = as.factor(Patient)) %>% 
  mutate(Dead = as.numeric(Dead)) #%>% 
 
# Load data with MS intensities 
load(here("outputs", "05_filter-out-B-samples", "05_data.filtered.Rdata")) 

# Load 45 cell-specific proteins
load(here("outputs", "11_RNAseq-comparison_II", "11_SEC-UC-long-df.RData"))
fac <- SEC_long$rowname %>% 
  unique() %>% 
  c("MRC1", "IDH2", "FAS", "ITGB8")


# Filter MS intensities for 45 cell-specific proteins
data.filtered.45 <- data.filtered %>% 
  filter(Suggested.Symbol %in% fac)
data.filtered.45[is.na(data.filtered.45)] <- 0

# Create a new directory for processed data:
dir.create(here("output", "14_survival-analysis"))


############################# UC SAMPLES ######################################

# Select UC samples
ms.uc <- data.filtered.45 %>% 
  dplyr::select(c("Suggested.Symbol", starts_with("U"))) %>%   
  column_to_rownames(var = "Suggested.Symbol") %>%   
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Patient") %>% 
  mutate(Patient = as.factor(str_remove(string = Patient, pattern = "U"))) %>% 
  full_join(x = os, y = ., by = "Patient") 

# Create survival object
ms.uc$OS <- Surv(ms.uc$OS.days, ms.uc$Dead) 

# Create table for individual cut-offs
tab <- data.frame(fac = NA, cut = NA, p.value = NA, HR = NA, LCL = NA, UCL = NA, bellow.cut = NA, above.cut = NA)

# Create table for best cut-offs in individual proteins
tab_sel_uc <- data.frame(fac = NA, cut = NA, p.value = NA, HR = NA, LCL = NA, UCL = NA, bellow.cut = NA, above.cut = NA)

#  Compute survival analysis for each cut-off in each individual protein
for (j in c(1:(length(fac)))) {
  
  # Extract data for an individual protein
  f_str <- fac[j]
  f <-  eval(parse(text = paste0("ms.uc$",f_str)))
  f <- as.numeric(as.character(factor(f)))   
  
  if (length(f) > 0)
  {
    p = NA
    cut = NA
    HR = NA
    LCL = NA
    UCL = NA
    bellow.cut = NA
    above.cut = NA
    
    # Sort unique MS intensities
    f_sort <- f %>% as.data.frame() %>% distinct()
    f.sort = sort(f_sort$.)
    
    # Compute survival analysis for individual cut-offs
    i = 1
    for (i in c(1:(length(f.sort) - 1)))
    {
      groups.f = rep(NA, length(f))
      c <- (f.sort[i] + f.sort[i + 1]) / 2
      groups.f[which(f < c)] = paste0("<", as.character(c))
      groups.f[which(f >= c)] = paste0(">", as.character(c))
      
      fit.f = survfit(OS ~ groups.f, data = ms.uc)
      test.f = survdiff(OS ~ groups.f, rho = 0, data = ms.uc)
      cox.f = coxph(OS ~ groups.f, data = ms.uc)
      
      p[i] = 1 - pchisq(test.f[[5]], df = 1)
      
      cut[i] = (f.sort[i] + f.sort[i + 1]) / 2
      
      HR[i] = summary(cox.f)$conf.int[1]
      LCL[i] = summary(cox.f)$conf.int[3]
      UCL[i] = summary(cox.f)$conf.int[4]
      
       bellow.cut[i] = fit.f[[1]][1]
       above.cut[i] = fit.f[[1]][2]
    }
    
    # Save results for each protein
    tab1 <- data.frame(fac = f_str, cut = cut, p.value = p, HR = HR, LCL = LCL, UCL = UCL, bellow.cut = bellow.cut, above.cut = above.cut) 
    tab <<- rbind(tab, tab1, c(NA))
    
    # Filter cut-off dividing into groups with minimal size  = 3 patients and set minimal cut-off at 100 000
    tab2 <- tab1 %>% 
      filter(bellow.cut > 2 & bellow.cut < 9) %>% 
      filter(cut >= 100000)
  
    #  Extract best cut-off for each protein
    tab_min <- tab2[which(min(tab2$p.value) == tab2$p.value),] 
    tab_sel_uc <<- rbind(tab_sel_uc, tab_min) 
  }
}

# Omit first row with NA values 
tab_sel_uc <- tab_sel_uc[-1,] 


#  Export results
write.csv(tab_sel_uc, file = here("outputs", "14_survival-analysis", "14_45-proteins-best-cut-offs.csv"))
