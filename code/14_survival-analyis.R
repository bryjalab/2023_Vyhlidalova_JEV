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
load(here("analysis", "05_data.filtered.Rdata")) 

# Load 157 proteins specific for cancer patients
load(here("output", "09_controls-addition", "09_patient-specific-data.Rdata"))

# Filter MS intensities for 157 patient-specific proteins
data.filtered.157 <- data.filtered %>% 
  filter(Suggested.Symbol %in% patient.specific$Suggested.Symbol)
data.filtered.157[is.na(data.filtered.157)] <- 0

# Create a new directory for processed data:
dir.create(here("output", "14_survival-analysis"))

############################# SEC SAMPLES ######################################

# Select SEC samples
ms.sec <- data.filtered.157  %>% 
  dplyr::select(c(starts_with("S"), "ID", -"Suggested.Symbol")) %>% 
  column_to_rownames(var = "ID") %>%  
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Patient") %>% 
  mutate(Patient = as.factor(str_remove(string = Patient, pattern = "S"))) %>% 
  full_join(x = os, y = ., by = "Patient")

# Extract protein names
fac <- ms.sec %>% 
  dplyr::select(-c(Patient:Dead)) %>% 
  names()

# Create survival object
ms.sec$OS <- Surv(ms.sec$OS.days, ms.sec$Dead) 

# Create table for individual cut-offs
tab <- data.frame(fac = NA, cut = NA, p.value = NA, HR = NA, LCL = NA, UCL = NA)

# Create table for best cut-offs in individual proteins
tab_sel_sec <- data.frame(fac = NA, cut = NA, p.value = NA, HR = NA, LCL = NA, UCL = NA)

# Compute survival analysis for each cut-off in each individual protein
for (j in c(1:(length(fac)))) {
  
  # Extract data for an individual protein
  f_str <- fac[j]
  f <-  eval(parse(text = paste0("ms.sec$",f_str))) 
  f <- as.numeric(as.character(factor(f)))  
  
  if (length(f) > 0)
  {
    p = NA
    cut = NA
    HR = NA
    LCL = NA
    UCL = NA

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
      
      fit.f = survfit(OS ~ groups.f, data = ms.sec)
      test.f = survdiff(OS ~ groups.f, rho = 0, data=ms.sec)
      cox.f = coxph(OS ~ groups.f, data=ms.sec)
      
      p[i] = 1 - pchisq(test.f[[5]], df = 1)
      
      cut[i] = (f.sort[i] + f.sort[i + 1]) / 2

      HR[i] = summary(cox.f)$conf.int[1]
      LCL[i] = summary(cox.f)$conf.int[3]
      UCL[i] = summary(cox.f)$conf.int[4]
    }
    
    # Save results for each protein
    tab1 <- data.frame(fac = f_str, cut = cut, p.value = p, HR = HR, LCL = LCL, UCL = UCL) 
    tab <<- rbind(tab, tab1, c(NA))
    
    #  Extract best cut-off for each protein
    tab_min <- tab1[which(min(tab1$p.value) == tab1$p.value),] 
    tab_sel_sec <<- rbind(tab_sel_sec, tab_min) 
  }
}

# Omit first row with NA values and compute adjusted P-value
tab_sel_sec <- tab_sel_sec[-1,] 
tab_sel_sec$adj.p.val <- p.adjust(tab_sel_sec$p.value, method = "fdr", n = length(tab_sel_sec$p.value)) 
  

############################# UC SAMPLES ######################################

# Select UC samples
ms.uc <- data.filtered.157 %>% 
  dplyr::select(c("ID", starts_with("U"))) %>% 
  column_to_rownames(var = "ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Patient") %>% 
  mutate(Patient = as.factor(str_remove(string = Patient, pattern = "U"))) %>% 
  full_join(x = os, y = ., by = "Patient") 

# Create survival object
ms.uc$OS <- Surv(ms.uc$OS.days, ms.uc$Dead) 

# Create table for individual cut-offs
tab <- data.frame(fac = NA, cut = NA, p.value = NA, HR = NA, LCL = NA, UCL = NA)

# Create table for best cut-offs in individual proteins
tab_sel_uc <- data.frame(fac = NA, cut = NA, p.value = NA, HR = NA, LCL = NA, UCL = NA)

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
      
      
    }
    
    # Save results for each protein
    tab1 <- data.frame(fac = f_str, cut = cut, p.value = p, HR = HR, LCL = LCL, UCL = UCL) 
    tab <<- rbind(tab, tab1, c(NA))
    
    #  Extract best cut-off for each protein
    tab_min <- tab1[which(min(tab1$p.value) == tab1$p.value),] 
    tab_sel_uc <<- rbind(tab_sel_uc, tab_min) 
  }
}

# Omit first row with NA values and compute adjusted P-value
tab_sel_uc <- tab_sel_uc[-1,] 
tab_sel_uc$adj.p.val <- p.adjust(tab_sel_uc$p.value, method = "fdr", n = length(tab_sel_uc$p.value)) 

######################### COMBINE RESULTS ######################################

# Join results from SEC and UC samples
sec.uc <- inner_join(x = tab_sel_sec, y = tab_sel_uc, by = "fac", suffix = c("_sec", "_uc")) %>% 
  rename(ID = fac) 

# Add gene names
sec.uc <- data.filtered %>% 
  select(c(Suggested.Symbol, ID)) %>% 
  inner_join(x = ., y = sec.uc, by = "ID")

#  Export results
write.csv(sec.uc, file = here("output", "14_survival-analysis", "14_157-proteins-best-cut-offs.csv"))
