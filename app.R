# Seizure Event Prediction Following Brain Injury
# Originally written April 2020 by Dominic DiSanto, Arvon Clemens, Felix Proessl as a Final Project for Advanced R Computing 
# Currently updated/maintained by: Dominic DiSanto (jdominicdisanto@gmail.com) and Nabil Awan (naa86@pitt.edu)
# Last Updated: 4/13/2020
# Version: 0.7

# For more information, please visit 
# https://github.com/domdisanto/PostTBI_SeizurePrediction_App/blob/master/README.md 


#////////////////////////////////////////////////////////////////////////////////////////////////////////////#

# Loading Packages
library(shiny)
library(magrittr)
library(shinyjs)
library(DT)
library(ggplot2)
library(ggthemes)
library(plotly)
library(shinythemes)


#///////////////////////////////////////////////////////////////////////////#

# Behind the Scenes Functions & Data Import -------------------------------

#///////////////////////////////////////////////////////////////////////////#


# _1 indicated model 1 (primary model)
# _2 indicated model 2 (RCT model)
# _3 indicated model 3 (Year-2 seizure model)
# _4 indicated model 4 (Year-2 new seizure incidence model)


# Importing Data
phats_1 <- read.csv("Model_1_Baseline_Phats.csv")
phats_graph_1 <- phats_1 %>% mutate(`Seizure Event Status` = ifelse(Observed_Event==1, "Seizure Event", "No Seizure Event"))

phats_2 <- read.csv("Model_2_RCT_Phats.csv")
phats_graph_2 <- phats_2 %>% mutate(`Seizure Event Status` = ifelse(Observed_Event==1, "Seizure Event", "No Seizure Event"))

phats_3 <- read.csv("Model_3_Yr2_Phats.csv")
phats_graph_3 <- phats_3 %>% mutate(`Seizure Event Status` = ifelse(Observed_Event==1, "Seizure Event", "No Seizure Event"))

phats_4 <- read.csv("Model_3_Yr2_NewSz_Phats.csv")
phats_graph_4 <- phats_4 %>% mutate(`Seizure Event Status` = ifelse(Observed_Event==1, "Seizure Event", "No Seizure Event"))



phats_event_1 <- phats_1[phats_1$Observed_Event==1,]
phats_nonevent_1 <- phats_1[phats_1$Observed_Event==0,]

phats_event_2 <- phats_2[phats_2$Observed_Event==1,]
phats_nonevent_2 <- phats_2[phats_2$Observed_Event==0,]

phats_event_3 <- phats_3[phats_3$Observed_Event==1,]
phats_nonevent_3 <- phats_3[phats_3$Observed_Event==0,]

phats_event_4 <- phats_4[phats_4$Observed_Event==1,]
phats_nonevent_4 <- phats_4[phats_4$Observed_Event==0,]




# Calculating Phat and Calculating Status
pred_prob_func_1 <- function(input) {
  pred_prob_1 <- 1 / (1 + exp(-(input$craniectomy_1*1.192+
                                  input$acute_sz_1*1.172+
                                  (input$contusion_num_1=="1")*0.212 + (input$contusion_num_1=="2")*0.452 + (input$contusion_num_1=="3")*0.556 + (input$contusion_num_1=="4")*0.574 +
                                  input$craniotomy_1*0.496+
                                  input$psychhosp_1*0.482+
                                  input$incarcerate_1*0.427+
                                  input$ctfrag_1*0.364+
                                  input$sdh_1*0.359+
                                  input$drug_1*0.263+
                                  input$hosp_alcohol_1*0.188+
                                  input$injury_severity_1*0.171+
                                  input$sah_1*0.114+
                                  as.numeric(input$acute_los_1)*0.003+
                                  -2.773
                                )))
  classification_1 <- ifelse(pred_prob_1 >= input$thresholdslider_1, "PTS Event", "No Event")
  return(list(pred_prob_1=pred_prob_1, classification_1=classification_1))
}


pred_prob_func_2 <- function(input) {
  pred_prob_2 <- 1 / (1 + exp(-(input$craniectomy_2*1.232+
                                  (input$contusion_num_2=="1")*0.23 + (input$contusion_num_2=="2")*0.491 + (input$contusion_num_2=="3")*0.638 + (input$contusion_num_2=="4")*0.668 +
                                  input$craniotomy_2*0.544+
                                  input$psychhosp_2*0.53+
                                  (input$race_2=="African-American")*0.449 + (input$race_2=="Other")*0.000 +
                                  input$sdh_2*0.362+
                                  input$incarcerate_2*0.327+
                                  input$ctfrag_2*0.321+
                                  input$drug_2*0.231+
                                  input$hosp_alcohol_2*0.209+
                                  input$injury_severity_2*0.206+
                                  -2.607
  )))
  classification_2 <- ifelse(pred_prob_2 >= input$thresholdslider_2, "PTS Event", "No Event")
  return(list(pred_prob_2=pred_prob_2, classification_2=classification_2))
}


pred_prob_func_3 <- function(input) {
  pred_prob_3 <- 1 / (1 + exp(-(input$yr1_sz_3*3.067+
                                  input$craniectomy_3*0.981+
                                  input$acute_sz_3*0.684+
                                  input$drug_3*0.629+
                                  input$ctfrag_3*0.434+
                                  input$craniotomy_3*0.414+
                                  as.numeric(input$acute_los_3)*0.007+
                                  -3.371
  )))  
  classification_3 <- ifelse(pred_prob_3 >= input$thresholdslider_3, "PTS Event", "No Event")
  return(list(pred_prob_3=pred_prob_3, classification_3=classification_3))
}


pred_prob_func_4 <- function(input) {
  pred_prob_4 <- 1 / (1 + exp(-(input$craniectomy_4*1.047+
                                  input$acute_sz_4*0.892+
                                  (input$contusion_num_4=="1")*0 + (input$contusion_num_4=="2")*0.603 + (input$contusion_num_4=="3")*0.846 + (input$contusion_num_4=="4")*0.866 +
                                  input$drug_4*0.653+
                                  input$craniotomy_4*0.563+
                                  as.numeric(input$acute_los_4)*0.005+
                                  -3.801
  )))  
  classification_4 <- ifelse(pred_prob_4 >= input$thresholdslider_4, "PTS Event", "No Event")
  return(list(pred_prob_4=pred_prob_4, classification_4=classification_4))
}


# Determining theoretical patient quantile
quantile_fun_1 <- function(input=NULL, value=NULL) {
  if((is.null(input) & is.null(value)) | (!is.null(input) & !is.null(value))) stop("Please enter information for either only input OR value. Not both")
  if(is.null(input)) pred_prob_val_1 <- value
  if(is.null(value)) pred_prob_val_1 <-  pred_prob_func_1(input)$pred_prob_1
  
  total_percent_1 <- nrow(phats_1[phats_1[,"Phat"] < pred_prob_val_1,]) / nrow(phats_1)  
  pte_percent_1 <- nrow(phats_event_1[phats_event_1[,"Phat"] < pred_prob_val_1,]) / nrow(phats_event_1)  
  nonpte_percent_1 <- nrow(phats_nonevent_1[phats_nonevent_1[,"Phat"] < pred_prob_val_1,]) / nrow(phats_nonevent_1)  
  
  return(list(total_percent_1=total_percent_1, pte_percent_1=pte_percent_1, nonpte_percent_1=nonpte_percent_1))
}

quantile_fun_2 <- function(input=NULL, value=NULL) {
  if((is.null(input) & is.null(value)) | (!is.null(input) & !is.null(value))) stop("Please enter information for either only input OR value. Not both")
  if(is.null(input)) pred_prob_val_2 <- value
  if(is.null(value)) pred_prob_val_2 <-  pred_prob_func_2(input)$pred_prob_2
  
  total_percent_2 <- nrow(phats_2[phats_2[,"Phat"] < pred_prob_val_2,]) / nrow(phats_2)  
  pte_percent_2 <- nrow(phats_event_2[phats_event_2[,"Phat"] < pred_prob_val_2,]) / nrow(phats_event_2)  
  nonpte_percent_2 <- nrow(phats_nonevent_2[phats_nonevent_2[,"Phat"] < pred_prob_val_2,]) / nrow(phats_nonevent_2)  
  
  return(list(total_percent_2=total_percent_2, pte_percent_2=pte_percent_2, nonpte_percent_2=nonpte_percent_2))
}

quantile_fun_3 <- function(input=NULL, value=NULL) {
  if((is.null(input) & is.null(value)) | (!is.null(input) & !is.null(value))) stop("Please enter information for either only input OR value. Not both")
  if(is.null(input)) pred_prob_val_3 <- value
  if(is.null(value)) pred_prob_val_3 <-  pred_prob_func_3(input)$pred_prob_3
  
  total_percent_3 <- nrow(phats_3[phats_3[,"Phat"] < pred_prob_val_3,]) / nrow(phats_3)  
  pte_percent_3 <- nrow(phats_event_3[phats_event_3[,"Phat"] < pred_prob_val_3,]) / nrow(phats_event_3)  
  nonpte_percent_3 <- nrow(phats_nonevent_3[phats_nonevent_3[,"Phat"] < pred_prob_val_3,]) / nrow(phats_nonevent_3)  
  
  return(list(total_percent_3=total_percent_3, pte_percent_3=pte_percent_3, nonpte_percent_3=nonpte_percent_3))
}

quantile_fun_4 <- function(input=NULL, value=NULL) {
  if((is.null(input) & is.null(value)) | (!is.null(input) & !is.null(value))) stop("Please enter information for either only input OR value. Not both")
  if(is.null(input)) pred_prob_val_4 <- value
  if(is.null(value)) pred_prob_val_4 <-  pred_prob_func_4(input)$pred_prob_4
  
  total_percent_4 <- nrow(phats_4[phats_4[,"Phat"] < pred_prob_val_4,]) / nrow(phats_4)  
  pte_percent_4 <- nrow(phats_event_4[phats_event_4[,"Phat"] < pred_prob_val_4,]) / nrow(phats_event_4)  
  nonpte_percent_4 <- nrow(phats_nonevent_4[phats_nonevent_4[,"Phat"] < pred_prob_val_4,]) / nrow(phats_nonevent_4)  
  
  return(list(total_percent_4=total_percent_4, pte_percent_4=pte_percent_4, nonpte_percent_4=nonpte_percent_4))
}


# Determining threshold specific sensitivity & specificity manually using only the phat file
model_per_1 <- function(input) {
  new_sens_1 <- round(sum((phats_1$Phat > input$thresholdslider_1) & phats_1$Observed_Event==1) / (sum((phats_1$Phat > input$thresholdslider_1) & phats_1$Observed_Event==1) + 
                                                                                             sum((phats_1$Phat < input$thresholdslider_1) & phats_1$Observed_Event==1)), 3) # TP / (TP + FN) 
  new_spec_1 <- round(sum((phats_1$Phat < input$thresholdslider_1) & phats_1$Observed_Event==0) / (sum((phats_1$Phat < input$thresholdslider_1) & phats_1$Observed_Event==0) + 
                                                                                             sum((phats_1$Phat > input$thresholdslider_1) & phats_1$Observed_Event==0)), 3) # TN / (TN + FP) 
  # sens automatic using perform_mat    # round(perform_mat[which.min(abs(perform_mat$threshold-input$thresholdslider)),3], 4)
  new_ppv_1 <- round(sum((phats_1$Phat > input$thresholdslider_1) & phats_1$Observed_Event==1) / (sum((phats_1$Phat > input$thresholdslider_1) & phats_1$Observed_Event==1) + 
                                                                                            sum((phats_1$Phat > input$thresholdslider_1) & phats_1$Observed_Event==0)), 3) # TP / (TP + FP) 
  new_npv_1 <- round(sum((phats_1$Phat < input$thresholdslider_1) & phats_1$Observed_Event==0) / (sum((phats_1$Phat < input$thresholdslider_1) & phats_1$Observed_Event==0) + 
                                                                                            sum((phats_1$Phat < input$thresholdslider_1) & phats_1$Observed_Event==1)), 3) # TN / (TN + FN) 
  return(list(sens_1=new_sens_1, spec_1=new_spec_1, ppv_1=new_ppv_1, npv_1=new_npv_1))
}

model_per_2 <- function(input) {
  new_sens_2 <- round(sum((phats_2$Phat > input$thresholdslider_2) & phats_2$Observed_Event==1) / (sum((phats_2$Phat > input$thresholdslider_2) & phats_2$Observed_Event==1) + 
                                                                                                     sum((phats_2$Phat < input$thresholdslider_2) & phats_2$Observed_Event==1)), 3) # TP / (TP + FN) 
  new_spec_2 <- round(sum((phats_2$Phat < input$thresholdslider_2) & phats_2$Observed_Event==0) / (sum((phats_2$Phat < input$thresholdslider_2) & phats_2$Observed_Event==0) + 
                                                                                                     sum((phats_2$Phat > input$thresholdslider_2) & phats_2$Observed_Event==0)), 3) # TN / (TN + FP) 
  # sens automatic using perform_mat    # round(perform_mat[which.min(abs(perform_mat$threshold-input$thresholdslider)),3], 4)
  new_ppv_2 <- round(sum((phats_2$Phat > input$thresholdslider_2) & phats_2$Observed_Event==1) / (sum((phats_2$Phat > input$thresholdslider_2) & phats_2$Observed_Event==1) + 
                                                                                                    sum((phats_2$Phat > input$thresholdslider_2) & phats_2$Observed_Event==0)), 3) # TP / (TP + FP) 
  new_npv_2 <- round(sum((phats_2$Phat < input$thresholdslider_2) & phats_2$Observed_Event==0) / (sum((phats_2$Phat < input$thresholdslider_2) & phats_2$Observed_Event==0) + 
                                                                                                    sum((phats_2$Phat < input$thresholdslider_2) & phats_2$Observed_Event==1)), 3) # TN / (TN + FN) 
  return(list(sens_2=new_sens_2, spec_2=new_spec_2, ppv_2=new_ppv_2, npv_2=new_npv_2))
}

model_per_3 <- function(input) {
  new_sens_3 <- round(sum((phats_3$Phat > input$thresholdslider_3) & phats_3$Observed_Event==1) / (sum((phats_3$Phat > input$thresholdslider_3) & phats_3$Observed_Event==1) + 
                                                                                                     sum((phats_3$Phat < input$thresholdslider_3) & phats_3$Observed_Event==1)), 3) # TP / (TP + FN) 
  new_spec_3 <- round(sum((phats_3$Phat < input$thresholdslider_3) & phats_3$Observed_Event==0) / (sum((phats_3$Phat < input$thresholdslider_3) & phats_3$Observed_Event==0) + 
                                                                                                     sum((phats_3$Phat > input$thresholdslider_3) & phats_3$Observed_Event==0)), 3) # TN / (TN + FP) 
  # sens automatic using perform_mat    # round(perform_mat[which.min(abs(perform_mat$threshold-input$thresholdslider)),3], 4)
  new_ppv_3 <- round(sum((phats_3$Phat > input$thresholdslider_3) & phats_3$Observed_Event==1) / (sum((phats_3$Phat > input$thresholdslider_3) & phats_3$Observed_Event==1) + 
                                                                                                    sum((phats_3$Phat > input$thresholdslider_3) & phats_3$Observed_Event==0)), 3) # TP / (TP + FP) 
  new_npv_3 <- round(sum((phats_3$Phat < input$thresholdslider_3) & phats_3$Observed_Event==0) / (sum((phats_3$Phat < input$thresholdslider_3) & phats_3$Observed_Event==0) + 
                                                                                                    sum((phats_3$Phat < input$thresholdslider_3) & phats_3$Observed_Event==1)), 3) # TN / (TN + FN) 
  return(list(sens_3=new_sens_3, spec_3=new_spec_3, ppv_3=new_ppv_3, npv_3=new_npv_3))
}

model_per_4 <- function(input) {
  new_sens_4 <- round(sum((phats_4$Phat > input$thresholdslider_4) & phats_4$Observed_Event==1) / (sum((phats_4$Phat > input$thresholdslider_4) & phats_4$Observed_Event==1) + 
                                                                                                     sum((phats_4$Phat < input$thresholdslider_4) & phats_4$Observed_Event==1)), 3) # TP / (TP + FN) 
  new_spec_4 <- round(sum((phats_4$Phat < input$thresholdslider_4) & phats_4$Observed_Event==0) / (sum((phats_4$Phat < input$thresholdslider_4) & phats_4$Observed_Event==0) + 
                                                                                                     sum((phats_4$Phat > input$thresholdslider_4) & phats_4$Observed_Event==0)), 3) # TN / (TN + FP) 
  # sens automatic using perform_mat    # round(perform_mat[which.min(abs(perform_mat$threshold-input$thresholdslider)),3], 4)
  new_ppv_4 <- round(sum((phats_4$Phat > input$thresholdslider_4) & phats_4$Observed_Event==1) / (sum((phats_4$Phat > input$thresholdslider_4) & phats_4$Observed_Event==1) + 
                                                                                                    sum((phats_4$Phat > input$thresholdslider_4) & phats_4$Observed_Event==0)), 3) # TP / (TP + FP) 
  new_npv_4 <- round(sum((phats_4$Phat < input$thresholdslider_4) & phats_4$Observed_Event==0) / (sum((phats_4$Phat < input$thresholdslider_4) & phats_4$Observed_Event==0) + 
                                                                                                    sum((phats_4$Phat < input$thresholdslider_4) & phats_4$Observed_Event==1)), 3) # TN / (TN + FN) 
  return(list(sens_4=new_sens_4, spec_4=new_spec_4, ppv_4=new_ppv_4, npv_4=new_npv_4))
}


#/////////////////////////////////////////////#
###### User Interface & Input Section ########
#/////////////////////////////////////////////#

ui <-  fluidPage(
  navbarPage("Seizure Prediction Following TBI", theme=shinytheme("lumen"),
             
             
             ### Instructional tab
             
             tabPanel("Instructional", fluid=T, 
                      
                      h2(strong("Predicting Epilepsy among Individuals with Moderate-to-Severe Traumatic Brain Injury")),

                      h3(strong("Background")),
                      h4("The present application uses results from Amy Wagner's research group from the University of Pittsburgh School of Medicine's Department of Physical Medicine & Rehabilitation. The prediction model used comes from the unpublished work \"A Prognostic Model of Seizure Events Following Moderate-to-Severe Traumatic Brain Injury\" (Awan et al, 2024; In Review)."),
                      
                      h4("The prediction calculation is derived from a group LASSO derived logistic regression, fit in a cohort of over 6,089  individuals with moderate-to-severe traumatic brain injury, who enrolled in the Traumatic Brain Injury Model Systems from 2010 to 2018. This page contains information related to the model development, input data, and the required data necessary to apply this tool to", 
                         em("your specific patient.")), 
                      
                      h3(strong("Current Case Uses and Considerations")),
                      h4(strong("External Validation Research Collaborations:"), " This PTE calculator is currently being made available for research use with potential collaborators for the purposes of external validation of the model and model refinement."),
                      
                      h4(strong("Default Classification Threshold Considerations:"), " The default PTE Classification Thresholds are set to maximize Negative Predictive Value using a Specificity Cut-point of ~60% for each model.  This choice in default PTE Classification is recommended when considering prognostication and/or developing/testing clinical decision algorithms wherein high false negative rates (i.e. someone develops PTE who was not identified as such with the PTE calculator) would have a consequential impact on recovery and management. However, the calculator tool allows for investigators to examine how changes in the PTE classification threshold affect NPV, PPV and other model indicators, which may be appropriate for some research questions under consideration for development of clinical decision trees or algorithms using this tool."),
                      
                      h4(strong("Missing Data Considerations:"),  " The Data Elements incorporated into this calculator are considered to be commonly available patient/research participant clinical history or dataset.  Further model testing is being pursued to systematically evaluate the impact of missing data on fidelity of this PTE calculator and test multiple methods of imputation aimed at minimizing the impact of missinginess on PTE classification for moderate-to-severe TBI survivors."),
                      
                      h3(strong("Abbreviations:")),
                      h4(strong("CT"), " = Computed Tomography", br(),
                         strong("EDH"), " = Epidural Hematoma", br(), 
                         strong("GCS"), " = Glasgow Coma Scale", br(),
                         strong("LASSO"), " = Least Absolute Shrinkage and Selection Operator", br(),
                         strong("NPV"), " = Negative Predictive Value", br(),
                         strong("PPV"), " = Positive Predictive Value", br(),
                         strong("PTE"), " = Post-Traumatic Epilepsy", br(),
                         strong("SDH"), " = Subdural Hematoma", br(),
                         strong("TBI"), " = Traumatic Brain Injury", br(),
                         strong("TBIMS"), " = Traumatic Brain Injury Model System"),
                      
                      h3(strong("Variable Descriptions & Data Sources:")),
                      h4(em("Below is a brief summary of the information required for the calculation of an individual's predicted probability and possible sources from which the data can be obtained. 
                                 The specific variables and their assessment via patient interview and/or medical record review are more thoroughly discussed on the ", 
                            a("TBIMS website.", href="https://www.tbindsc.org/Syllabus.aspx"))),
                      
                      # Variables & Data Sources 
                      datatable(data.frame(
                        Variable = c("Race", "SDH", "Intracranial Fragments", "EDH", "Contusions",
                                     "Alcohol at Hospitalization", "Acute Seizures", "Craniotomy/Craniectomy",
                                     "Pre-Injury History of Incarceration", "Pre-Injury History of Psychiatric Hospitalization"),
                        `Data Source` = c("Self-reported", # Race
                                          "Review of CT imaging and/or radiologist summary report", # SDH
                                          "Review of CT imaging and/or radiologist summary report; Alternatively present of penetrating brain injury suffices as evidence of present fragments", # Intracranial Fragments
                                          "Review of CT imaging and/or radiologist summary report", # EDH
                                          "Review of CT imaging and/or radiologist summary report; Presence of 4 or more contusions (even if total number is unknown) are included in the \"4+\" category ", # Contusions
                                          "ICD-9 and ICD-10 coding", # ETOH/Alcohol
                                          "ICD-9 and ICD-10 coding", # Acute Seizures
                                          "Acute care medical records (TBIMS Form-1)", # Craniotomy/Craniectomy # Ask Raj/Dan/Team to confirm.
                                          "Pre-injury history (participant or proxy)", # Incarceration
                                          "Pre-injury history (participant or proxy)" # Psychiatric Hospitalization
                        )),
                        
                        rownames = F, options=list(dom='t', ordering=F, headerCallback = DT::JS(
                          "function(thead) {",
                          "  $(thead).css('font-size', '25px');",
                          "}"))) %>% 
                        formatStyle(columns = 1, fontWeight = 'bold') %>% 
                        formatStyle(columns = c(1:2), `font-size`="17px")
                      
             ),  # end to first tabpanel
             
### Baseline prediction model tab
             
             tabPanel("Baseline Prediction Model", fluid=TRUE, icon=icon("chart-bar"),
                      # Title and Headers
                      h1(strong("Predicting Seizure Events in Moderate-to-Severe Traumatic Brain Injury")),
                      h2(strong("Baseline Model")),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      
                      sidebarLayout(
                        
                        
                        # Sidebar (Input Section)
                        sidebarPanel(
                          # Imaging Inputs (SDH, SAH, EDH, CT Fragment)    
                          h2(strong("Input Dashboard")),
                          shinyjs::useShinyjs(),
                          id = "side-panel_1",
                          
                          actionButton("reset_button_yr1_1", "New Patient/Reset to Defaults"),
                          

                          h3(strong("Imaging Data")),
                          checkboxInput("sdh_1", label = ("Subdural Hematoma"), value = NULL),
                          checkboxInput("sah_1", label = ("Subarachnoid Hemorrhage"), value = NULL),
                          checkboxInput("ctfrag_1", label = "Intracranial Fragments", value = NULL),

                          
                          # Contusions Dropdown
                          selectInput("contusion_num_1", h3(strong("Number of Contusions:")),
                                      c("0" = 0,
                                        "1" = 1,
                                        "2" = 2,
                                        "3" = 3,
                                        "4+" = 4)),
                          
                          
                          # Hospitalization/Injury Data (injury severity [PTA, TTFC, GCS], Alcohol, Acute Seizures, Craniotomy/Craniectomy)
                          h3(strong("Injury/Hospitalization Information")),
                          textInput("acute_los_1", label = "Acute LOS (days)", value = NULL),
                          checkboxInput("acute_sz_1", label = "Acute Seizures", value = NULL),
                          checkboxInput("craniotomy_1", label="Craniotomy", value=NULL),
                          checkboxInput("craniectomy_1", label="Craniectomy", value=NULL),
                          checkboxInput("injury_severity_1", label="Injury Severity (Check If Severe, Don't Check If Moderate)", value=NULL),
                          
                          # Pre-injury medical history information (Incarceration, Psych Hospitalization, Neurodegenerative Disease)
                          h3(strong("Pre-Injury Information")),
                          checkboxInput("incarcerate_1", label = "Pre-Injury History of Incarcaration", value = NULL),
                          checkboxInput("drug_1", label = "Pre-Injury Drug Use", value = NULL),
                          checkboxInput("hosp_alcohol_1", label = "Alcohol Use at Hospitalization", value = NULL),
                          checkboxInput("psychhosp_1", label = "Pre-Injury Psychiatric Hospitalization/Institutionalization", value = NULL)

                        ), # closure to sidebarPanel(
                        
                        
                        
                        # mainPanel: output section for tabular and plotted results
                        mainPanel(
                          h3(strong("Tabular Summary")),
                          tags$head((tags$style("#tabular_header{font-size:22px; font-weight: bold;}"))),
                          
                          h4(em(strong("Patient Level Results"))),
                          dataTableOutput("patient_table_1"),
                          h4(em(strong("Model Performance Information"))),
                          dataTableOutput("model_table_1"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          sliderInput("thresholdslider_1", "Manually Set Classification Threshold", min=0, max=1,value = 0.1242, width = "600px", round = -4, step = 0.001),
                          actionButton("reset_threshold_1", "Reset to Default Classification Threshold"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          h3(strong("Visual Summary")),
                          tags$head((tags$style("#graph_header_1{font-size:22px; font-weight: bold;}"))),#,
                          
                          plotlyOutput('int_plot_1')
                          
                        ) # closure to mainPanel(
                      ), # closure to sidebarLayout(    
                      
                      
                      
                      # Footer Equations/Notes
                      hr(),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h3("Calculations & Analytic Background"),
                      tags$head((tags$style("#calculation_header{font-size:22px; font-weight: bold;}"))),
                      
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Logistic regression formulation of the log odds-ratio (or logit):"),
                      tags$head((tags$style("#logit_formula_1_text{font-size:16px;}"))),
                      
                      withMathJax(paste("$$logit(\\hat{p}) = -2.773 + (Craniectomy) \\times 1.192 + 
                       (Acute Seizures) \\times 1.172 +
                       (Contusions (1)) \\times 0.212 + (Contusions (2)) \\times 0.452 +
                       (Contusions (3)) \\times 0.556 + (Contusions (4+)) \\times 0.574 + $$
                       \n
                       $$(Intracranial Fragments) \\times 0.364 +
                       (Psych Hospitalization) \\times 0.482 +
                       (Craniotomy) \\times 0.496 +  
                       (SDH) \\times 0.359 +
                       (SAH) \\times 0.114 + 
                       (Injury Severity) \\times 0.171 + $$
                       \n
                       $$(Incarceration) \\times 0.427 +
                       (Preinjury Drug Use) \\times 0.263 + 
                       (Alcohol at Hospitalization) \\times 0.188 +
                       (Acute LOS) \\times 0.003 $$")),
                                            
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Individual predicted probability values then calculated using the following formula:"),
                      tags$head((tags$style("#phat_text{font-size:16px;}"))),
                      
                      withMathJax(paste("$$\\hat{p} = \\frac{1}{1+exp(-logit(\\hat{p}))}$$")),
                      
             ), # end to second tabpanel
             
             
             
### RCT prediction model tab
             
             tabPanel("RCT Prediction Model", fluid=TRUE, icon=icon("chart-bar"),
                      # Title and Headers
                      h1(strong("Predicting Seizure Events in Moderate-to-Severe Traumatic Brain Injury")),
                      h2(strong("RCT Model")),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      
                      sidebarLayout(
                        
                        
                        # Sidebar (Input Section)
                        sidebarPanel(
                          # Imaging Inputs (SDH, EDH, CT Fragment)    
                          h2(strong("Input Dashboard")),
                          shinyjs::useShinyjs(),
                          id = "side-panel_2",
                          
                          actionButton("reset_button_yr1_2", "New Patient/Reset to Defaults"),
                          
                          # Race
                          selectInput("race_2", h3(strong("Race")),
                                      c("White" = 0,
                                        "African-American" = 1,
                                        "Other" = 2)),
                          
                          
                          h3(strong("Imaging Data")),
                          checkboxInput("sdh_2", label = ("Subdural Hematoma"), value = NULL),
                          checkboxInput("ctfrag_2", label = "Intracranial Fragments", value = NULL),
                          
                          
                          # Contusions Dropdown
                          selectInput("contusion_num_2", h3(strong("Number of Contusions:")),
                                      c("0" = 0,
                                        "1" = 1,
                                        "2" = 2,
                                        "3" = 3,
                                        "4+" = 4)),
                          
                          
                          # Hospitalization/Injury Data (injury severity [PTA, TTFC, GCS], Alcohol, Acute Seizures, Craniotomy/Craniectomy)
                          h3(strong("Injury/Hospitalization Information")),
                          checkboxInput("craniotomy_2", label="Craniotomy", value=NULL),
                          checkboxInput("craniectomy_2", label="Craniectomy", value=NULL),
                          checkboxInput("injury_severity_2", label="Injury Severity (Check If Severe, Don't Check If Moderate)", value=NULL),

                          
                          # Pre-injury medical history information (Incarceration, Psych Hospitalization, Neurodegenerative Disease)
                          h3(strong("Pre-Injury Information")),
                          checkboxInput("incarcerate_2", label = "Pre-Injury History of Incarcaration", value = NULL),
                          checkboxInput("drug_2", label = "Pre-Injury Drug Use", value = NULL),
                          checkboxInput("hosp_alcohol_2", label="Alcohol at Hospitalization", value=NULL),
                          checkboxInput("psychhosp_2", label = "Pre-Injury Psychiatric Hospitalization/Institutionalization", value = NULL)
                          
                        ), # closure to sidebarPanel(
                        
                        
                        
                        # mainPanel: output section for tabular and plotted results
                        mainPanel(
                          h3(strong("Tabular Summary")),
                          tags$head((tags$style("#tabular_header{font-size:22px; font-weight: bold;}"))),
                          
                          h4(em(strong("Patient Level Results"))),
                          dataTableOutput("patient_table_2"),
                          h4(em(strong("Model Performance Information"))),
                          dataTableOutput("model_table_2"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          sliderInput("thresholdslider_2", "Manually Set Classification Threshold", min=0, max=1,value = 0.1284, width = "600px", round = -4, step = 0.001),
                          actionButton("reset_threshold_2", "Reset to Default Classification Threshold"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          h3(strong("Visual Summary")),
                          tags$head((tags$style("#graph_header_2{font-size:22px; font-weight: bold;}"))),
                          
                          plotlyOutput('int_plot_2')
                          
                        ) # closure to mainPanel(
                      ), # closure to sidebarLayout(    
                      
                      
                      
                      # Footer Equations/Notes
                      hr(),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h3("Calculations & Analytic Background"),
                      tags$head((tags$style("#calculation_header{font-size:22px; font-weight: bold;}"))),
                      
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Logistic regression formulation of the log odds-ratio (or logit):"),
                      tags$head((tags$style("#logit_formula_2_text{font-size:16px;}"))),
                      
                      
                      withMathJax(paste("$$logit(\\hat{p}) = -2.607 + (Craniectomy) \\times 1.232 + 
                       (Contusions (1)) \\times 0.23 + (Contusions (2)) \\times 0.491 +
                       (Contusions (3)) \\times 0.638 + (Contusions (4+)) \\times 0.668 + $$
                       \n
                       $$(Intracranial Fragments) \\times 0.321 +
                       (Psych Hospitalization) \\times 0.53 +
                       (Race (African-American)) \\times 0.381 +
                       (Race (Other)) \\times 0.000 +
                       (Craniotomy) \\times 0.544 +  
                       (SDH) \\times 0.362 + $$
                       \n
                       $$(Incarceration) \\times 0.327 +
                       (Injury Severity) \\times 0.206 +
                       (Alcohol at Hospitalization) \\times 0.209 +
                       (Preinjury Drug Use) \\times 0.231  $$")),
                      
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Individual predicted probability values then calculated using the following formula:"),
                      tags$head((tags$style("#phat_text{font-size:16px;}"))),
                      
                      withMathJax(paste("$$\\hat{p} = \\frac{1}{1+exp(-logit(\\hat{p}))}$$")),
                      
             ), # end to third tabpanel
             
             
             
### Year-2 seizure prediction model tab
             
             tabPanel("Year-2 Seizure Prediction Model", fluid=TRUE, icon=icon("chart-bar"),
                      # Title and Headers
                      h1(strong("Predicting Seizure Events in Moderate-to-Severe Traumatic Brain Injury")),
                      h2(strong("Year-2 Seizure (New or Recurrent) Model")),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      
                      sidebarLayout(
                        
                        
                        # Sidebar (Input Section)
                        sidebarPanel(
                          # Imaging Inputs (SDH, EDH, CT Fragment)    
                          h2(strong("Input Dashboard")),
                          shinyjs::useShinyjs(),
                          id = "side-panel_3",
                          
                          actionButton("reset_button_yr1_3", "New Patient/Reset to Defaults"),

                          h3(strong("Previous Seizure")),
                          checkboxInput("yr1_sz_3", label = "Seizure in Year-1", value = NULL),
                          
                          h3(strong("Imaging Data")),
                          checkboxInput("ctfrag_3", label = "Intracranial Fragments", value = NULL),

                          
                          # Hospitalization/Injury Data (injury severity [PTA, TTFC, GCS], Alcohol, Acute Seizures, Craniotomy/Craniectomy)
                          h3(strong("Injury/Hospitalization Information")),
                          textInput("acute_los_3", label = "Acute LOS (days)", value = NULL),
                          checkboxInput("acute_sz_3", label = "Acute Seizures", value = NULL),
                          checkboxInput("craniectomy_3", label="Craniectomy", value=NULL),
                          checkboxInput("craniotomy_3", label="Craniotomy", value=NULL),

                          # Pre-injury medical history information (Incarceration, Psych Hospitalization, Neurodegenerative Disease)
                          h3(strong("Pre-Injury Information")),
                          checkboxInput("drug_3", label = "Pre-Injury Drug Use", value = NULL),
                          
                        ), # closure to sidebarPanel(
                        
                        
                        
                        # mainPanel: output section for tabular and plotted results
                        mainPanel(
                          h3(strong("Tabular Summary")),
                          tags$head((tags$style("#tabular_header{font-size:22px; font-weight: bold;}"))),
                          
                          h4(em(strong("Patient Level Results"))),
                          dataTableOutput("patient_table_3"),
                          h4(em(strong("Model Performance Information"))),
                          dataTableOutput("model_table_3"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          sliderInput("thresholdslider_3", "Manually Set Classification Threshold", min=0, max=1,value = 0.04845, width = "600px", round = -4, step = 0.001),
                          actionButton("reset_threshold_3", "Reset to Default Classification Threshold"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          h3(strong("Visual Summary")),
                          tags$head((tags$style("#graph_header_3{font-size:22px; font-weight: bold;}"))),
                          
                          plotlyOutput('int_plot_3')
                          
                        ) # closure to mainPanel(
                      ), # closure to sidebarLayout(    
                      
                      
                      
                      # Footer Equations/Notes
                      hr(),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h3("Calculations & Analytic Background"),
                      tags$head((tags$style("#calculation_header{font-size:22px; font-weight: bold;}"))),
                      
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Logistic regression formulation of the log odds-ratio (or logit):"),
                      tags$head((tags$style("#logit_formula_3_text{font-size:16px;}"))),
                      
                      
                      withMathJax(paste("$$logit(\\hat{p}) = -3.371 + (Seizure at Year-1) \\times 3.067 +
                       (Craniectomy) \\times 0.981 + 
                       (Craniotomy) \\times 0.414 +
                       (Acute Seizures) \\times 0.684 + $$
                       \n
                       $$(Intracranial Fragments) \\times 0.434 +
                       (Preinjury Drug Use) \\times 0.629 + 
                       (Acute LOS) \\times 0.007 $$")),
                      
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Individual predicted probability values then calculated using the following formula:"),
                      tags$head((tags$style("#phat_text{font-size:16px;}"))),
                      
                      withMathJax(paste("$$\\hat{p} = \\frac{1}{1+exp(-logit(\\hat{p}))}$$")),
                      

             ), # end to fourth tabpanel
             
             
             
### Year-2 new seizure incidence prediction model tab
             
             tabPanel("Year-2 New Seizure Incidence Prediction Model", fluid=TRUE, icon=icon("chart-bar"),
                      # Title and Headers
                      h1(strong("Predicting Seizure Events in Moderate-to-Severe Traumatic Brain Injury")),
                      h2(strong("Year-2 New Seizure Model")),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      
                      sidebarLayout(
                        
                        
                        # Sidebar (Input Section)
                        sidebarPanel(
                          # Imaging Inputs (SDH, EDH, CT Fragment)    
                          h2(strong("Input Dashboard")),
                          shinyjs::useShinyjs(),
                          id = "side-panel_4",
                          
                          actionButton("reset_button_yr1_4", "New Patient/Reset to Defaults"),
                          
                          
                          h3(strong("Imaging Data")),
                          # Contusions Dropdown
                          selectInput("contusion_num_4", h3(strong("Number of Contusions:")),
                                      c("0" = 0,
                                        "1" = 1,
                                        "2" = 2,
                                        "3" = 3,
                                        "4+" = 4)),
                          
                          # Hospitalization/Injury Data (injury severity [PTA, TTFC, GCS], Alcohol, Acute Seizures, Craniotomy/Craniectomy)
                          h3(strong("Injury/Hospitalization Information")),
                          textInput("acute_los_4", label = "Acute LOS (days)", value = NULL),
                          checkboxInput("acute_sz_4", label = "Acute Seizures", value = NULL),
                          checkboxInput("craniectomy_4", label="Craniectomy", value=NULL),
                          checkboxInput("craniotomy_4", label="Craniotomy", value=NULL),
                          
                          # Pre-injury medical history information (Incarceration, Psych Hospitalization, Neurodegenerative Disease)
                          h3(strong("Pre-Injury Information")),
                          checkboxInput("drug_4", label = "Pre-Injury Drug Use", value = NULL),
                          
                        ), # closure to sidebarPanel(
                        
                        
                        
                        # mainPanel: output section for tabular and plotted results
                        mainPanel(
                          h3(strong("Tabular Summary")),
                          tags$head((tags$style("#tabular_header{font-size:22px; font-weight: bold;}"))),
                          
                          h4(em(strong("Patient Level Results"))),
                          dataTableOutput("patient_table_4"),
                          h4(em(strong("Model Performance Information"))),
                          dataTableOutput("model_table_4"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          sliderInput("thresholdslider_4", "Manually Set Classification Threshold", min=0, max=1,value = 0.04941, width = "600px", round = -4, step = 0.001),
                          actionButton("reset_threshold_4", "Reset to Default Classification Threshold"),
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          
                          div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                          h3(strong("Visual Summary")),
                          tags$head((tags$style("#graph_header_4{font-size:22px; font-weight: bold;}"))),
                          
                          plotlyOutput('int_plot_4')
                          
                        ) # closure to mainPanel(
                      ), # closure to sidebarLayout(    
                      
                      
                      
                      # Footer Equations/Notes
                      hr(),
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h3("Calculations & Analytic Background"),
                      tags$head((tags$style("#calculation_header{font-size:22px; font-weight: bold;}"))),
                      
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Logistic regression formulation of the log odds-ratio (or logit):"),
                      tags$head((tags$style("#logit_formula_4_text{font-size:16px;}"))),
                      
                      
                      withMathJax(paste("$$logit(\\hat{p}) = -3.801 + (Craniectomy) \\times 1.047 + 
                       (Acute Seizures) \\times 0.892 + 
                       (Contusions (1)) \\times 0.000 + (Contusions (2)) \\times 0.603 +
                       (Contusions (3)) \\times 0.846 + (Contusions (4+)) \\times 0.866 + $$
                       \n
                       $$(Craniotomy) \\times 0.563 +  
                       (Preinjury Drug Use) \\times 0.653 + 
                       (Acute LOS) \\times 0.005 $$")),
                      
                      div(style="vertical-align:top; width: 100px;",HTML("<br>")),
                      h4("Individual predicted probability values then calculated using the following formula:"),
                      tags$head((tags$style("#phat_text{font-size:16px;}"))),
                      
                      withMathJax(paste("$$\\hat{p} = \\frac{1}{1+exp(-logit(\\hat{p}))}$$")),
                      
             ) # end to fifth tabpanel
             
  ) # closure to entire tabsetpanel
) # closure to fluid page(


#//////////////////////////////////////////////////////#
########### Server and Output Information ############# 
#//////////////////////////////////////////////////////#

server= function(input, output) { 
  
  # Table 
  output$patient_table_1 <- renderDataTable({ 
    datatable(data.frame(
      `PTS Prediction` = pred_prob_func_1(input)$classification_1,
      `Predicted Probability` = round(pred_prob_func_1(input)$pred_prob_1, 3),
      Threshold = input$thresholdslider_1,
      `Total Percentile` = round(quantile_fun_1(input=input)$total_percent_1, 4)*100,
      `PTS Percentile` = round(quantile_fun_1(input=input)$pte_percent_1, 4)*100,
      `No PTS Percentile` = round(quantile_fun_1(input=input)$nonpte_percent_1, 4)*100),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}"))) %>% formatStyle(columns = c(1:2), fontWeight = 'bold')
  })

  output$patient_table_2 <- renderDataTable({ 
    datatable(data.frame(
      `PTS Prediction` = pred_prob_func_2(input)$classification_2,
      `Predicted Probability` = round(pred_prob_func_2(input)$pred_prob_2, 3),
      Threshold = input$thresholdslider_2,
      `Total Percentile` = round(quantile_fun_2(input=input)$total_percent_2, 4)*100,
      `PTS Percentile` = round(quantile_fun_2(input=input)$pte_percent_2, 4)*100,
      `No PTS Percentile` = round(quantile_fun_2(input=input)$nonpte_percent_2, 4)*100),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}"))) %>% formatStyle(columns = c(1:2), fontWeight = 'bold')
  })

  output$patient_table_3 <- renderDataTable({ 
    datatable(data.frame(
      `PTS Prediction` = pred_prob_func_3(input)$classification_3,
      `Predicted Probability` = round(pred_prob_func_3(input)$pred_prob_3, 3),
      Threshold = input$thresholdslider_3,
      `Total Percentile` = round(quantile_fun_3(input=input)$total_percent_3, 4)*100,
      `PTS Percentile` = round(quantile_fun_3(input=input)$pte_percent_3, 4)*100,
      `No PTS Percentile` = round(quantile_fun_3(input=input)$nonpte_percent_3, 4)*100),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}"))) %>% formatStyle(columns = c(1:2), fontWeight = 'bold')
  })
  
  output$patient_table_4 <- renderDataTable({ 
    datatable(data.frame(
      `PTS Prediction` = pred_prob_func_4(input)$classification_4,
      `Predicted Probability` = round(pred_prob_func_4(input)$pred_prob_4, 3),
      Threshold = input$thresholdslider_4,
      `Total Percentile` = round(quantile_fun_4(input=input)$total_percent_4, 4)*100,
      `PTS Percentile` = round(quantile_fun_4(input=input)$pte_percent_4, 4)*100,
      `No PTS Percentile` = round(quantile_fun_4(input=input)$nonpte_percent_4, 4)*100),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}"))) %>% formatStyle(columns = c(1:2), fontWeight = 'bold')
  })
  
  
    
  output$model_table_1 <- renderDataTable({ 
    datatable(data.frame(
      Sensitivity = model_per_1(input)$sens_1,
      Specificity = model_per_1(input)$spec_1, 
      PPV = model_per_1(input)$ppv_1,
      NPV = model_per_1(input)$npv_1),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}")))
  })
  
  output$model_table_2 <- renderDataTable({ 
    datatable(data.frame(
      Sensitivity = model_per_2(input)$sens_2,
      Specificity = model_per_2(input)$spec_2, 
      PPV = model_per_2(input)$ppv_2,
      NPV = model_per_2(input)$npv_2),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}")))
  })  
  
  output$model_table_3 <- renderDataTable({ 
    datatable(data.frame(
      Sensitivity = model_per_3(input)$sens_3,
      Specificity = model_per_3(input)$spec_3, 
      PPV = model_per_3(input)$ppv_3,
      NPV = model_per_3(input)$npv_3),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}")))
  })
  
  output$model_table_4 <- renderDataTable({ 
    datatable(data.frame(
      Sensitivity = model_per_4(input)$sens_4,
      Specificity = model_per_4(input)$spec_4, 
      PPV = model_per_4(input)$ppv_4,
      NPV = model_per_4(input)$npv_4),
      rownames = F, 
      options = list(dom = 't', 
                     ordering=F,
                     initComplete = htmlwidgets::JS(
                       "function(settings, json) {",
                       paste0("$(this.api().table().container()).css({'font-size': '", "17px", "'});"),
                       "}")))
  })  
  
  
  observeEvent(input$reset_button_1, {
    shinyjs::reset("side-panel_1")
  })

  observeEvent(input$reset_button_2, {
    shinyjs::reset("side-panel_2")
  })

  observeEvent(input$reset_button_3, {
    shinyjs::reset("side-panel_3")
  })
  
  observeEvent(input$reset_button_4, {
    shinyjs::reset("side-panel_4")
  })

    
    
  observeEvent(input$reset_threshold_1, {
    shinyjs::reset("thresholdslider_1")
  })

  observeEvent(input$reset_threshold_2, {
    shinyjs::reset("thresholdslider_2")
  })
  
  observeEvent(input$reset_threshold_3, {
    shinyjs::reset("thresholdslider_3")
  })

  observeEvent(input$reset_threshold_4, {
    shinyjs::reset("thresholdslider_4")
  })
  
  
  
  # Graphics/Visualization
  
  output$graph_header_1 <- renderText({
    "Visual Summary"
  })

  output$graph_header_2 <- renderText({
    "Visual Summary"
  })
  
  output$graph_header_3 <- renderText({
    "Visual Summary"
  })

  output$graph_header_4 <- renderText({
    "Visual Summary"
  })
  
  
  output$int_plot_1 <- renderPlotly({
    
    phats_graph_1$`Total Percentile` <- unlist(t(sapply(phats_graph_1$Phat, function(x) quantile_fun_1(value=x)))[,1])
    phats_graph_1$`PTS Percentile` <- unlist(t(sapply(phats_graph_1$Phat, function(x) quantile_fun_1(value=x)))[,2])
    phats_graph_1$`No PTS Percentile` <- unlist(t(sapply(phats_graph_1$Phat, function(x) quantile_fun_1(value=x)))[,3])
    
    int_plot_1 <- ggplot(phats_graph_1, aes(x=Phat, fill = `Seizure Event Status`)) + geom_density(alpha=0.5) +
      geom_vline(xintercept = input$thresholdslider_1, linetype = 'dashed') +
      geom_vline(xintercept = pred_prob_func_1(input)$pred_prob_1) +
      xlab('Threshold Percentage') + ylab('Density') +
      theme_minimal() + scale_fill_manual(values=c("#5D3A9B", "#E66100"), name="")
    
    
    ggplotly(int_plot_1, tooltip=c("x", "fill"))
    
  })

  output$int_plot_2 <- renderPlotly({
    
    phats_graph_2$`Total Percentile` <- unlist(t(sapply(phats_graph_2$Phat, function(x) quantile_fun_2(value=x)))[,1])
    phats_graph_2$`PTS Percentile` <- unlist(t(sapply(phats_graph_2$Phat, function(x) quantile_fun_2(value=x)))[,2])
    phats_graph_2$`No PTS Percentile` <- unlist(t(sapply(phats_graph_2$Phat, function(x) quantile_fun_2(value=x)))[,3])
    
    int_plot_2 <- ggplot(phats_graph_2, aes(x=Phat, fill = `Seizure Event Status`)) + geom_density(alpha=0.5) +
      geom_vline(xintercept = input$thresholdslider_2, linetype = 'dashed') +
      geom_vline(xintercept = pred_prob_func_2(input)$pred_prob_2) +
      xlab('Threshold Percentage') + ylab('Density') +
      theme_minimal() + scale_fill_manual(values=c("#5D3A9B", "#E66100"), name="")
    
    
    ggplotly(int_plot_2, tooltip=c("x", "fill"))
    
  })

  output$int_plot_3 <- renderPlotly({
    
    phats_graph_3$`Total Percentile` <- unlist(t(sapply(phats_graph_3$Phat, function(x) quantile_fun_3(value=x)))[,1])
    phats_graph_3$`PTS Percentile` <- unlist(t(sapply(phats_graph_3$Phat, function(x) quantile_fun_3(value=x)))[,2])
    phats_graph_3$`No PTS Percentile` <- unlist(t(sapply(phats_graph_3$Phat, function(x) quantile_fun_3(value=x)))[,3])
    
    int_plot_3 <- ggplot(phats_graph_3, aes(x=Phat, fill = `Seizure Event Status`)) + geom_density(alpha=0.5) +
      geom_vline(xintercept = input$thresholdslider_3, linetype = 'dashed') +
      geom_vline(xintercept = pred_prob_func_3(input)$pred_prob_3) +
      xlab('Threshold Percentage') + ylab('Density') +
      theme_minimal() + scale_fill_manual(values=c("#5D3A9B", "#E66100"), name="")
    
    
    ggplotly(int_plot_3, tooltip=c("x", "fill"))
    
  })
  
  output$int_plot_4 <- renderPlotly({
    
    phats_graph_4$`Total Percentile` <- unlist(t(sapply(phats_graph_4$Phat, function(x) quantile_fun_4(value=x)))[,1])
    phats_graph_4$`PTS Percentile` <- unlist(t(sapply(phats_graph_4$Phat, function(x) quantile_fun_4(value=x)))[,2])
    phats_graph_4$`No PTS Percentile` <- unlist(t(sapply(phats_graph_4$Phat, function(x) quantile_fun_4(value=x)))[,3])
    
    int_plot_4 <- ggplot(phats_graph_4, aes(x=Phat, fill = `Seizure Event Status`)) + geom_density(alpha=0.5) +
      geom_vline(xintercept = input$thresholdslider_4, linetype = 'dashed') +
      geom_vline(xintercept = pred_prob_func_4(input)$pred_prob_4) +
      xlab('Threshold Percentage') + ylab('Density') +
      theme_minimal() + scale_fill_manual(values=c("#5D3A9B", "#E66100"), name="")
    
    
    ggplotly(int_plot_4, tooltip=c("x", "fill"))
    
  })
  
}


shinyApp(ui, server)
