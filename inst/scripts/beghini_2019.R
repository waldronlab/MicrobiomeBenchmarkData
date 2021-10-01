
# Create phyloseq object

nychanes_otu_table_url <- "https://github.com/waldronlab/nychanesmicrobiome/raw/master/inst/extdata/otu_table_mc10_w_tax.biom"
nychanes_otu_table <- tempfile()
download.file(nychanes_otu_table_url, nychanes_otu_table)

nychanes_taxa_tree_url <- "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.tre"
nychanes_taxa_tree <- tempfile()
download.file(nychanes_taxa_tree_url, nychanes_taxa_tree)

metadata <- sas7bdat::read.sas7bdat("http://nychanes.org/wp-content/uploads/sites/6/2018/01/public_v2_010518.sas7bdat") %>%
  tibble::as_tibble()
original_map <- readr::read_tsv("https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/original_map.tsv")
sample_selection <- readr::read_tsv("https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/smokingsampleselection.tsv")
colnames(sample_selection) <- c("key", "alternative_only", "never", "former", "secondhand", "cigarette")
sample_metadata <- sample_selection %>%
  dplyr::select(-1) %>%
  tidyr::pivot_longer(cols = tidyselect::everything(), names_to = "smoking_status", values_to = "status") %>%
  dplyr::filter(status == TRUE) %>%
  dplyr::select("smoking_status") %>%
  dplyr::bind_cols(sample_selection) %>%
  dplyr::full_join(original_map, by = c("key" = "KEY")) %>%
  dplyr::left_join(metadata, by = c("key" = "KEY")) %>%
  dplyr::relocate(Burklab_ID, key) %>%
  # magrittr::set_colnames(tolower(colnames(.))) %>%
  tibble::column_to_rownames(var = "Burklab_ID") %>%
  as.data.frame()

remove_cols <- c(
  ## controls
  "20151013CZTME1","20151020CZE3","20151020CZE4","20151020TME1","NC1","NC2",
  ## duplicates
  "NYDH0036","NYDH0051","NYDH0060","NYDH0152","NYDH0213","NYDH0487","NYDH0492","NYDH0522","NYDH0527","NYDH0545R","NYDH0649c", "NYDH0661", "NYDH0691", "NYDH0893","NYDH0931","NYDH0988","NYDH1042","NYDH1460","NYDH1353"
)

ps <- phyloseq::import_biom(BIOMfilename = nychanes_otu_table, treefilename = nychanes_taxa_tree)
colnames(phyloseq::tax_table(ps)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
ps <- phyloseq::prune_samples(!phyloseq::sample_names(ps) %in% remove_cols, ps)
phyloseq::sample_names(ps) <- gsub('c|R', '', phyloseq::sample_names(ps))
phyloseq::sample_data(ps) <- sample_metadata
ps <- phyloseq::prune_samples(phyloseq::sample_sums(ps) > 1000, ps)
ps <- phyloseq::subset_taxa(ps, !Class %in% c("D_2__Chloroplast") & !Family %in% c("D_4__Mitochondria"))
phyloseq::tax_table(ps)[,"Genus"] <- sub(" [1-9]", "", phyloseq::tax_table(ps)[,"Genus"])
phyloseq::tax_table(ps)[, "Phylum"] <- sub("^.+__", "", as.character(phyloseq::tax_table(ps)[,"Phylum"]))
phyloseq::tax_table(ps)[, "Genus"] <- sub("^.+__", "", as.character(phyloseq::tax_table(ps)[,"Genus"]))


# Annotate factors

phyloseq::sample_data(ps)$SMOKER3CAT <- factor(phyloseq::sample_data(ps)$SMOKER3CAT,
                                               labels  = c("Never smoker",
                                                           "Current smoker",
                                                           "Former smoker"))
phyloseq::sample_data(ps)$COTININE <- as.numeric(phyloseq::sample_data(ps)$COTININE)
phyloseq::sample_data(ps)$COTININECAT <- factor(phyloseq::sample_data(ps)$COTININECAT, labels = c("Active smoker",
                                                                                                  "Non-smoker + secondhand smoke",
                                                                                                  "Non-smoker no secondhand smoke",
                                                                                                  "Skipped"))
phyloseq::sample_data(ps)$RACE <- factor(phyloseq::sample_data(ps)$RACE, labels = c("Non-Hispanic White",
                                                                                    "Non-Hispanic Black",
                                                                                    "Hispanic",
                                                                                    "Asian",
                                                                                    "Other"))
phyloseq::sample_data(ps)$GENDER <- factor(phyloseq::sample_data(ps)$GENDER, labels = c("Male","Female"))
phyloseq::sample_data(ps)$AGEGRP5C <- factor(phyloseq::sample_data(ps)$AGEGRP5C, labels = c("20-29",
                                                                                            "30-39",
                                                                                            "40-49",
                                                                                            "50-59",
                                                                                            "60 AND OVER"))
phyloseq::sample_data(ps)$SR_ACTIVE <- factor(phyloseq::sample_data(ps)$SR_ACTIVE, labels = c("Very active",
                                                                                              "Somewhat active",
                                                                                              "Not very active/not active at all"))
phyloseq::sample_data(ps)$EDU4CAT <- relevel(factor(phyloseq::sample_data(ps)$EDU4CAT, labels = c( "Less than High school diploma",
                                                                                                   "High school graduate/GED",
                                                                                                   "Some College or associate's degree",
                                                                                                   "College graduate or more")), "College graduate or more")
phyloseq::sample_data(ps)$OHQ_1 <- factor(phyloseq::sample_data(ps)$OHQ_1, levels=1:7,
                                          labels=c("6 mos or less",
                                                   ">6mos, <= 1 yr",
                                                   ">1 yr, <= 2 yr",
                                                   ">2 yr, <= 3 yr",
                                                   ">3 yr, <=5 yr",
                                                   ">5 yr", "never"))
phyloseq::sample_data(ps)$OHQ_2 <- factor(phyloseq::sample_data(ps)$OHQ_2, levels=1:6,
                                          labels = c("More than once a day",
                                                     "Once a day",
                                                     "Every few days",
                                                     "Every few weeks",
                                                     "Never",
                                                     "No teeth or dentures"))
phyloseq::sample_data(ps)$OHQ_3 <- factor(phyloseq::sample_data(ps)$OHQ_3, labels = c("Yes","No",NA))
phyloseq::sample_data(ps)$OHQ_5 <- as.numeric(phyloseq::sample_data(ps)$OHQ_5)
phyloseq::sample_data(ps)$OHQ_5_3CAT <- cut(phyloseq::sample_data(ps)$OHQ_5, breaks=c(-1,0,5,7),
                                            labels = c("0","1-5","6-7"), include.lowest = FALSE)
phyloseq::sample_data(ps)$DBTS_NEW <- factor(phyloseq::sample_data(ps)$DBTS_NEW, labels = c("Yes",
                                                                                            "No",NA))
phyloseq::sample_data(ps)$smokingstatus <- factor(phyloseq::sample_data(ps)$smokingstatus, levels = c("cigarette",
                                                                                                      "never",
                                                                                                      "former",
                                                                                                      "alternativeonly",
                                                                                                      "secondhand"),
                                                  labels = c("Cigarette","Never smoker","Former smoker","Alternative smoker","Secondhand"))
phyloseq::sample_data(ps)$CURRENT_SMKER <- factor(phyloseq::sample_data(ps)$CURRENT_SMKER, labels = c("Yes","No"))
phyloseq::sample_data(ps)$CIGARETTES <- factor(phyloseq::sample_data(ps)$CIGARETTES, labels = c("Yes","No",NA))
phyloseq::sample_data(ps)$BMI <- as.double(phyloseq::sample_data(ps)$BMI)
phyloseq::sample_data(ps)$SPAGE <- as.integer(phyloseq::sample_data(ps)$SPAGE)
phyloseq::sample_data(ps)$MERCURYU <- as.double(phyloseq::sample_data(ps)$MERCURYU )
phyloseq::sample_data(ps)$SMOKER <- factor(phyloseq::sample_data(ps)$SMOKER, labels = c("Smoker",
                                                                                        "Non smoker"))
phyloseq::sample_data(ps)$SMQ_13_1_1 <- factor(phyloseq::sample_data(ps)$SMQ_13_1_1, labels = c("Cigarettes",
                                                                                                "Cigars/Cigarillos",
                                                                                                #"Chewing tobacco",
                                                                                                #"Snuff",
                                                                                                "Hookah Pipe",
                                                                                                "E-Cigsarettes",
                                                                                                #"Nicotine patches/Gum/Other",
                                                                                                "NA"))
phyloseq::sample_data(ps)$SMQ_7 <- as.integer(phyloseq::sample_data(ps)$SMQ_7)
phyloseq::sample_data(ps)$INC25KMOD <- factor(phyloseq::sample_data(ps)$INC25KMOD, labels = c("Less Than $20,000",
                                                                                              "$20,000-$49,999",
                                                                                              "$50,000-$74,999",
                                                                                              "$75,000-$99,999",
                                                                                              "$100,000 or More",
                                                                                              "NA"))
phyloseq::sample_data(ps)$POVGROUP6_0812CT <- factor(phyloseq::sample_data(ps)$POVGROUP6_0812CT, labels = c("0 to 5%",
                                                                                                            "5 to 10%",
                                                                                                            "10 to 20%",
                                                                                                            "20 to 30%",
                                                                                                            "30 to 40%",
                                                                                                            "40 to 100%"))
phyloseq::sample_data(ps)$SMQ_4UNIT <- factor(phyloseq::sample_data(ps)$SMQ_4UNIT, labels = c(1,4,52,NA))
phyloseq::sample_data(ps)$US_BORN <- factor(phyloseq::sample_data(ps)$US_BORN,
                                            levels=c("1","2","NA"),
                                            labels = c("US-Born, 50 States, DC, PR and Territories",
                                                       "Other", NA))
phyloseq::sample_data(ps)$DMQ_7YEAR <- as.integer(phyloseq::sample_data(ps)$DMQ_7YEAR)
phyloseq::sample_data(ps)$DMQ_2 <- factor(phyloseq::sample_data(ps)$DMQ_2,
                                          levels=1:6,
                                          labels=c("Married","Widowed","Divorced",
                                                   "Separated","Never married",
                                                   "Living with partner"))
phyloseq::sample_data(ps)$AGEGRP3C <- factor(phyloseq::sample_data(ps)$AGEGRP3C, levels=1:3,
                                             labels=c("20-34", "35-64", "65 and over"))
phyloseq::sample_data(ps)$EDU3CAT <- relevel(factor(phyloseq::sample_data(ps)$EDU3CAT, levels=1:3,
                                                    labels=c("High School Diploma or Less", "Some College or Associate's Degree",
                                                             "College Graduate or More")),"College Graduate or More")
phyloseq::sample_data(ps)$INC10K <- as.numeric(as.character(phyloseq::sample_data(ps)$INC10K))
phyloseq::sample_data(ps)$INC3C <- relevel(cut(phyloseq::sample_data(ps)$INC10K,
                                               breaks = stats::quantile(phyloseq::sample_data(ps)$INC10K, probs=c(0, 1/3, 2/3, 1), na.rm=TRUE),
                                               include.lowest = TRUE,
                                               labels=c("Less Than $30,000", "$30,000 - $60,000", "$60,000 or more")), "$60,000 or more")
phyloseq::sample_data(ps)$A1C <- as.integer(phyloseq::sample_data(ps)$A1C)
phyloseq::sample_data(ps)$GLUCOSE <- as.integer(phyloseq::sample_data(ps)$GLUCOSE)

format_dbq <- function(hanes_df, n) {
  #function to format DBQ (dietary behavior questionnaire) questions
  varname <- paste0("DBQ_", n)
  var <- as.numeric(hanes_df[[varname]])

  #For unit variables, 1=per day, 2=per week, 3=per month
  #this is to get number per week, so if the unit is 1 (per day), multiply by 7;
  #2 (week), by 1; 3 (month), by 1/4.333333 = 3/13
  mult <- c(7,1,3/13)[hanes_df[[paste0(varname, "UNIT")]]]
  mult[is.na(mult)] <- 0

  var * mult
}

for(i in 3:10) phyloseq::sample_data(ps)[[paste0("DBQ_",i)]] <- format_dbq(phyloseq::sample_data(ps), i)
phyloseq::sample_data(ps)$DBQ_10_3CAT <- cut(phyloseq::sample_data(ps)$DBQ_10, breaks=c(-1,.999,5,999),
                                             labels = c("0-<1","1-5","6 or more"))


colnames(phyloseq::sample_data(ps)) <- tolower(colnames(phyloseq::sample_data(ps)))

ps_subset <- phyloseq::subset_samples(ps, smoking_status %in% c("never", "cigarette"))
ps_subset <- phyloseq::prune_taxa(phyloseq::taxa_sums(x) > 0, x) # same number of taxa


beghini_count_matrix <- as.matrix(as.data.frame(phyloseq::otu_table(ps_subset)))
