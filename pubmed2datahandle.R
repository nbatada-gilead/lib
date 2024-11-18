pubmed2datahandle <- function(input_id) {

  library(rvest)
  library(stringr)
  

  cancer_types <- list(
    LAML = c("acute myeloid leukemia"),
    ACC = c("adrenocortical carcinoma"),
    BLCA = c("bladder urothelial carcinoma", "bladder cancer", "bladder"),
    BCC = c("basal cell carcinoma", "squamous cell carcinoma"),
    LGG = c("brain lower grade glioma"),
    BRCA = c("breast invasive carcinoma", "breast cancer", "breast tumor"),
    CESC = c("cervical squamous cell carcinoma", "endocervical adenocarcinoma"),
    CHOL = c("cholangiocarcinoma"),
    LCML = c("chronic myelogenous leukemia"),
    COAD = c("colon adenocarcinoma", "colorectal cancer", "colorectal carcinoma", "colon cancer"),
    ESCA = c("esophageal"),
    LIHC = c("liver hepatocellular carcinoma", "hepatocellular carcinoma", "liver cancer"),
    LUAD = c("lung adenocarcinoma", "lung cancer"),
    LUSC = c("lung squamous cell carcinoma"),
    DLBC = c("lymphoid neoplasm diffuse large b-cell lymphoma"),
    OV = c("ovarian serous cystadenocarcinoma", "ovarian"),
    PAAD = c("pancreatic adenocarcinoma", "pancreatic"),
    PRAD = c("prostate adenocarcinoma", "prostate")
  )
  

  get_cancer_type <- function(text) {
    text <- tolower(text)
    for (acronym in names(cancer_types)) {
      for (name in cancer_types[[acronym]]) {
        if (grepl(name, text)) {
          return(acronym)
        }
      }
    }
    return("UNKNOWN")
  }
  

  classify_unknown <- function(text) {
    text <- tolower(text)
    if (grepl("pancancer|pan-cancer", text)) {
      return("PANCANCER")
    }
    if (grepl("tumor|tumour|carcinoma|cancer", text)) {
      return("CANCER")
    }
    return("UNKNOWN")
  }
  

  parse_pubmed <- function(pubmed_id) {
    url <- paste0("https://pubmed.ncbi.nlm.nih.gov/", pubmed_id, "/")
    page <- read_html(url)
    

    title <- page %>% html_node("h1.heading-title") %>% html_text(trim = TRUE)
    abstract <- page %>% html_node("div.abstract") %>% html_text(trim = TRUE)
    citation_text <- page %>% html_node("span.cit") %>% html_text(trim = TRUE)
    first_author_name <- page %>% html_node("a.full-name") %>% html_text(trim = TRUE)
    
    year <- str_extract(citation_text, "\\b\\d{4}\\b")
    last_name <- toupper(word(first_author_name, -1))
    
    # Determine cancer type
    cancer_type_title <- get_cancer_type(title)
    cancer_type_abstract <- get_cancer_type(abstract)
    cancer_type <- ifelse(cancer_type_abstract != "UNKNOWN", cancer_type_abstract, cancer_type_title)
    
    if (cancer_type == "UNKNOWN") {
      cancer_type <- classify_unknown(abstract)
    }
    

    result <- paste(cancer_type, year, pubmed_id, last_name, sep = "_")
    return(result)
  }
  
  # Function to fetch PubMed ID from GEO ID
  fetch_geo_id <- function(geo_id) {
    geo_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo_id)
    page <- read_html(geo_url)
    pubmed_id <- page %>% html_nodes("a") %>% html_text(trim = TRUE) %>%
      str_subset("^\\d+$") %>% first()
    
    if (!is.na(pubmed_id)) {
      return(parse_pubmed(pubmed_id))
    } else {
      return("No PMID found for GEO ID")
    }
  }
  
  # Main logic
  if (startsWith(input_id, "GSE")) {
    return(fetch_geo_id(input_id))
  } else if (grepl("^\\d+$", input_id)) {
    return(parse_pubmed(input_id))
  } else {
    stop("Invalid input. Please provide a valid PubMed ID or GEO ID.")
  }
}


