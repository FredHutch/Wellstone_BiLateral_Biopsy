# BSS
# all analyses exclude 01-0022
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
library(corrr)

pkg_dir <- "/Users/cwon2/CompBio/Wellstone_BiLateral_Biopsy"
load(file.path(pkg_dir, "data", "bs_methyl.rda"))
load(file.path(pkg_dir, "data", "comprehensive_df.rda"))
comprehensive_df <- comprehensive_df %>%
  dplyr::rename_with(~str_remove(., '-logSum'), ends_with('-logSum'))

# which one has 4A161S and 4A161L? exclude 01-0022
#
bss <- bs_methyl %>% 
  dplyr::mutate(`4A161S` = if_else(grepl("4A161S", `haplotype_chr4`),
                                   TRUE, FALSE)) %>% 
  dplyr::mutate(`4A161L` = if_else(grepl("4A161L", `haplotype_chr4`),
                                   TRUE, FALSE)) %>%
  dplyr::filter(! Subject == "01-0022") %>%
  dplyr::mutate(allele=if_else(`use BSSL`, "4qA161L", "4qA161S")) %>%
  dplyr::mutate(allele=factor(allele, levels=c("4qA161S", "4qA161L")))

#
# distribution of BSS, group by use BSSL
#
bss %>%
  group_by(`use BSSL`) %>%
  summarise(mean=mean(`BSS Value`),  min=min(`BSS Value`), max=max(`BSS Value`))


#
# boxplot with dots group by 161S and 161L; no 01-0022
#
bss %>% dplyr::filter(!Subject=="01-0022") %>%
  dplyr::mutate(allele=if_else(`use BSSL`, "4qA161L", "4qA161S")) %>%
  dplyr::mutate(allele=factor(allele, levels=c("4qA161S", "4qA161L"))) %>%
  ggplot(aes(x=allele, y=`BSS Value`)) +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    geom_jitter(width=0.2, alpha=0.3, size=1, color="grey50") +
    theme_minimal() +
    labs(x="", y="% methylation") +
    theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust=1))
ggsave(file=file.path(pkg_dir, "manuscript", "figures", "BSS-boxplot-group-by-161L.pdf"),
       height=2.2, width=1.5)

#
# L and R (pearson and scatter plot, incluing 13-0022)
#  
bss_cor <- bss %>% dplyr::select(Subject, location, `BSS Value`) %>%
  spread(key=location, value=`BSS Value`) %>%
  dplyr::select(-Subject) %>%
  corrr::correlate() %>% dplyr::filter(!is.na(L)) %>% pull(L) 

bss %>% dplyr::select(Subject, location, `BSS Value`) %>%
  spread(key=location, value=`BSS Value`) %>%
  ggplot(aes(x=L, y=R)) +
    geom_point() +
    geom_smooth(method="lm", se=FALSE, linewidth=0.7, color="grey75", alpha=0.3) +
    annotate("text", x=Inf, y=-Inf, color="gray25",
             size=3, hjust=1, vjust=-2,
             label=paste0("Pearson=", format(bss_cor, digit=2))) +
    theme_minimal() +
    labs(x="L (% methylation)", y="R (% methylation)")
ggsave(file=file.path(pkg_dir, "manuscript", "figures", "BSS-left-and-right.pdf"), 
       width=2.5, height=2)    

#
# bss and other measurements
#


# 1. bss vs. STIR status
bss %>% dplyr::select(-Subject) %>%
  left_join(comprehensive_df, by="sample_id") %>% 
  ggplot(aes(x=STIR_status, y=`BSS Value`)) +
    geom_boxplot(width=0.5, outlier.shape=NA) +
    geom_jitter(width=0.3, alpha=0.3, size=1, color="grey50") +
    facet_wrap(~allele) +
    theme_bw() +
    labs(x="", y="% methylation")
ggsave(file=file.path(pkg_dir, "manuscript", "figures", "BSS-STIR-by-allele.pdf"), 
       width=2.5, height=2)  

# 2. Dux4 and other, seperate two allele

## 4qA161S
bss %>% dplyr::select(-Subject) %>%
  dplyr::filter(`allele`=="4qA161S") %>%
  left_join(comprehensive_df, by="sample_id") %>% 
  dplyr::select(`BSS Value`, `DUX4-M6`, `Foot Dorsiflexors`,
                `CSS`, `Fat_Infilt_Percent`, `Cumulative Score`) %>%
  GGally::ggcorr(label = TRUE, label_round=1, 
                 label_size = 3, size=2, hjust=0.9, 
                 legend.position="none", layout.exp=3)

# 4qA161L
bss %>% dplyr::select(-Subject) %>%
  dplyr::filter(`allele`=="4qA161L") %>%
  left_join(comprehensive_df, by="sample_id") %>% 
  dplyr::select(`BSS Value`, `DUX4-M6`, `Foot Dorsiflexors`,
                `CSS`, `Fat_Infilt_Percent`, `Cumulative Score`) %>%
  GGally::ggcorr(label = TRUE, label_round=1, 
                 label_size = 3, size=2, hjust=0.9, 
                 legend.position="none", layout.exp=3)  

# 161S
cor_161S <- bss %>% dplyr::select(-Subject) %>%
  dplyr::filter(`allele`=="4qA161S") %>%
  left_join(comprehensive_df, by="sample_id") %>% 
  dplyr::select(`BSS Value`, `DUX4-M6`, `Foot Dorsiflexors`,
                `CSS`, `Fat_Infilt_Percent`, `Cumulative Score`) %>%
  corrr::correlate() %>% 
  select(term, `BSS Value`)
cor_161L <- bss %>% dplyr::select(-Subject) %>%
  dplyr::filter(`allele`=="4qA161L") %>%
  left_join(comprehensive_df, by="sample_id") %>% 
  dplyr::select(`BSS Value`, `DUX4-M6`, `Foot Dorsiflexors`,
                `CSS`, `Fat_Infilt_Percent`, `Cumulative Score`) %>%
  corrr::correlate() %>% 
  select(term, `BSS Value`) 

cor_161S %>% 
  full_join(cor_161L, by="term", suffix=c("4qA161S", "4qA161L")) %>%
  drop_na() %>%
  rename(`4qA161S`=`BSS Value4qA161S`, `4qA161L`=`BSS Value4qA161L`) %>%
  gather(key=allele, value=cor, -term) %>%
  ggplot(aes(x=allele, y=term)) +
    geom_tile(aes(fill=cor), colour = "grey50") +
    geom_text(aes(label=format(cor, digit=1)), size=3) +
    scale_fill_steps2(low="grey75", mid="grey95", high="grey75", 
                     breaks=seq(-1, 1, length.out=10)) +
    theme_minimal() +
    labs(x="", y="", title="Pearson correlation") +
    theme(legend.position = "none", 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave(file.path(pkg_dir, "manuscript", "figures", "BSS-correlation.pdf"), 
       width=2.8, height=2)
