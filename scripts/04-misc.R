# simulation distribution of the STIR+/+ (1,1), STIR-/- (0,0), and STIR-/+ (0,1) paire 


# 0 = STIR-
# 1 = STIR+

n_simulation <- 10000
sim_list <- vector("list", n_simulation)

for (j in 1:n_simulation) {
  n_pos <- 25 # 1: STIR+
  n_neg <- 43 # 0: STIR-
  n_pairs <- (n_pos + n_neg) / 2

  set.seed <- j
 
  res <- vector("integer", n_pairs)
  for (i in 1:n_pairs) {
    stir <- c(rep(1, n_pos), rep(0, n_neg))
    res[i] <- sum(sample(x=stir, 2, replace=FALSE))
    n_pos <- n_pos - res[i]
     n_neg <- n_neg - (2-res[i])
  }

  sim_list[[j]] <- table(res)
}

sim_res <- do.call(rbind, sim_list) %>%
  as.data.frame() %>%
  rename(`STIR-/-`= `0`, `STIR+/-`=`1`, `STIR+/+`=`2`) %>%
  gather(key="pairs", value="n")
annot <- sim_res %>% group_by(pairs) %>%
  summarise(mu = mean(n)) %>%
  dplyr::mutate(x_mu = mu, y_mu = 0) %>%
  dplyr::mutate(x_n=c(19, 5, 10), y_n=0) %>%
  dplyr::mutate(mu_label = format(mu, digit=2))


 ggplot(sim_res, aes(x=n)) +
   geom_density(adjust=3.5, n=256) +
   #geom_histogram(bins=10) +
   facet_wrap( ~ pairs, nrow=3) +
   geom_point(data=annot, aes(x=x_n, y=y_n), size=2, color="red", shape=4) +
   geom_text(data=annot, aes(x=x_mu, y=y_mu, label=paste0("mu = ", mu_label)), 
             hjust=0, vjust=0, color="gray50") +
   geom_vline(data=annot, aes(xintercept=mu),linetype="dashed", color="gray50") +
   labs(x="number of pairs", y="density") +
   theme_bw()
ggsave(file.path(pkg_dir, "figures", "STIR-2000-simulation-pairs.pdf"), 
       width=4, height=4)   