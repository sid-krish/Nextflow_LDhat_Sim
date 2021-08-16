library(ggplot2)

ldhat_results <- read.csv("ldhat_collected_results_sweep_1.csv")

# ggplot(faithful, aes(x=eruptions, y=waiting)) + 
#   geom_point()


# ggplot(data=ldhat_results, aes(x=factor(genome_size_sim), y=(max_rho - scaled_rho_sim)/scaled_rho_sim)) + 
#   geom_boxplot(aes(fill=factor(sample_size_sim), color=factor(sample_size_sim))) + 
#   facet_wrap(~ scaled_rho_sim, label=label_both, ncol=2) + 
#   geom_hline(yintercept=0)

ggplot(data=ldhat_results, aes(x=factor(max_rho), y=scaled_rho_sim)) +
  geom_boxplot(aes(fill=factor(genome_size_sim), color=factor(genome_size_sim)))
