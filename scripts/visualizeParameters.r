args <- commandArgs(TRUE)
output <- args[1]
files <- args [-1]

library(yaml)
library(ggplot2)
config <- yaml.load_file("config/config.yaml")

nfactors = config$nfactors
my_target = config$my_target
target = unlist(config$target)
names <- config$name_parameter

load("EvoStatus/EvoStatus.RData")
iteration <- evolutionary_log$n.iteration

#################################################################################
############################# Parameters convergence ############################
pdf(output)
par(mfrow=c(nrows, ncols))

if(my_target){
  for(index in 1:nfactors){
    if (my_target) {
      ggplot() +
        geom_line(aes(x = 1:iteration, y = evolutionary_log$convergence[index, 1:iteration]), size = 1, linetype = "solid") +
        ggtitle(names[index]) +
        theme_minimal() + 
        xlab("Iteration") + ylab("Optimal Value") +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 13, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, size = 1.2)
        ) +
        scale_x_continuous(breaks = seq(1, iteration)) +
        geom_hline(yintercept = target[index], color = "black", linetype = "solid", size = 1) +
        geom_point(aes(x = 1:iteration, y = evolutionary_log$convergence[index, 1:iteration]), color = "red", size = 3) +
        geom_path(aes(x = 1:iteration, y = evolutionary_log$convergence[index, 1:iteration]), color = "red", size = 1)
    } else {
      ggplot() +
        geom_line(aes(x = 1:iteration, y = evolutionary_log$convergence[index, 1:iteration]), size = 1, linetype = "solid") + 
        ggtitle(names[index]) +
        theme_minimal() + 
        xlab("Iteration") + ylab("Optimal Value") +
        theme(
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 13, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, size = 1.2)
        ) +
        scale_x_continuous(breaks = seq(1, iteration)) +
        geom_point(aes(x = 1:iteration, y = evolutionary_log$convergence[index, 1:iteration]), color = "red", size = 3) +
        geom_path(aes(x = 1:iteration, y = evolutionary_log$convergence[index, 1:iteration]), color = "red", size = 1)
    }
  }
}
#################################################################################
# track how much the objectvie function improved from frst iteration
plot(evolutionary_log$est_approx_avg - evolutionary_log$est_approx_avg[1],
     type = "o",
     ylab = "Objective Function",
     xlab = "Iteration"
)
title(main = "Evolution of Target Function Over Iterations using kernel regression")
#################################################################################
# track how the variation of objective function over iteration
# # high variance indicates that individuals are in different part o search space and we are still not in a refinement phase
plot(evolutionary_log$est_approx_sd,
     type = "o",
     ylab = "Variance",
     xlab = "Iteration"
)
title(main = "Variance of Objective Function Over Iterations")
#################################################################################
dev.off()
