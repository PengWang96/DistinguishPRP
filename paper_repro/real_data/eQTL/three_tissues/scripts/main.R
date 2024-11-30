pvec = c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(pvec, inverse_P_mis)

pvec <- seq(0, 0.05, length = 1000)
k_vec <- sapply(pvec, inverse_P_mis)

rm(list = ls())
library(DistinguishPRP)
library(ggplot2)
library(foreach)
library(doParallel)
library(rstudioapi)
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())

t1 <- Sys.time()
set.seed(234)
data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
head(data)


subsample <- 1:nrow(data)
# subsample <- which(data$Gene == "ENSG00000176945.16")
gene <- data[subsample, 1]
hat_beta <- as.matrix(data[subsample, c(2, 4, 6)])
hat_sigma_sq <- as.matrix(data[subsample, c(3, 5, 7)]^2)

P_mis_vec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(P_mis_vec, inverse_P_mis)
gc()

N <- 10000
r <- 0.05 # burn-in rate
m <- ncol(hat_beta)
num <- nrow(hat_beta)

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

parts <- 1
rows_per_part <- floor(num / parts)
remainder <- num %% parts
results_list <- vector("list", length(P_mis_vec))  # 保存每个P_mis_vec对应的结果


for(p_idx in 1:length(P_mis_vec)) {
  heterogeneity_level <- k_vec[p_idx]
  p_values_per_pmis <- vector("list", parts)  # 每个parts的结果将存储在这个list中

  for(i in 1:parts) {
    start_row <- ((i - 1) * rows_per_part) + 1
    end_row <- i * rows_per_part
    if (i == parts) {
      end_row <- end_row + remainder
    }
    part_hat_beta <- hat_beta[start_row:end_row, ]
    part_hat_sigma_sq <- hat_sigma_sq[start_row:end_row, ]

    p_values_per_pmis[[i]] <- foreach(x = 1:nrow(part_hat_beta), .combine = c, .packages = "DistinguishPRP") %dopar% {
      metropolis_hastings(N, r, m, part_hat_beta[x, ], part_hat_sigma_sq[x, ], heterogeneity_level = heterogeneity_level)$p_value
    }
  }
  results_list[[p_idx]] <- do.call(c, p_values_per_pmis)
}
stopCluster(cl)

p_values_matrix <- do.call(cbind, results_list)
colnames(p_values_matrix) <- paste0("P_mis_", as.character(P_mis_vec))
p_values_df <- as.data.frame(p_values_matrix)
head(p_values_df)
#   xlab("Frequentist p-values") +


result <- as.data.frame(cbind(gene, p_values_df))
save(result, file = "../output/results_PRP_matrix.rda") # without fixed heterogeneity_level
write.table(result, file = "../output/results_PRP_matrix.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
t2 <- Sys.time()
(t2 - t1)
ggplot(result, aes(x = P_mis_0.01)) +
  geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  xlab("Posterior p-values under distinguishability criterion") +
  ylab("Count")





################################################################################
#######  Compare hist of p-values with PRP and eqtlbma  ########################
################################################################################
rm(list = ls())
library(ggplot2)
library(readxl)
library(rstudioapi)
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}

setwd(dirname(current_file))
print(getwd())
load("../output/results_PRP_Pmis0.05(MH-MCMC).rda")
PRP <- result
load("../output/results_Cochran_Q.rda")
Q_pval <- result
load("../output/results_PRP_Pmis0(MH-MCMC).rda")
PRP_k0 <- result
eqtlbma <- read.csv("../eqtlbma_output/three_tissues/data/marginal_p_value.csv", header=T)
# eqtlbma <- as.data.frame(read_excel("../eqtlbma_output/three_tissues/data/output_prefix_avg_bfs.xlsx"))


col_names <- as.character(eqtlbma[1, ]); eqtlbma <- eqtlbma[-1,]; colnames(eqtlbma) <- col_names

index <- match(PRP[, 1], eqtlbma[, 1])
eqtlbma <- eqtlbma[index, ]

if(identical(PRP[, 1], eqtlbma[, 1])) {
  print("The first column of a and b are identical.")
} else {
  print("The first column of a and b are not identical.")
}

if(identical(PRP[, 2], eqtlbma[, 2])) {
  print("The second column of a and b are identical.")
} else {
  print("The second column of a and b are not identical.")
}

eqtlbma[ , 3:ncol(eqtlbma)] <- lapply(eqtlbma[ , 3:ncol(eqtlbma)], as.numeric)



# eqtlbmaaa <- eqtlbma
# eqtlbmaaa$`config.0-0-0` <- 1 - (eqtlbmaaa$`config.1` + eqtlbmaaa$`config.2` + eqtlbmaaa$`config.3` +
#                                    eqtlbmaaa$`config.1-2` + eqtlbmaaa$`config.1-3` + eqtlbmaaa$`config.2-3` +
#                                    eqtlbmaaa$`config.1-2-3`)
# eqtlbmaaa$`config.0-0-0` <- pmin(pmax(eqtlbmaaa$`config.0-0-0`, 0), 1)
#
# eqtlbmaaa <- eqtlbmaaa[, c("gene", "snp", "gene.post", "config.0-0-0",
#                        "config.1", "config.2", "config.3",
#                        "config.1-2", "config.1-3", "config.2-3", "config.1-2-3")]
# head(eqtlbmaaa)
# save(eqtlbmaaa, file = "../output/results_eqtlbma.rda")


eqtlbmaa <- rowSums(eqtlbma[ , 4:(ncol(eqtlbma)-1)], na.rm = TRUE)

eqtlbmaa <- 1 - as.numeric(eqtlbmaa) #################
# eqtlbmaa <- as.numeric(1 - eqtlbma[, ncol(eqtlbma)])


PRPP <- PRP$p_values
2 * mean(PRPP)
sum(PRPP <= 0.05)/length(PRPP)
sum(eqtlbmaa <= 0.05)/length(eqtlbmaa)

cor.test(PRPP, eqtlbmaa)
cor.test(-log10(PRPP), -log10(eqtlbmaa))




data <- data.frame(
  value = c(PRPP, eqtlbmaa),
  Method = factor(c(rep("posterior-PRP", length(PRPP)), rep("eqtlbma", length(eqtlbmaa))),
                    levels = c("posterior-PRP", "eqtlbma"))
)
# data <- data.frame(
#   value = c(eqtlbma, PRP),
#   Method = factor(c(rep("eqtlbma", length(eqtlbma)), rep("posterior-PRP", length(PRP))),
#                   levels = c("eqtlbma", "posterior-PRP"))
# )


ggplot(data.frame(eqtlbmaa), aes(x = eqtlbmaa)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.6) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.00, 1.00))


plot <- ggplot(data.frame(PRPP), aes(x = PRPP)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  scale_x_continuous(limits = c(-0.00, 1)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, hjust = 0.5),
    axis.title.x = element_text(size = 17), # Increase x-axis title font size
    axis.title.y = element_text(size = 17)  # Increase y-axis title font size
  ) +
  xlab("Posterior-PRPs under Distinguishability Criterion") +
  ylab("Count") +
  ylim(0, 1500)
plot
ggsave("../plot/eqtl_pval.pdf", plot)



which(eqtlbmaa > 0.75 & PRPP < 0.001)


library(forestplot)
data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
data <- na.omit(data)
selected_gene <- data[subsample, ]

plot_matrix <- matrix(c(
  selected_gene$bhat_Artery_Aorta, selected_gene$bhat_Liver, selected_gene$bhat_Muscle_Skeletal,
  selected_gene$bhat_Artery_Aorta - 1.96 * selected_gene$se_Artery_Aorta,
  selected_gene$bhat_Liver - 1.96 * selected_gene$se_Liver,
  selected_gene$bhat_Muscle_Skeletal - 1.96 * selected_gene$se_Muscle_Skeletal,
  selected_gene$bhat_Artery_Aorta + 1.96 * selected_gene$se_Artery_Aorta,
  selected_gene$se_Liver + 1.96 * selected_gene$se_Liver,
  selected_gene$bhat_Muscle_Skeletal + 1.96 * selected_gene$se_Muscle_Skeletal
), nrow=3, byrow=TRUE)

labels <- c('Artery Aorta', 'Liver', 'Muscle Skeletal')

forestplot(
  labeltext = labels,
  mean = plot_matrix[1,],
  lower = plot_matrix[2,],
  upper = plot_matrix[3,],
  xlab = "Effect Sizes",
  lwd.ci = 4,
  txt_gp = fpTxtGp(
    label = gpar(cex = 1.2),      # 调整y轴标签字体大小
    ticks = gpar(cex = 1.2),      # 调整x轴刻度标签字体大小
    xlab = gpar(cex = 1.4)        # 调整x轴标签字体大小
  ),
  new_page = FALSE
)
















# Forest plot
rm(list = ls())
library(forestplot)
library(ggplotify)
library(patchwork)
library(rstudioapi)
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}
setwd(dirname(current_file))
print(getwd())
data <- read.table("../data/simple_data.3tissue.summary", header = TRUE)
data <- na.omit(data)
head(data)

# selected_genes <- c('ENSG00000214226.8', 'ENSG00000116205.12',
#                     'ENSG00000166002.6', 'ENSG00000219200.11')
# selected_genes <- data$Gene[which(eqtlbmaa > 0.75 & PRPP < 0.001)[c(2, 13, 1, 3)]]
selected_genes <- data$Gene[which(eqtlbmaa > 0.75 & PRPP < 0.0001)[c(1, 2, 3, 6)]]
plots <- list()
for (i in 1:4) {
  selected_gene <- data[data$Gene == selected_genes[i],]

  plot_matrix <- matrix(c(
    selected_gene$bhat_Artery_Aorta, selected_gene$bhat_Liver, selected_gene$bhat_Muscle_Skeletal,
    selected_gene$bhat_Artery_Aorta - 1.96 * selected_gene$se_Artery_Aorta,
    selected_gene$bhat_Liver - 1.96 * selected_gene$se_Liver,
    selected_gene$bhat_Muscle_Skeletal - 1.96 * selected_gene$se_Muscle_Skeletal,
    selected_gene$bhat_Artery_Aorta + 1.96 * selected_gene$se_Artery_Aorta,
    selected_gene$bhat_Liver + 1.96 * selected_gene$se_Liver,
    selected_gene$bhat_Muscle_Skeletal + 1.96 * selected_gene$se_Muscle_Skeletal
  ), nrow=3, byrow=TRUE)

  # Create labels for the plot
  labels <- c('Artery Aorta', 'Liver', 'Muscle Skeletal')

  p <-
    forestplot(labeltext=labels,
               mean=plot_matrix[1,],
               lower=plot_matrix[2,],
               upper=plot_matrix[3,],
               xlab=paste("Effect Size for Gene", selected_genes[i]),
               lwd.ci = 3, # Make the lines bolder
               txt_gp = fpTxtGp(
                 label = gpar(cex = 1.2),      # 调整y轴标签字体大小
                 ticks = gpar(cex = 1.2),      # 调整x轴刻度标签字体大小
                 xlab = gpar(cex = 1.4)        # 调整x轴标签字体大小
               ),
               boxsize = 0.15, # 调整图形元素大小
               lineheight = unit(1.5, "cm"), # 增加行高
               # col = fpColors(box = "royalblue", line = "darkblue", summary = "royalblue"),
               graphwidth = unit(80, "mm"), # 增加图形部分的宽度
               new_page = F)

  plots[[i]] <- p
}

p1_grob <- grid2grob(print(plots[[1]]))
p2_grob <- grid2grob(print(plots[[2]]))
p3_grob <- grid2grob(print(plots[[3]]))
p4_grob <- grid2grob(print(plots[[4]]))

p_combined <- wrap_elements(p1_grob) + wrap_elements(p2_grob) +
  wrap_elements(p3_grob) + wrap_elements(p4_grob)
p_combined
















################################################################################
####### Compare varying Pmis ###################################################
################################################################################
rm(list = ls())
library(ggplot2)
setwd("D:/R/Replicability/Peng/real data/eQTL/three tissues/output")
load_data <- function(k) {
  load(paste0("reslutQ_k", k, ".rda"))
  data.frame(k = k, p_value = result_all$p_values)
}
k_values <- c(0, 0.2726814, 0.3349839, 0.4626489, 0.5420572, 100)
P_mis_values <- c(0, 0.05, 0.1, 0.2, 0.25, 0.5)
data <- do.call(rbind, lapply(k_values, load_data))

library(dplyr)
groups <- c("P_mis = 0, # of sig. p-value = 8439",
            "P_mis = 0.05, # of sig. p-value = 8286",
            "P_mis = 0.1, # of sig. p-value = 8053",
            "P_mis = 0.2, # of sig. p-value = 7412",
            "P_mis = 0.25, # of sig. p-value = 6913",
            "P_mis = 0.5, # of sig. p-value = 0")
data <- data %>%
  mutate(group = factor(k, levels = k_values, labels = groups))



# # Define a function to calculate the count of p-values less than 0.05
# count_p_less_than_0.05 <- function(df) {
#   sum(df < 0.05)
# }
# counts <- tapply(data$p_value, data$k, count_p_less_than_0.05)

ggplot(data, aes(x = p_value)) +
  geom_histogram(color="white", fill="#6186ad", bins=20, boundary=0, alpha=0.8) +
  facet_wrap(~ group) +
  theme_bw() +
  scale_x_continuous(limits = c(-0.00, 1.00)) +
  labs(x = "p-value", y = "Count", title = "Histogram of p-values for different P_mis")
























rm(list = ls())
library(ggplot2)
library(foreach)
library(doParallel)
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}
library(rstudioapi)

if (!requireNamespace("rprojroot", quietly = TRUE)) {
  install.packages("rprojroot")
}
library(rprojroot)
if (interactive()) {
  current_file <- rstudioapi::getActiveDocumentContext()$path
} else {
  args <- commandArgs(trailingOnly = FALSE)
  current_file <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
}

setwd(dirname(current_file))
print(getwd())

t1 <- Sys.time()
set.seed(234)
source("./my_posterior.R")

data_MS <- read.table("../data/Muscle_Skeletal.out", header = F); data_MS <- na.omit(data_MS)
data_AA <- read.table("../data/Artery_Aorta.out", header = F); data_AA <- na.omit(data_AA)
data_L <- read.table("../data/Liver.out", header = F); data_L <- na.omit(data_L)
head(data_MS)
head(data_AA)
head(data_L)

# 删除第三列或第四列中存在 0 的行
data_MS_filtered <- data_MS[data_MS$V3 != 0 & data_MS$V4 != 0, ]
data_AA_filtered <- data_AA[data_AA$V3 != 0 & data_AA$V4 != 0, ]
data_L_filtered <- data_L[data_L$V3 != 0 & data_L$V4 != 0, ]
# 通过 gene (V1) 和 snp (V2) 找到交集
common_genes <- Reduce(intersect, list(data_MS_filtered$V1, data_AA_filtered$V1, data_L_filtered$V1))
# 保留交集行
data_MS_final <- data_MS_filtered[data_MS_filtered$V1 %in% common_genes, ]
data_AA_final <- data_AA_filtered[data_AA_filtered$V1 %in% common_genes, ]
data_L_final <- data_L_filtered[data_L_filtered$V1 %in% common_genes, ]
data_MS_final <- data_MS_final[, -2]
data_AA_final <- data_AA_final[, -2]
data_L_final <- data_L_final[, -2]
colnames(data_MS_final) <- c("gene", "bhat_Muscle_Skeletal", "se_Muscle_Skeletal")
colnames(data_AA_final) <- c("gene", "bhat_Artery_Aorta", "se_Artery_Aorta")
colnames(data_L_final) <- c("gene", "bhat_Liver", "se_Liver")

# 按 gene 列合并数据框
# data <- Reduce(function(x, y) merge(x, y, by = "gene"),
#                list(data_MS_final, data_AA_final, data_L_final))
data <- Reduce(function(x, y) merge(x, y, by = "gene"),
               list(data_AA_final, data_L_final, data_MS_final))
head(data)
# write.table(data, file = "../data/simple_data.3tissue.summary",
#             sep = "\t", row.names = FALSE, quote = FALSE)

# ind <- which(data$gene == "ENSG00000091262.15")
subsample <- 1:10 # nrow(data)
# subsample <- which(data$Gene == "ENSG00000176945.16")
gene <- data[subsample, 1]
hat_beta <- as.matrix(data[subsample, c(2, 4, 6)])
hat_sigma_sq <- as.matrix(data[subsample, c(3, 5, 7)]^2)
rm(data_MS, data_AA, data_L,
   data_MS_filtered, data_AA_filtered, data_L_filtered,
   data_MS_final, data_AA_final, data_L_final, data)
gc()

pvec <- c(10^seq(-10, log10(0.05), 0.01),0.05)
k_vec <- sapply(pvec, inverse_P_mis)
N <- 10000
r <- 0.05 # burn-in rate
m <- ncol(hat_beta)
num <- nrow(hat_beta)

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

parts <- 1
rows_per_part <- floor(num / parts)
remainder <- num %% parts
results_list <- vector("list", parts)

for(i in 1:parts) {
  start_row <- ((i - 1) * rows_per_part) + 1
  end_row <- i * rows_per_part
  if (i == parts) {
    end_row <- end_row + remainder
  }

  part_hat_beta <- hat_beta[start_row:end_row, ]
  part_hat_sigma_sq <- hat_sigma_sq[start_row:end_row, ]

  sim_results <- foreach(x = 1:nrow(part_hat_beta), .combine = c) %dopar% {
    metropolis_hastings(N, r, m, part_hat_beta[x,], part_hat_sigma_sq[x,], k_vec = k_vec)$p_value #
    # fixed_effect(N, m, part_hat_beta[x,], part_hat_sigma_sq[x,])$p_value # MC in FE
    # frequency_pvalue(part_hat_beta[x,], part_hat_sigma_sq[x,]) # Q test in FE
    # calc_PRP_FE_NumInt(part_hat_beta[x,], part_hat_sigma_sq[x,])$Integral_Value # numerically calc PRP in FE
  }
  results_list[[i]] <- sim_results
}

stopCluster(cl)
p_values <- unlist(results_list)

ggplot(data.frame(p_value = p_values), aes(x = p_value)) +
  geom_histogram(color = "white", fill = "#6186ad", bins = 20, boundary = 0, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  xlab("Posterior p-values under distinguishability criterion") +
  ylab("Count")

result <- data.frame(gene = gene, p_values = p_values)
sum(result$p_values <= 0.05)
# save(result, file = "../output/results_PRP_Pmis0.05(MH-MCMC).rda")
# save(result, file = "../output/results_PRP_Pmis0(MH-MCMC).rda")
# save(result, file = "../output/results_PRP_Pmis0(MC).rda")
# save(result, file = "../output/results_Cochran_Q.rda")
# write.table(result, file = "../output/results_Cochran_Q.txt",
#             sep = "\t", row.names = FALSE, quote = FALSE)
# save(result, file = "../output/results_PRP.rda") # without fixed heterogeneity_level
# write.table(result, file = "../output/results_PRP.txt",
#             sep = "\t", row.names = FALSE, quote = FALSE)
# save(result, file = "../output/results_PRP_Pmis0(NumInt).rda")
t2 <- Sys.time()
(t2 - t1)
a <- result$P_mis_0
b <- result$p_values
plot(a, b)














