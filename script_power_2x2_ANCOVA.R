library(MASS); library(lme4); library(car); library(lmerTest); library(tidyverse)

# Expecting to Teach Group ----
set.seed(1)
teach_size <- 1000 # size of the simulated population                                       
teach_meanvector <- c(94, 51, 56) # means for pre and post tests

teach_correl_matrix <- matrix(c(1, 0.6, 0.6, 
                                0.6, 1, 0.8,
                                0.6, 0.8, 1), # values for the sample correlation matrix
                              ncol = 3)
teach_sds <- c(31, 17, 19) # between subject SDs on each test
D <- diag(teach_sds)

teach_covariance_matrix <- D%*%teach_correl_matrix%*%D

teach_distribution <- mvrnorm(n = teach_size,
                              mu = teach_meanvector, 
                              Sigma = teach_covariance_matrix)
head(teach_distribution)
teach_distribution <- teach_distribution %>% as_tibble() %>% 
  mutate(group = "teach") %>%
  rename(pre_test=V1, low_pressure=V2, high_pressure=V3)

head(teach_distribution)

# Expecting to Test Group ----
set.seed(1)
test_size <- 1000 # size of the simulated population                                       
test_meanvector <- c(92, 60, 60) # means for and post tests

test_correl_matrix <- matrix(c(1, 0.6, 0.6, 
                               0.6, 1, 0.8,
                               0.6, 0.8, 1), # values for the sample correlation matrix
                             ncol = 3)
test_sds <- c(35, 23, 22) # between subject SDs on each test
D <- diag(test_sds)

test_covariance_matrix <- D%*%test_correl_matrix%*%D

test_distribution <- mvrnorm(n = test_size,
                             mu = test_meanvector, 
                             Sigma = test_covariance_matrix)
head(test_distribution)
test_distribution <- test_distribution %>% as_tibble() %>% 
  mutate(group = "test") %>%
  rename(pre_test=V1, low_pressure=V2, high_pressure=V3)


# Using the aov() function versus the lmer() function ----
# Below, I fit the model using the lmer() function, 
# but you could build a similar model using aov(), t.test(), lm(), etc.
# I will use lmer() for it's flexibility across designs, but you could 
# replace this with whatever analysis best matches your method. 
test_sample <- sample_n(test_distribution, size=10)
teach_sample <- sample_n(teach_distribution, size=10)
dat <- rbind(test_sample, teach_sample)
dat_long <- dat %>% rownames_to_column(var="subID") %>%
  pivot_longer(cols = low_pressure:high_pressure, 
               names_to = "test", values_to = "score")
head(dat_long)

mod0 <- lmer(score~pre_test+group*test+(1|subID), data=dat_long, REML=TRUE)
Anova(mod0, type="III")
anova(mod0)
anova(mod0)[[6]][1] # covariate p-value
anova(mod0)[[6]][2] # group p-value
anova(mod0)[[6]][3] # test p-value
anova(mod0)[[6]][4] # GxT p-value


dat_long %>% pivot_wider(names_from = "test", values_from = "score") %>%
  mutate(diff = high_pressure-low_pressure) %>%
  group_by(group) %>%
  summarize(pre_ave=mean(pre_test), 
            pre_sd=sd(pre_test),
            low_ave=mean(low_pressure), 
            low_sd=sd(low_pressure), 
            high_ave=mean(high_pressure), 
            high_sd=sd(high_pressure),
            diff_ave = mean(diff),
            diff_sd = sd(diff)) %>%
  pivot_wider(values_from=pre_ave:diff_sd, names_from = group)


# Re-sampling to Estimate Statistical Power ----
ns = c(seq(from=30, to=100, by=10)) # sample sizes for the simulations

k = 200 # number of iterations at each sample size



sim_data <- data.frame(matrix(NA, nrow = length(ns)*k, 
                              ncol = 22))
sim_data <- sim_data %>% rename(n=X1, 
                                iter=X2, 
                                Cov_p=X3, 
                                Group_p=X4, 
                                Test_p=X5, 
                                GxT_p=X6,
                                teach_pre_ave=X7, 
                                teach_pre_sd=X8,
                                teach_low_ave=X9, 
                                teach_low_sd=X10, 
                                teach_high_ave=X11, 
                                teach_high_sd=X12,
                                teach_diff_ave = X13,
                                teach_diff_sd = X14,
                                test_pre_ave=X15, 
                                test_pre_sd=X16, 
                                test_low_ave=X17, 
                                test_low_sd=X18, 
                                test_high_ave=X19, 
                                test_high_sd=X20,
                                test_diff_ave = X21,
                                test_diff_sd = X22)



row = 0
set.seed(1)
for (n in ns) {
  print(n)
  
  for (i in 1:k) {
    row = row + 1
    print(row)
    
    test_sample <- sample_n(test_distribution, size=n)
    teach_sample <- sample_n(teach_distribution, size=n)
    
    dat <- rbind(test_sample, teach_sample)
    dat_long <- dat %>% rownames_to_column(var="subID") %>%
      pivot_longer(cols = low_pressure:high_pressure, 
                   names_to = "test", values_to = "score")
    
    options(contrasts=c("contr.sum","contr.poly"))
    mod0 <- lmer(score~pre_test+group*test+(1|subID), data=dat_long, REML=TRUE)
    sim_data[row, "n"] <- n
    sim_data[row, "iter"] <- i
    
    sim_data[row, "Cov_p"] <- anova(mod0)[[6]][1] # covariate p-value
    sim_data[row, "Group_p"] <- anova(mod0)[[6]][2] # group p-value
    sim_data[row, "Test_p"] <- anova(mod0)[[6]][3] # test p-value
    sim_data[row, "GxT_p"] <- anova(mod0)[[6]][4] # GxT p-value
    
    # Summary Stats
    SUM_DATA<-dat_long %>% pivot_wider(names_from = "test", values_from = "score") %>%
      mutate(diff = high_pressure-low_pressure) %>%
      group_by(group) %>%
      summarize(pre_ave=mean(pre_test), 
                pre_sd=sd(pre_test),
                low_ave=mean(low_pressure), 
                low_sd=sd(low_pressure), 
                high_ave=mean(high_pressure), 
                high_sd=sd(high_pressure),
                diff_ave = mean(diff),
                diff_sd = sd(diff)) %>%
      pivot_wider(values_from=pre_ave:diff_sd, names_from = group)
    
    sim_data[row, "teach_pre_ave"] <- SUM_DATA[, "pre_ave_teach"]
    sim_data[row, "teach_pre_sd"] <- SUM_DATA[, "pre_sd_teach"]
    sim_data[row, "teach_low_ave"] <- SUM_DATA[, "low_ave_teach"]
    sim_data[row, "teach_low_sd"] <- SUM_DATA[, "low_sd_teach"]
    sim_data[row, "teach_high_ave"] <- SUM_DATA[, "high_ave_teach"]
    sim_data[row, "teach_high_sd"] <- SUM_DATA[, "high_sd_teach"]
    sim_data[row, "teach_diff_ave"] <- SUM_DATA[, "diff_ave_teach"]
    sim_data[row, "teach_diff_sd"] <- SUM_DATA[, "diff_sd_teach"]
    
    sim_data[row, "test_pre_ave"] <- SUM_DATA[, "pre_ave_test"]
    sim_data[row, "test_pre_sd"] <- SUM_DATA[, "pre_sd_test"]
    sim_data[row, "test_low_ave"] <- SUM_DATA[, "low_ave_test"]
    sim_data[row, "test_low_sd"] <- SUM_DATA[, "low_sd_test"]
    sim_data[row, "test_high_ave"] <- SUM_DATA[, "high_ave_test"]
    sim_data[row, "test_high_sd"] <- SUM_DATA[, "high_sd_test"]
    sim_data[row, "test_diff_ave"] <- SUM_DATA[, "diff_ave_test"]
    sim_data[row, "test_diff_sd"] <- SUM_DATA[, "diff_sd_test"]
    
    }
}


sim_summary <- sim_data %>% group_by(n) %>%
  summarize(Cov_sig = sum(Cov_p<0.05)/n(),
            Group_sig = sum(Group_p<0.05)/n(),
            Test_sig = sum(Test_p<0.05)/n(),
            GxT_sig = sum(GxT_p<0.05)/n(),) %>%
  pivot_longer(cols=Cov_sig:GxT_sig, names_to = "Effect", values_to = "Power") %>%
  mutate(Effect = fct_relevel(Effect, "Cov_sig", "Group_sig", "Test_sig", "GxT_sig"))
sim_summary

# Color blind friendly palette 
cbPalette <- c("#000000", "#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(data=sim_summary, aes(x=n, y=Power)) +
  geom_line(aes(lty=Effect, col=Effect)) +
  geom_point(aes(shape=Effect, col=Effect)) +
  scale_x_continuous(name = "n/Group") +
  scale_y_continuous(name = "Statistical Power", limits=c(0,1),
                     breaks=c(seq(from=0, to=1, by=0.2))) +
  theme_bw()+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")



head(sim_data)
sim_data_long <- sim_data %>% select(-Cov_p, -Group_p, -Test_p, -GxT_p) %>%
  pivot_longer(cols=teach_pre_ave:test_diff_sd, 
               names_to=c("Group", "Test", "Variable"),
               names_sep="_",
               values_to = "Score") %>%
  pivot_wider(names_from = "Variable", values_from = "Score")

head(sim_data_long)

ggplot(data=sim_data_long %>% filter(Test == "diff"), 
       aes(x=ave)) +
  geom_histogram(aes(fill=Group), col="black")+
  facet_wrap(~n)+
  scale_x_continuous(name = "Difference (High - Low Pressure)") +
  scale_y_continuous(name = "Count") +
  theme_bw()+
  scale_fill_manual(values=cbPalette)+
  scale_colour_manual(values=cbPalette)+
  theme(axis.text=element_text(size=12, color="black"), 
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=12, face="bold"),
        plot.title=element_text(size=12, face="bold", hjust=0.5),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        legend.position = "bottom")
