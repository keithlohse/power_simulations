library(MASS); library(lme4); library(car); library(lmerTest); library(tidyverse)

# Expecting to Teach Group ----
set.seed(1)
size <- 1000 # size of the simulated population                                       
meanvector <- c(51, 56) # means for pre and post tests

correl_matrix <- matrix(c(1, 0.8,
                                0.8, 1), # values for the sample correlation matrix
                              ncol = 2)
sds <- c(17, 19) # between subject SDs on each test
D <- diag(sds)

covariance_matrix <- D%*%correl_matrix%*%D

distribution <- mvrnorm(n = size,
                              mu = meanvector, 
                              Sigma = covariance_matrix)
summary(distribution)

head(distribution)
distribution <- distribution %>% as_tibble() %>% 
  mutate(group = "teach") %>%
  rename(low_pressure=V1, high_pressure=V2)

hist(distribution$low_pressure, xlim=c(0,150))
hist(distribution$high_pressure, xlim=c(0,150))


# Using the aov() function versus the lmer() function ----
# Below, I fit the model using the lmer() function, 
# but you could build a similar model using aov(), t.test(), lm(), etc.
# I will use lmer() for it's flexibility across designs, but you could 
# replace this with whatever analysis best matches your method. 

set.seed(1)
teach_sample <- sample_n(distribution, size=10)

dat_long <- teach_sample %>% rownames_to_column(var="subID") %>%
  pivot_longer(cols = low_pressure:high_pressure, 
               names_to = "test", values_to = "score")
head(dat_long)

mod0 <- lmer(score~test+(1|subID), data=dat_long, REML=TRUE)
Anova(mod0, type="III")
anova(mod0)
anova(mod0)[[6]][1] # test p-value

modelAOV <- aov(score~test+Error(subID), data = dat_long)
summary(modelAOV)
summary(modelAOV)[[2]][[1]][1,5]


mod_t_test <- t.test(score~test, data=dat_long, paired=TRUE)
print(mod_t_test)
mod_t_test$p.value


# Taking the descriptive statistics we might want from each sample
head(dat_long)
dat_long %>% pivot_wider(names_from = "test", values_from = "score") %>%
  mutate(diff = high_pressure-low_pressure) %>%
  summarize(low_ave = mean(low_pressure),
            low_sd = sd(low_pressure),
            high_ave = mean(high_pressure),
            high_sd = sd(high_pressure),
            diff_ave = mean(diff),
            diff_sd = sd(diff))


# Re-sampling to Estimate Statistical Power ----
ns = c(seq(from=10, to=100, by=10)) # sample sizes for the simulations

k = 200 # number of iterations at each sample size



sim_data <- data.frame(matrix(NA, nrow = length(ns)*k, 
                              ncol = 9))
sim_data <- sim_data %>% rename(n=X1, 
                                iter=X2, 
                                WS_p = X3,
                                low_ave=X4, 
                                low_sd=X5, 
                                high_ave=X6, 
                                high_sd=X7,
                                diff_ave = X8,
                                diff_sd = X9)

row = 0
set.seed(1)
for (n in ns) {
  print(n)
  
  for (i in 1:k) {
    row = row + 1
    print(row)
    
    teach_sample <- sample_n(distribution, size=n)
    
    dat_long <- teach_sample %>% rownames_to_column(var="subID") %>%
      pivot_longer(cols = low_pressure:high_pressure, 
                   names_to = "test", values_to = "score")
    
    options(contrasts=c("contr.sum","contr.poly"))
    mod0 <- lmer(score~test+(1|subID), data=dat_long, REML=TRUE)
    
    sim_data[row, 1] <- n
    sim_data[row, 2] <- i
    sim_data[row, 3] <- anova(mod0)[[6]][1] # WS p-value
    
    # Summary Statistics
    SUM_DATA <- dat_long %>% pivot_wider(names_from = "test", values_from = "score") %>%
      mutate(diff = high_pressure-low_pressure) %>%
      summarize(low_ave = mean(low_pressure),
                low_sd = sd(low_pressure),
                high_ave = mean(high_pressure),
                high_sd = sd(high_pressure),
                diff_ave = mean(diff),
                diff_sd = sd(diff))
    
    
    sim_data[row, 4] <- SUM_DATA[, "low_ave"]
    sim_data[row, 5] <- SUM_DATA[, "low_sd"]
    
    sim_data[row, 6] <- SUM_DATA[, "high_ave"]
    sim_data[row, 7] <- SUM_DATA[, "high_sd"]
  
    sim_data[row, 8] <- SUM_DATA[, "diff_ave"]
    sim_data[row, 9] <- SUM_DATA[, "diff_sd"]
    
    }
}

head(sim_data)
sim_summary <- sim_data %>% group_by(n) %>%
  summarize(WS_p = sum(WS_p<0.05)/n(),
            ave_dz = mean(diff_ave/diff_sd),
            ave_d = mean(diff_ave/sqrt(((low_ave^2)+(high_sd^2))/2))) %>%
  pivot_longer(cols=WS_p, names_to = "Effect", values_to = "Power")
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


write.csv(sim_data, "./rm_simulation.csv")
