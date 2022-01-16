# ------------------------------------------------------------------------------
# synthetic example exploring DM with 1 dimensional data

library(hydroGOF)
library(dplyr)
library(ggplot2)

set.seed(2021)
x <- sin(seq(0,7.5,length.out = 100))+2.5
x2 <- x*0.25

mean(x)
mean(x2)
sd(x)
sd(x2)
cor(x,x2)

plot(x)
points(x2, col='red')

gof(x2,x)
gof(x2/0.25, x)


set.seed(1e6)
x3 <- x2 + rnorm(100,0,0.05)

cor(x,x3)
mean(x3)
prop <- mean(x3)/mean(x)
gof(x3,x)
gof(x3/prop, x)
plot(x)
points(x3, col='red')
points(x3/prop, col='blue')



# grouped example
tbl <- tibble(id = 1:100,
              group = rep(1:5, each=20),
              true = x,
              x2 = x2,
              x3 = x3)

set.seed(1e8)
subtr <- tibble(group= 1:5,
                sub = rnorm(5,0,0.2))

tbl <- left_join(tbl, subtr, by="group") %>% 
    mutate(x4 = x3+sub) %>% 
    select(-sub) %>% 
    group_by(group) %>% 
    mutate(bias_cor_x4 = x4/(mean(x4)/mean(true)),
           bias_cor_x3 = x3/(mean(x3)/mean(true)),
           bias_cor = x3/(mean(x3)/mean(true)),
           bias_cor2 = x3 * (mean(true)/mean(x3)),
           dm = sum(true) * (x3/sum(x3)))

# test coefficient of variance
tbl %>% 
    summarise(cv_true = sd(true) / mean(true),
              cv_x3 = sd(x3) / mean(x3),
              cv_bias = sd(bias_cor_x3) / mean(bias_cor_x3))

# overall gof
cor(tbl$true, tbl$x4) # dynamics
sd(tbl$x4)/sd(tbl$true) # variability
pbias(tbl$x4, tbl$true) # bias

cor(tbl$true, tbl$bias_cor_x4) # dynamics
sd(tbl$bias_cor_x4)/sd(tbl$true) # variability
pbias(tbl$bias_cor_x4, tbl$true) # bias


# group-wise gof
( group_gof <- tbl %>% 
    group_by(group) %>% 
    summarise(dyn_nb = cor(true, x4),
              var_nb = sd(x4)/sd(true),
              bias_nb = pbias(x4,true),
              prop_nb = mean(x4)/mean(true),
              # dyn_bc_x3 = cor(true, bias_cor_x3),
              # var_bc_x3 = sd(bias_cor_x3)/sd(true),
              # bias_bc_x3 = pbias(bias_cor_x3,true),
              # prop_bc_x3 = mean(bias_cor_x3)/mean(true),
              dyn_bc_x4 = cor(true, bias_cor_x4),
              var_bc_x4 = sd(bias_cor_x4)/sd(true),
              bias_bc_x4 = pbias(bias_cor_x4,true),
              prop_bc_x4 = mean(bias_cor_x4)/mean(true))) 


# means
m <- tbl %>% 
    ungroup() %>% 
    mutate(n = 1:n()) %>% 
    group_by(group) %>% 
    summarise(m_x4 = mean(x4),
              m_tr = mean(true),
              start = min(n),
              end = max(n),
              n = n()) 


# ------------------------------------------------------------------------------
# Plot FIGURE 1 A and B


ggplot(tbl) +
    geom_line(aes(id, true, col='1 True function'), size = 2, linetype = 2) +
    geom_point(aes(id, x4, col = '2 Estimate'), size = 1) +
    geom_point(aes(id, bias_cor_x4, col = '3 Bias corrected estimate'), size = 1) +
    geom_vline(xintercept = seq(0,101,length.out = 6), col = "grey") +
    theme_bw() +
    geom_segment(data = m,
                 aes(x=start, xend = end, y=m_tr, yend = m_tr,
                     color = '1 True function'),
                 size = 1) +
    geom_segment(data = m,
                 aes(x=start, xend = end, y=m_x4, yend = m_x4,
                     color = '2 Estimate'),
                 size = 1) +
    scale_color_manual("", values = c("red", "blue", "black")) +
    labs(y = "", x = "") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())

# ggsave("figure 1A.pdf")



### autocorrelation
tbl_acf <- tibble(lag = 0:20,
                  true = acf(tbl$true, lag.max = 20)$acf,
                  x4 = acf(tbl$x4, lag.max = 20)$acf,
                  x4_bc = acf(tbl$bias_cor_x4, lag.max = 20)$acf)

ggplot(tbl_acf) +
    geom_line(aes(lag, true, color = "True function"), size = 1) +
    geom_line(aes(lag, x4, color = "Meta-model estimate"), size = 1) +
    geom_line(aes(lag, x4_bc, color = "Bias corrected estimate"), size = 1) +
    scale_colour_manual(name = "", values = c("black", "blue", "red")) +
    theme_bw() +
    labs(y = "Autocorrelation", x = "Lag")
# ggsave("autocorrelation figure 1B.pdf")


