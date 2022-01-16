# ------------------------------------------------------------------------------
# Upper Bhima Basin, downscaling runoff from 30min CWATM results to 5min


library(dasymetric) # install devtools::install_github("mkkallio/dasymetric")
library(lubridate)
library(gstat)
library(raster)
library(dplyr)
library(units)
library(sf)
library(imputeTS)
library(exactextractr)
library(ggplot2)
library(patchwork)
library(spdep)
source("R/kge_non_parametric.R")

# ------------------------------------------------------------------------------
# load source and target zone data
source <- readRDS("data/source_zones.RDS")
target <- readRDS("data/target_zones.RDS")



# ------------------------------------------------------------------------------
# Train meta-models

# extract data from source zones for model training
formula <- runoff ~ .
data <- create_model_frame(source, formula)

# create cumulative precipitation sums for the timeseries
temp <- L(data$precipitation, 1:3)
temp <- cbind(temp, zoo::rollsum(data$precipitation, 7, align = "right", na.pad = TRUE),
              zoo::rollsum(data$precipitation, 30, align = "right", na.pad = TRUE),
              zoo::rollsum(data$precipitation, 60, align = "right", na.pad = TRUE),
              zoo::rollsum(data$precipitation, 180, align = "right", na.pad = TRUE))
colnames(temp) <- c("p1","p2","p3", "pc7", "pc30", "pc60", "pc180")
data <- cbind(data, temp)

# create cumulative temperature means for the timeseries
temp <- L(data$Tavg, 1:2)
temp <- cbind(temp, zoo::rollmean(data$Tavg, 7, align = "right", na.pad = TRUE))
colnames(temp) <- c("t1","t2", "tc7")
data <- cbind(data, temp)

data <- data[complete.cases(data),]


# meta-model formula
formula <- runoff ~ precipitation + dune + elev + 
                    p1 + p2 + p3 + pc7 + pc30 + pc180 + days_since_prec + 
                    Tavg + t1 + t2 + tc7 


# train OLS meta-model

set.seed(6854641)
train <- sample(1:nrow(data), nrow(data)*0.75) # random 75% obs for training

linmod <- lm(formula, data = data[train,])
summary(linmod)

# train random forest meta-model
library(ranger)
rf <- ranger(formula, data = data[train,], num.trees =  500, importance = "permutation")
rf


# ------------------------------------------------------------------------------
# Goodness-of-fit metrics (for table 1)


# RF test set
testi <- predict(rf, data = data[-train,])
round(hydroGOF::gof(testi$predictions, data$runoff[-train]),2)
RNP(testi$predictions, data$runoff[-train])
sd(testi$predictions) / sd(data$runoff[-train])

# OLS test set
testi <- predict(linmod, newdata = data[-train,])
round(hydroGOF::gof(testi, data$runoff[-train]),2)
RNP(unname(testi), data$runoff[-train])
sd(testi) / sd(data$runoff[-train])

# RF train set
testi <- rf$predictions
round(hydroGOF::gof(testi, data$runoff[train]),2)
RNP(unname(testi), data$runoff[train])
sd(testi) / sd(data$runoff[train])

# OLS train set
testi <- linmod$fitted.values
round(hydroGOF::gof(testi, data$runoff[train]),2)
RNP(unname(testi), data$runoff[train])
sd(testi) / sd(data$runoff[train])



# ------------------------------------------------------------------------------
# apply meta-models on target zones

target$variable_ts <- lapply(seq_along(target$variable_ts), function(i) {
    x <- target[i,]
    data2 <- create_model_frame(x, formula = runoff ~.)
    temp <- cbind(L(data2$precipitation, 1:3), 
                  zoo::rollsum(data2$precipitation, 7, align = "right", na.pad = TRUE),
                  zoo::rollsum(data2$precipitation, 30, align = "right", na.pad = TRUE),
                  zoo::rollsum(data2$precipitation, 60, align = "right", na.pad = TRUE),
                  zoo::rollsum(data2$precipitation, 180, align = "right", na.pad = TRUE))
    colnames(temp) <- c("p1","p2","p3", "pc7", "pc30", "pc60", "pc180")
    x2 <- cbind(data2, temp)
    temp <- cbind(L(x2$Tavg, 1:2),
                  zoo::rollmean(data2$Tavg, 7, align = "right", na.pad = TRUE))
    colnames(temp) <- c("t1","t2", "tc7")
    x2 <- cbind(x2, temp)
    x2 <- x2[complete.cases(x2),]
    x2 <- rename(x2, fracgrassland.y = fracgrassland)
    
    pred <- predict(rf, data = x2)
    pred <- c(rep(NA, 179), pred$predictions)
    
    pred2 <- predict(linmod, newdata = x2)
    pred2 <- c(rep(NA, 179), pred2)
    
    
    y <- target$variable_ts[[i]]
    y$dm_rf <- pred
    y$dm_lm <- pred2
    
    return(y)
    
})

# save if needed
# saveRDS(target, "target_zones.RDS")


# ------------------------------------------------------------------------------
# evaluate performance of meta-models at target zone level (for table 1)

# RF
data <- lapply(target$variable_ts, function(x) {
    return(dplyr::select(x, runoff, dm_rf))
}) %>% do.call(rbind, .) %>% as_tibble()

hydroGOF::gof(data$dm_rf, drop_units(data$runoff))
RNP(data$dm_rf, drop_units(data$runoff))
sd(data$dm_rf, na.rm=TRUE) / sd(drop_units(data$runoff), na.rm=TRUE)



# OLS
data <- lapply(target$variable_ts, function(x) {
    return(dplyr::select(x, runoff, dm_lm))
}) %>% do.call(rbind, .) %>% as_tibble()

round(hydroGOF::gof(data$dm_lm, 
                    drop_units(data$runoff)),2)
RNP(data$dm_lm, drop_units(data$runoff))
sd(data$dm_lm, na.rm=TRUE) / sd(drop_units(data$runoff), na.rm=TRUE)







# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Interpolate (Dasymetric mapping)


# due to https://github.com/r-spatial/sf/issues/1762, spherical computations off
# this has small, but negligible influence on the computations.
sf_use_s2(FALSE) 


# load data and prepare
is.formula <- function(x){
    inherits(x,"formula")
}


# shift negative meta-model timeseries (affects LM)
target$variable_ts <- lapply(target$variable_ts, function(x) {
    
    x$lm_surrogate <- x$dm_lm
    dm_lm_shifted <- x$dm_lm
    dm_lm_trunc <- x$dm_lm
    
    test <- any(dm_lm_shifted < 0, na.rm=TRUE)
    if(test) {
        ind <- dm_lm_shifted < 0
        dm_lm_shifted <- dm_lm_shifted + abs(min(dm_lm_shifted,na.rm=TRUE))
        dm_lm_trunc[ind] <- 0 
    }  
    
    x$dm_lm_shifted <- dm_lm_shifted
    x$dm_lm_trunc <- dm_lm_trunc
    return(x)
})



# ------------------------------------------------------------------------------
# interpolation using {dasymetric}



ds <- list()

#### AREA WEIGHTING (using source zone values directly)
system.time({
    ds[[ "HS_aw" ]] <- interpolate_aw(source,
                                      target,
                                      variable = "runoff")
}) # ~9 sec


#### Dasymetric mapping with RF meta-model
system.time({
    ds[[ "HS_dm_rf" ]] <- interpolate_dm(source,
                                         target,
                                         variable = "runoff",
                                         dasymetric = "dm_rf",
                                         verbose = TRUE)
}) # ~12.5 sec


#### Dasymetric mapping with OLS meta-model
system.time({
    ds[[ "HS_dm_lm" ]] <- interpolate_dm(source,
                                         target,
                                         variable = "runoff",
                                         dasymetric = "dm_lm_shifted",
                                         verbose = TRUE)
}) # ~12.5 sec






# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Evaluate downscaled performance


## prepare for computing spatial autocorrelation of errors (Moran's I)
dmat <- 1/st_distance(st_centroid(target))^2
diag(dmat) <- 0
listw <- mat2listw(drop_units(dmat), style = "M")

# compute global goodness-of-fit
# global measures
errors <- lapply(seq_along(ds), function(i) {
    x <- ds[[i]]
    name <- names(ds)[i]
    
    temp <- lapply(1:nrow(x), function(ii) {
        # runoff = CWatM simulation
        trgt <- pull(x$variable_ts[[ii]], "runoff") 
        
        # runoff_int = interpolation output
        int <- set_units(pull(x$variable_ts[[ii]],"runoff_int"), "mm/s") 
        err <-  trgt - int
        
        out <- tibble(CWATM = trgt, HS = int, error = err)
    }) %>% do.call(rbind, .)
    
    # goodness-of-fit
    gof <- hydroGOF::gof(units::drop_units(temp$HS), 
                         units::drop_units(temp$CWATM))
    sdev <- sd(temp$HS, na.rm=TRUE) / sd(temp$CWATM, na.rm=TRUE)
    rnp <- RNP(units::drop_units(temp$HS), units::drop_units(temp$CWATM))
    n_neg <- sum(drop_units(temp$HS) < 0, na.rm=TRUE) / sum(!is.na(temp$HS))
    
    #morans i
    mean_error <- sapply(1:nrow(x), function(ii) {
        trgt <- pull(x$variable_ts[[ii]], "runoff")
        
        int <- set_units(pull(x$variable_ts[[ii]],"runoff_int"), "mm/s")
        err <-  trgt - int
        
        out <- mean(err, na.rm=TRUE)
    })
    
    
    moran <- moran.test(mean_error, listw)
    
    gof <- rbind(gof, sdev, matrix(rnp, nrow=4), n_neg, moran$estimate[1])
    
    rownames(gof)[21:26] <- c("sdev", names(rnp), "neg", "Morans_I")
    return(gof)
}) %>% do.call(cbind,.)
colnames(errors) <- names(ds)
errors <- data.frame(errors)

## add meta-model performance to errors-object
data <- lapply(target$variable_ts, function(x) {
    return(dplyr::select(x, runoff, dm_rf))
}) %>% do.call(rbind, .) %>% as_tibble()
errors$RF <- c(hydroGOF::gof(data$dm_rf, drop_units(data$runoff))[,1],
               sd(data$dm_rf, na.rm=TRUE) / sd(drop_units(data$runoff), na.rm=TRUE),
               RNP(data$dm_rf, drop_units(data$runoff)),
               sum(data$dm_rf < 0, na.rm=TRUE) / sum(!is.na(data$dm_rf)),
               NA)
data <- lapply(target$variable_ts, function(x) {
    return(dplyr::select(x, runoff, dm_lm))
}) %>% do.call(rbind, .) %>% as_tibble()
errors$LM <- c(hydroGOF::gof(data$dm_lm, drop_units(data$runoff))[,1],
               sd(data$dm_lm, na.rm=TRUE) / sd(drop_units(data$runoff), na.rm=TRUE),
               RNP(data$dm_lm, drop_units(data$runoff)),
               sum(data$dm_lm < 0, na.rm=TRUE) / sum(!is.na(data$dm_lm)),
               NA)


# Compute Moran's I for meta-models
data <- lapply(target$variable_ts, function(x) {
    tmp <- dplyr::select(x, runoff, dm_rf)
    return(mean(tmp$dm_rf - drop_units(tmp$runoff), na.rm=TRUE))
}) %>% do.call(rbind, .) %>% as_tibble()
errors$RF[27] <- moran.test(data$V1, listw)$estimate[1]


data <- lapply(target$variable_ts, function(x) {
    tmp <- dplyr::select(x, runoff, dm_lm)
    return(mean(tmp$dm_lm - drop_units(tmp$runoff), na.rm=TRUE))
}) %>% do.call(rbind, .) %>% as_tibble()
errors$LM[27] <- moran.test(data$V1, listw)$estimate[1]


data <- lapply(target$variable_ts, function(x) {
    tmp <- dplyr::select(x, runoff, dm_lm_shifted)
    return(mean(tmp$dm_lm_shifted - drop_units(tmp$runoff), na.rm=TRUE))
}) %>% do.call(rbind, .) %>% as_tibble()
errors$LM_shifted[27] <- moran.test(data$V1, listw)$estimate[1]



## COMPUTE KGE SKILL SCORE of downscaled data 
# 0.92 is the theoretical max KGE here because of 8% bias between CWatM 
# simulations at 30- and 5-arcmin
ref <- errors$HS_aw[25]
c(RF = (errors$RF[25] - ref) / (0.92-ref), 
  DM_RF = (errors$HS_dm_rf[25] - ref) / (0.92-ref),
  LM = (errors$LM[25] - ref) / (0.92-ref),
  DM_LM = (errors$HS_dm_lm[25] - ref) / (0.92-ref))



# ------------------------------------------------------------------------------
# Compute goodness-of-fit for 3x3 moving window

# get neighbourhood
touching <- st_touches(target)
ngb <- lapply(seq_along(touching), function(i) {
    x <- touching[[i]]
    ngbs <- lapply(x, function(xx) {
        return(touching[[xx]])
    })
    ngbs <- unique(c(i, unlist(ngbs)))
    return(ngbs)
})

# local gofs
local_errors <- lapply(seq_along(ds), function(i) {
    x <- ds[[i]]
    name <- names(ds)[i]
    
    temp <- lapply(1:nrow(x), function(ii) {
        ngbs <- ngb[[ii]]
        
        tmp <- lapply(ngbs, function(iii) {
            trgt <- units::drop_units(x$variable_ts[[iii]]$runoff)
            int <- units::drop_units(x$variable_ts[[iii]]$runoff_int)
            err <-  trgt - int
            
            out <- tibble(CWATM = trgt, HS = int, error = err)
        }) %>% do.call(rbind, .)
        
        gof <- hydroGOF::gof(tmp$HS, 
                             tmp$CWATM)
        gof <- t(gof[c(5,6,16,19),1])
        # rownames(gof) <- x$targetID[ii]
        
        sdev <- sd(tmp$HS, na.rm=TRUE) / sd(tmp$CWATM, na.rm=TRUE)
        rnp <- RNP(tmp$HS, tmp$CWATM)
        gof <- c(gof, sdev, matrix(rnp, nrow=4))
        names(gof) <- c("NRMSE %", "PBIAS %", "r", "KGE", "sdev", names(rnp))    
        return(gof)
    }) %>% do.call(rbind, .)
    
    return(temp)
})
names(local_errors) <- names(ds)


## add meta-models to local_errors

# RF
x <- target
name <- "RF"

temp <- lapply(1:nrow(x), function(ii) {
    ngbs <- ngb[[ii]]
    
    tmp <- lapply(ngbs, function(iii) {
        trgt <- units::drop_units(x$variable_ts[[iii]]$runoff)
        int <- x$variable_ts[[iii]]$dm_rf
        err <-  trgt - int
        
        out <- tibble(CWATM = trgt, HS = int, error = err)
    }) %>% do.call(rbind, .)
    
    gof <- hydroGOF::gof(tmp$HS, 
                         tmp$CWATM)
    gof <- t(gof[c(5,6,16,19),1])
    # rownames(gof) <- x$targetID[ii]
    
    sdev <- sd(tmp$HS, na.rm=TRUE) / sd(tmp$CWATM, na.rm=TRUE)
    rnp <- RNP(tmp$HS, tmp$CWATM)
    gof <- c(gof, sdev, matrix(rnp, nrow=4))
    names(gof) <- c("NRMSE %", "PBIAS %", "r", "KGE", "sdev", names(rnp))  
    
    return(gof)
}) %>% do.call(rbind, .)
local_errors[["RF"]] <- temp


# OLS
x <- target
name <- "LM"

temp <- lapply(1:nrow(x), function(ii) {
    ngbs <- ngb[[ii]]
    
    tmp <- lapply(ngbs, function(iii) {
        trgt <- units::drop_units(x$variable_ts[[iii]]$runoff)
        int <- x$variable_ts[[iii]]$dm_lm
        err <-  trgt - int
        
        out <- tibble(CWATM = trgt, HS = int, error = err)
    }) %>% do.call(rbind, .)
    
    gof <- hydroGOF::gof(tmp$HS, 
                         tmp$CWATM)
    gof <- t(gof[c(5,6,16,19),1])
    # rownames(gof) <- x$targetID[ii]
    
    sdev <- sd(tmp$HS, na.rm=TRUE) / sd(tmp$CWATM, na.rm=TRUE)
    rnp <- RNP(tmp$HS, tmp$CWATM)
    gof <- c(gof, sdev, matrix(rnp, nrow=4))
    names(gof) <- c("NRMSE %", "PBIAS %", "r", "KGE", "sdev", names(rnp))  
    
    return(gof)
}) %>% do.call(rbind, .)
local_errors[["LM"]] <- temp






# ------------------------------------------------------------------------------
# PLOT LOCAL BIAS between source and target zone simulations (figure 4) 
target$local_bias <- local_errors[["HS_aw"]][,"PBIAS %"]

ggplot() +
    geom_sf(data = target, aes(fill = local_bias),
            color = 'black', size = 0.2) +
    geom_sf(data = source, fill = "transparent", size = 0.5, color = 'black') +
    scale_fill_gradient2(name = "Bias %") +
    theme_minimal()

#ggsave("bias_fig.pdf", width = 8, height = 6)



# ------------------------------------------------------------------------------
# PLOT LOCAL NON-PARAMETRIC KGE (RNP) (figure 5)

tmp <- target[,"targetID"]
for(i in seq_along(local_errors)) {
    name <- names(local_errors)[i]
    tmp[[name]] <- local_errors[[i]][,"RNP"]
}

plotdata <- tmp %>% 
    st_as_sf() %>% 
    select(-targetID) %>% 
    tidyr::gather(ds, value, -geometry) %>% 
    mutate(ds = factor(ds, levels = c("HS_aw", "LM", "HS_dm_lm",
                                      "", "RF", "HS_dm_rf")))

ggplot(plotdata) + 
    geom_sf(aes(fill = (value)), color = 'transparent') + 
    facet_wrap(~ds, ncol = 3, drop = FALSE) +
    scale_fill_gradient2(low = 'blue', mid = "white", high = 'red', 
                         midpoint = -0.41, limits = c(-2,1),
                         oob = scales::squish,
                         name = "RNP") +
    theme_minimal()

# ggsave(filename = "kge_np_map.pdf")








# ------------------------------------------------------------------------------
# Compute spatial correlations

# ------------------------------------------------------------------------------
# spatial correlations

distmat <- st_centroid(target) %>% 
    st_distance() %>% 
    drop_units()
# cormat_cwatm <- distmat
# cormat_cwatm[,] <- NA
# cormat_rf <- cormat_cwatm
# cormat_dmrf <- cormat_cwatm

tbl <- tibble(targetID = rep(target$targetID, each = nrow(target)),
              dist = as.vector(distmat),
              cwatm_5min = NA,
              aw = NA,
              lm = NA,
              rf = NA,
              dm_lm = NA,
              dm_rf = NA)


pb <- txtProgressBar(0, nrow(distmat), style = 3)
for(i in 1:nrow(target)) {
    for(j in i:nrow(target)) {
        if(i == j) next
        
        ind <- (nrow(target) * (i-1)) + j
        
        r1 <- target$variable_ts[[i]]$runoff
        r2 <- target$variable_ts[[j]]$runoff
        pearson <- cor(r1, r2, use = "complete.obs")
        tbl$cwatm_5min[ind] <- pearson
        
        r1 <- target$variable_ts[[i]]$dm_rf
        r2 <- target$variable_ts[[j]]$dm_rf
        pearson <- cor(r1, r2, use = "complete.obs")
        tbl$rf[ind] <- pearson
        
        r1 <- target$variable_ts[[i]]$dm_lm
        r2 <- target$variable_ts[[j]]$dm_lm
        pearson <- cor(r1, r2, use = "complete.obs")
        tbl$lm[ind] <- pearson
        
        r1 <- ds$HS_dm_rf$variable_ts[[i]]$runoff_int
        r2 <- ds$HS_dm_rf$variable_ts[[j]]$runoff_int
        pearson <- cor(r1, r2, use = "complete.obs")
        tbl$dm_rf[ind] <- pearson
        
        r1 <- ds$HS_dm_lm$variable_ts[[i]]$runoff_int
        r2 <- ds$HS_dm_lm$variable_ts[[j]]$runoff_int
        pearson <- cor(r1, r2, use = "complete.obs")
        tbl$dm_lm[ind] <- pearson
        
        r1 <- ds$HS_aw$variable_ts[[i]]$runoff_int
        r2 <- ds$HS_aw$variable_ts[[j]]$runoff_int
        pearson <- cor(r1, r2, use = "complete.obs")
        tbl$aw[ind] <- pearson
    }
    setTxtProgressBar(pb, i)
}
close(pb)
tbl <- tbl[complete.cases(tbl),]


# ------------------------------------------------------------------------------
# PLOT spatial correlations (Figure 6)

temp <- tbl %>% 
    tidyr::gather(model, pearson, -dist, -targetID) %>% 
    mutate(binned_dist = floor(dist/20000)*20,
           # mutate(binned_dist = floor(dist/2000),
           model = factor(model, levels = c("cwatm_5min",
                                            "aw",
                                            "lm", 
                                            "dm_lm",
                                            "rf",
                                            "dm_rf"))) %>% 
    # filter(binned_dist < 50) %>% 
    group_by(binned_dist, model) %>% 
    summarise(quant = quantile(pearson, probs = c(0,0.05,0.25,0.5,0.75,0.95,1)),
              quantile = c(0,0.05,0.25,0.5,0.75,0.95,1)) %>% 
    tidyr::pivot_wider(names_from = quantile, values_from = quant)

colors <-  inlmisc::GetColors(8)[3:7]

p1 <- ggplot() +
    # CWATM ribbon
    # geom_ribbon(data = filter(temp, model == "cwatm_5min"), 
    #             aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "CWATM_5min"), 
    #             alpha = 0.1) +
    # geom_ribbon(data = filter(temp, model == "cwatm_5min"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "CWATM_5min"), 
    #             alpha = 0.2) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.05`, color = "cwatm_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.95`, color = "cwatm_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.5`, color = "cwatm_5min", linetype = "median")) +
    # geom_line(data = filter(temp, model == "cwatm_5min"),
    #           aes(binned_dist, `0.25`, color = "cwatm_5min", linetype = "50%")) +
    # geom_line(data = filter(temp, model == "cwatm_5min"),
    #           aes(binned_dist, `0.75`, color = "cwatm_5min", linetype = "50%")) +
    # AW ribbon
    geom_ribbon(data = filter(temp, model == "aw"), 
                aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "90%"), 
                alpha = 0.2) +
    geom_line(data = filter(temp, model == "aw"), 
              aes(binned_dist, `0.5`), colour = colors[1]) +
    # geom_ribbon(data = filter(temp, model == "aw"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "50%"), 
    #             alpha = 0.3) + 
    scale_fill_manual(name = "AW", values = c(colors[1], colors[1])) + 
    scale_color_manual(values = c("black")) +
    theme_minimal() +
    labs(y = "Pearson correlation coefficient", x = "Distance [km]")
p1

p2 <- ggplot() +
    # CWATM ribbon
    # geom_ribbon(data = filter(temp, model == "cwatm_5min"), 
    #             aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "CWATM_5min"), 
    #             alpha = 0.1) +
    # geom_ribbon(data = filter(temp, model == "cwatm_5min"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "CWATM_5min"), 
    #             alpha = 0.2) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.05`, color = "cwatm_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.95`, color = "cwatm_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.5`, color = "cwatm_5min", linetype = "median")) +
    # geom_line(data = filter(temp, model == "cwatm_5min"),
    #           aes(binned_dist, `0.25`, color = "cwatm_5min", linetype = "50%")) +
    # geom_line(data = filter(temp, model == "cwatm_5min"),
    #           aes(binned_dist, `0.75`, color = "cwatm_5min", linetype = "50%")) +
    # AW ribbon
    geom_ribbon(data = filter(temp, model == "lm"), 
                aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "lm_90%"), 
                alpha = 0.2) +
    geom_line(data = filter(temp, model == "lm"), 
              aes(binned_dist, `0.5`), colour = colors[2]) +
    # geom_ribbon(data = filter(temp, model == "lm"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "lm_50%"), 
    #             alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "dm_lm"), 
                aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%"), 
                alpha = 0.2) +
    geom_line(data = filter(temp, model == "dm_lm"), 
              aes(binned_dist, `0.5`), colour = colors[3]) +
    # geom_ribbon(data = filter(temp, model == "dm_lm"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%"), 
    #             alpha = 0.3) + 
    # scale_fill_manual(name = "LM", values = c(colors[3],colors[3],colors[2],colors[2])) +
    scale_fill_manual(name = "LM", values = c(colors[3],colors[2])) +
    scale_color_manual(values = c("black")) +
    theme_minimal() +
    labs(y = "Pearson correlation coefficient", x = "Distance [km]")
p2

p3 <- ggplot() +
    # CWATM ribbon
    # geom_ribbon(data = filter(temp, model == "cwatm_5min"), 
    #             aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "CWATM_5min"), 
    #             alpha = 0.1) +
    # geom_ribbon(data = filter(temp, model == "cwatm_5min"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "CWATM_5min"), 
    #             alpha = 0.2) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.05`, color = "cwatm_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.95`, color = "cwatm_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "cwatm_5min"),
              aes(binned_dist, `0.5`, color = "cwatm_5min", linetype = "median")) +
    # geom_line(data = filter(temp, model == "cwatm_5min"),
    #           aes(binned_dist, `0.25`, color = "cwatm_5min", linetype = "50%")) +
    # geom_line(data = filter(temp, model == "cwatm_5min"),
    #           aes(binned_dist, `0.75`, color = "cwatm_5min", linetype = "50%")) +
    # AW ribbon
    geom_ribbon(data = filter(temp, model == "rf"), 
                aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "rf_90%"), 
                alpha = 0.2) +
    geom_line(data = filter(temp, model == "rf"), 
              aes(binned_dist, `0.5`), colour = colors[4]) +
    # geom_ribbon(data = filter(temp, model == "rf"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "rf_50%"), 
    #             alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "dm_rf"), 
                aes(binned_dist, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%"), 
                alpha = 0.2) +
    geom_line(data = filter(temp, model == "dm_rf"), 
              aes(binned_dist, `0.5`), colour = colors[5]) +
    # geom_ribbon(data = filter(temp, model == "dm_rf"), 
    #             aes(binned_dist, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%"), 
    #             alpha = 0.3) + 
    # scale_fill_manual(name = "RF", values = c(colors[5],colors[5],colors[4],colors[4])) 
    scale_fill_manual(name = "RF", values = c(colors[5],colors[4])) +
    scale_color_manual(values = c("black")) +
    theme_minimal() +
    labs(y = "Pearson correlation coefficient", x = "Distance [km]")
p3

p1 / p2 / p3 + plot_layout(guides = "collect")

# ggsave(filename = "spatial_correlation2.pdf", width = 6, height = 12)





### ---------------------------------------------------------------------------
# plot timeseries (long code)   FIGURE 7


ind <- st_within(target, source) %>% unlist
target$sourceID <- ind

data <- lapply(unique(ind), function(i) {
    
    sourceID <- source$sourceID[[i]]
    source_index <- i
    
    t <- filter(target, sourceID == i)
    out <- lapply(1:nrow(t), function(ii) {
        temp <- dplyr::select(t$variable_ts[[ii]], 
                              Date, 
                              CWATM_5min = runoff, 
                              LM = dm_lm,
                              RF = dm_rf) %>% 
            mutate(targetID = t$targetID[ii],
                   CWATM_5min = drop_units(CWATM_5min),
                   .after = Date)
    }) %>% do.call("rbind", .)
    
    t <- ds$HS_aw %>% 
        mutate(sourceID = ind) %>%  
        filter(sourceID == i)
    temp <- lapply(1:nrow(t), function(ii) {
        temp <- dplyr::select(t$variable_ts[[ii]], AW = runoff_int) %>% 
            pull() %>% 
            drop_units()
    }) %>% do.call("c", .)
    out$AW <- temp
    
    t <- ds$HS_dm_lm %>% 
        mutate(sourceID = ind) %>%  
        filter(sourceID == i)
    temp <- lapply(1:nrow(t), function(ii) {
        temp <- dplyr::select(t$variable_ts[[ii]],DM_LM = runoff_int) %>% 
            pull() %>% 
            drop_units()
    }) %>% do.call("c", .)
    out$DM_LM <- temp
    
    t <- ds$HS_dm_rf %>% 
        mutate(sourceID = ind) %>%  
        filter(sourceID == i)
    temp <- lapply(1:nrow(t), function(ii) {
        temp <- dplyr::select(t$variable_ts[[ii]], DM_RF = runoff_int) %>% 
            pull() %>% 
            drop_units()
    }) %>% do.call("c", .)
    out$DM_RF <- temp
    
    out <- mutate(out,
                  sourceID = sourceID,
                  source_index = source_index,
                  .after = targetID)
    
    return(out)
    
})


# each plot showing different period, LONG CODE
i <- 12
datestart <- as.Date("2006-09-01")
dateend <- as.Date("2006-10-15")
plotdata <- data[[i]] %>% 
    filter(Date < dateend,
           Date > datestart) %>% 
    select(-targetID, -sourceID, -source_index) %>% 
    group_by(Date) %>% 
    summarise_all(.funs = quantile, probs = c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE) %>% 
    mutate(quant = c(0.05,0.25,0.5,0.75,0.95),
           facet = case_when(quant %in% c(0.05, 0.95) ~ "90%",
                             quant %in% c(0.25, 0.75) ~ "50%"))

temp <- plotdata %>% 
    tidyr::gather(model, value, -Date, -quant, -facet) %>% 
    tidyr::pivot_wider(names_from = quant, values_from = value) %>% 
    filter(!is.na(facet),
           facet == "90%")

p1 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # AW ribbon
    geom_ribbon(data = filter(temp, model == "AW"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "lm_90%", color = "lm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "AW"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "lm_50%", color = "lm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "AW", values = c(colors[1],colors[1])) + 
    scale_color_manual(values = c("black", colors[1], colors[1])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p2 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # LM ribbon
    geom_ribbon(data = filter(temp, model == "LM"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "lm_90%", color = "lm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "LM"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "lm_50%", color = "lm_50%"), 
                alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "DM_LM"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%", color = "dm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "DM_LM"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%", color = "dm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "LM", values = c(colors[3],colors[3],colors[2],colors[2])) + 
    scale_color_manual(values = c("black", colors[3],colors[3],colors[2],colors[2])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p3 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # RF ribbon
    geom_ribbon(data = filter(temp, model == "RF"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "rf_90%", color = "rf_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "RF"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "rf_50%", color = "rf_50%"), 
                alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "DM_RF"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%", color = "dm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "DM_RF"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%", color = "dm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "RF", values = c(colors[5],colors[5],colors[4],colors[4])) + 
    scale_color_manual(values = c("black", colors[5],colors[5],colors[4],colors[4])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p1 / p2 / p3 + plot_layout(guides = "collect")


## second period plots
i <- 2
datestart <- as.Date("2004-06-15")
dateend <- as.Date("2004-07-31")
plotdata <- data[[i]] %>% 
    filter(Date < dateend,
           Date > datestart) %>% 
    select(-targetID, -sourceID, -source_index) %>% 
    group_by(Date) %>% 
    summarise_all(.funs = quantile, probs = c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE) %>% 
    mutate(quant = c(0.05,0.25,0.5,0.75,0.95),
           facet = case_when(quant %in% c(0.05, 0.95) ~ "90%",
                             quant %in% c(0.25, 0.75) ~ "50%"))

temp <- plotdata %>% 
    tidyr::gather(model, value, -Date, -quant, -facet) %>% 
    tidyr::pivot_wider(names_from = quant, values_from = value) %>% 
    filter(!is.na(facet),
           facet == "90%")

p4 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # AW ribbon
    geom_ribbon(data = filter(temp, model == "AW"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "lm_90%", color = "lm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "AW"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "lm_50%", color = "lm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "AW", values = c(colors[1],colors[1])) + 
    scale_color_manual(values = c("black", colors[1], colors[1])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p5 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # LM ribbon
    geom_ribbon(data = filter(temp, model == "LM"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "lm_90%", color = "lm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "LM"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "lm_50%", color = "lm_50%"), 
                alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "DM_LM"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%", color = "dm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "DM_LM"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%", color = "dm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "LM", values = c(colors[3],colors[3],colors[2],colors[2])) + 
    scale_color_manual(values = c("black", colors[3],colors[3],colors[2],colors[2])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p6 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # RF ribbon
    geom_ribbon(data = filter(temp, model == "RF"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "rf_90%", color = "rf_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "RF"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "rf_50%", color = "rf_50%"), 
                alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "DM_RF"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%", color = "dm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "DM_RF"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%", color = "dm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "RF", values = c(colors[5],colors[5],colors[4],colors[4])) + 
    scale_color_manual(values = c("black", colors[5],colors[5],colors[4],colors[4])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p4 / p5 / p6 + plot_layout(guides = "collect")

p1 + p4 + p2 + p5 + p3 + p6 + plot_layout(guides = "collect", nrow = 3)


# ------------------------------------------------------------------------------
## third panel
i <- 20
datestart <- as.Date("2003-06-15")
dateend <- as.Date("2005-11-30")
plotdata <- data[[i]] %>% 
    filter(Date < dateend,
           Date > datestart) %>% 
    select(-targetID, -sourceID, -source_index) %>% 
    group_by(Date) %>% 
    summarise_all(.funs = quantile, probs = c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE) %>% 
    mutate(quant = c(0.05,0.25,0.5,0.75,0.95),
           facet = case_when(quant %in% c(0.05, 0.95) ~ "90%",
                             quant %in% c(0.25, 0.75) ~ "50%"))

temp <- plotdata %>% 
    tidyr::gather(model, value, -Date, -quant, -facet) %>% 
    tidyr::pivot_wider(names_from = quant, values_from = value) %>% 
    filter(!is.na(facet),
           facet == "90%")

p7 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # AW ribbon
    geom_ribbon(data = filter(temp, model == "AW"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "lm_90%", color = "lm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "AW"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "lm_50%", color = "lm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "AW", values = c(colors[1],colors[1])) + 
    scale_color_manual(values = c("black", colors[1], colors[1])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p8 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # LM ribbon
    geom_ribbon(data = filter(temp, model == "LM"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "lm_90%", color = "lm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "LM"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "lm_50%", color = "lm_50%"), 
                alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "DM_LM"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%", color = "dm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "DM_LM"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%", color = "dm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "LM", values = c(colors[3],colors[3],colors[2],colors[2])) + 
    scale_color_manual(values = c("black", colors[3],colors[3],colors[2],colors[2])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p9 <- ggplot() +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.05`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.95`, color = "CWATM_5min", linetype = "90%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.25`, color = "CWATM_5min", linetype = "50%")) +
    geom_line(data = filter(temp, model == "CWATM_5min"),
              aes(Date, `0.75`, color = "CWATM_5min", linetype = "50%")) +
    # RF ribbon
    geom_ribbon(data = filter(temp, model == "RF"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "rf_90%", color = "rf_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "RF"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "rf_50%", color = "rf_50%"), 
                alpha = 0.3) + 
    geom_ribbon(data = filter(temp, model == "DM_RF"), 
                aes(Date, ymin=`0.05`, ymax=`0.95`, fill = "dm_90%", color = "dm_90%"), 
                alpha = 0.2) +
    geom_ribbon(data = filter(temp, model == "DM_RF"), 
                aes(Date, ymin=`0.25`, ymax=`0.75`, fill = "dm_50%", color = "dm_50%"), 
                alpha = 0.3) + 
    scale_fill_manual(name = "RF", values = c(colors[5],colors[5],colors[4],colors[4])) + 
    scale_color_manual(values = c("black", colors[5],colors[5],colors[4],colors[4])) +
    theme_minimal() +
    labs(y = "Runoff [mm/s]", x = "") + 
    facet_wrap(~facet, ncol = 2)

p7 / p8 / p9 + plot_layout(guides = "collect")

p1 + p4 + p7 + p2 + p5 + p8 + p3 + p6 + p9 + plot_layout(guides = "collect", nrow = 3)

# ggsave(filename = "timeseries.pdf", width = 16, height = 12)
