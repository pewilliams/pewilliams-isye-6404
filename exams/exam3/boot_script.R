library(boot)
library(data.table)
library(np)
B <- 1000 #bootstrap sample count
#estimator function
adat <- readRDS('/Users/peterwilliams/Projects/pewilliams-isye-6404/exams/exam3/data/auto_mpg.rds')
cpt <- data.table(disp = median(adat$disp),cyl = median(adat$cyl), accel = median(adat$accel), type = 'center point')
opt <- data.table(disp = max(adat$disp) - 1,cyl = max(adat$cyl)-1, accel = max(adat$accel) - 1, type = 'outside point')
spl_boot_cpt <- function(data, indices){
  #indices <- sample(1:nrow(adat),replace = T)
  bs_dat <- data[indices]
  sfit_x1 <- smooth.spline(bs_dat$disp,bs_dat$mpg, cv= F,df = 3)
  sfit_x2 <- smooth.spline(bs_dat$cyl,residuals(sfit_x1), cv= F,df = 3)
  sfit_x3 <- smooth.spline(bs_dat$accel,residuals(sfit_x2), cv= F,df = 3)
  prediction <- as.numeric(predict(sfit_x1, data.frame(cpt$disp))$y +
                             predict(sfit_x2, data.frame(cpt$cyl))$y + predict(sfit_x3, data.frame(cpt$accel))$y)
  prediction
}

spl_boot_opt <- function(data, indices){
  #indices <- sample(1:nrow(adat),replace = T)
  bs_dat <- data[indices]
  sfit_x1 <- smooth.spline(bs_dat$disp,bs_dat$mpg, cv= F,df = 3)
  sfit_x2 <- smooth.spline(bs_dat$cyl,residuals(sfit_x1), cv= F,df = 3)
  sfit_x3 <- smooth.spline(bs_dat$accel,residuals(sfit_x2), cv= F,df = 3)
  prediction <- as.numeric(predict(sfit_x1, data.frame(opt$disp))$y +
                             predict(sfit_x2, data.frame(opt$cyl))$y + predict(sfit_x3, data.frame(opt$accel))$y)
  prediction
}

spl_cpt_res <- boot(adat, statistic = spl_boot_cpt,
                    R=B)
spl_opt_res <- boot(adat, statistic = spl_boot_opt,
                    R=B)

cpt_results <- boot.ci(boot.out = spl_cpt_res,type = 'bca',conf = 0.9)
opt_results <- boot.ci(boot.out = spl_opt_res,type = 'bca',conf = 0.9)

boot_spl_bca <- data.table(
  Point = c('Center Point','Outside Point'),
  Prediction = c(cpt_results[[2]],opt_results[[2]]),
  Lower = c(cpt_results$bca[4],opt_results$bca[4]),
  Upper = c(cpt_results$bca[5],opt_results$bca[5]))



kr_boot_opt <- function(data, indices){
  #indices <- sample(1:nrow(adat),replace = T)
  kdat <- data[indices]
  kfit_x1 <- npreg(mpg~disp,gradients = TRUE, data = kdat)
  kdat$e1 <- residuals(kfit_x1)
  kfit_x2 <- npreg(e1~cyl,gradients = TRUE, data = kdat)
  kdat$e2 <- residuals(kfit_x2)
  kfit_x3 <- npreg(e2~accel,gradients = TRUE, data = kdat)
  prediction <- predict(kfit_x1, newdata = data.frame(disp = opt$disp)) + 
    predict(kfit_x2, newdata = data.frame(cyl = opt$cyl)) +
    predict(kfit_x3, newdata = data.frame(accel = opt$accel))
  prediction
}

kr_boot_cpt <- function(data, indices){
  #indices <- sample(1:nrow(adat),replace = T)
  kdat <- data[indices]
  kfit_x1 <- npreg(mpg~disp,gradients = TRUE, data = kdat)
  kdat$e1 <- residuals(kfit_x1)
  kfit_x2 <- npreg(e1~cyl,gradients = TRUE, data = kdat)
  kdat$e2 <- residuals(kfit_x2)
  kfit_x3 <- npreg(e2~accel,gradients = TRUE, data = kdat)
  prediction <- predict(kfit_x1, newdata = data.frame(disp = cpt$disp)) + 
    predict(kfit_x2, newdata = data.frame(cyl = cpt$cyl)) +
    predict(kfit_x3, newdata = data.frame(accel = cpt$accel))
  prediction
}

kr_cpt_res <- boot(adat, statistic = kr_boot_cpt,
                   R=B,parallel = 'multicore',
                   ncpus = 11)
kr_opt_res <- boot(adat, statistic = kr_boot_opt,
                   R=B,parallel = 'multicore',
                   ncpus = 11)

cpt_kresults <- boot.ci(boot.out = kr_cpt_res,type = 'bca',conf = 0.9)
opt_kresults <- boot.ci(boot.out = kr_opt_res,type = 'bca',conf = 0.9)

boot_kr_bca <- data.table(
  Point = c('Center Point','Outside Point'),
  Prediction = c(cpt_kresults[[2]],opt_kresults[[2]]),
  Lower = c(cpt_kresults$bca[4],opt_kresults$bca[4]),
  Upper = c(cpt_kresults$bca[5],opt_kresults$bca[5]))

boot_kr_bca$Method <- 'Kernel Regression'
boot_spl_bca$Method <- 'Spline'

boot_sum <- rbindlist(list(boot_kr_bca,boot_spl_bca))

saveRDS(boot_sum, file = '/Users/peterwilliams/Projects/pewilliams-isye-6404/exams/exam3/data/boot_res.rds')
