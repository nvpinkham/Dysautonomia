source("human_functions.22.03.15.R")


stool.norm <- norm.human(metabolites = meta.stool, samples.pick = row.names(meta.stool))

stool.norm$Label <- NULL
s.d <- dist(stool.norm)


stool.plsda <- mixOmics::splsda(stool.norm, map.stool$Disease.state, ncomp = 3)

sigs <- round(stool.plsda$prop_expl_var$X * 100)
sigs <- paste0("X-variate ", 1:3, ": ", sigs, "%")

map.stool.plsda <- with(map.stool, vegan3d::ordiplot3d(stool.plsda$variates$X,
                                                       pch = 21,
                                                       ax.col = NULL, 
                                                       bg = map.stool$col, 
                                                       xlab = sigs[1], 
                                                       ylab = sigs[2], 
                                                       zlab = sigs[3]))


#setEPS()       # Set postscript arguments
#postscript("Stool_plsda.eps", width = 5, height = 5) 

scatterplot3d(stool.plsda$variates$X,
              pch = 21,
              ax.col = NULL, 
              bg = map.stool$col, 
              xlab = sigs[1], 
              ylab = sigs[2], 
              zlab = sigs[3])

loadings <- map.stool.plsda$points

num.pair <- NULL

for(i in map.stool$Family.ID){
  
  
  loadings.pick <- loadings[map.stool$Family.ID == i,]
  map.pick <- map.stool[map.stool$Family.ID == i,]
  
  if(length(unique(map.pick$Disease.state)) > 1){
    
    
    loadings.pick.sick <- matrix(loadings.pick[map.pick$Disease.state == "Patient",], ncol = 2, byrow = F)
    loadings.pick.control <- matrix(loadings.pick[map.pick$Disease.state == "Relative",],ncol = 2, byrow = F)
    
    for(j in 1 : nrow(loadings.pick.sick)){
      
      x0 <- loadings.pick.sick[j, 1]
      y0 <- loadings.pick.sick[j, 2]
      
      for(k in 1: nrow(loadings.pick.control)){
        
        num.pair <- c(num.pair, 1)
        arrows(x0,
               y0, 
               loadings.pick.control[k,1], 
               loadings.pick.control[k,2], 
               col = 8, angle = 0)
        
        
      }
    }
  }
}

points(loadings,
       pch = 21,
       bg = map.stool$col)

set.seed(42)
res <- adonis2(dist(stool.norm) ~ map.stool$Disease.state)
res <- paste("distance ~ disease state \n f stat", 
             round(res$F[1], 3), "p val", res$`Pr(>F)`[1])
title(res)

table(map.stool$Disease.state)
length(num.pair)
#dev.off()

table(map.serum$Disease.state)

ps <- lme.ps <- NULL
norm1 <- norm2 <- NULL

stool.norm <- stool.norm[!is.na(map.stool$pam) , ]
meta.stool <- meta.stool[!is.na(map.stool$pam) , ]
map.stool <- map.stool[!is.na(map.stool$pam) , ]




for(i in 1 : ncol(stool.norm)){
  
  colnames(stool.norm)
  
  norm1[i] <- shapiro.test(stool.norm[ map.stool$pam == 1, i])$p.value
  norm2[i] <- shapiro.test(stool.norm[ map.stool$pam == 2, i])$p.value

 # ps[i] <- t.test( stool.norm[,i] ~ map.stool$pam)$p.value
  ps[i] <- wilcox.test( stool.norm[,i] ~ map.stool$pam)$p.value
  
  mod <- lmer( stool.norm[,i] ~ map.stool$pam + map.stool$Family.ID + (1 | map.stool$DonorID), REML = F)
  mod.null <- lmer( stool.norm[,i] ~ map.stool$Family.ID + (1 | map.stool$DonorID), REML = F)
  
  res <- anova(  mod,   mod.null)
  lme.ps[i] <- res$`Pr(>Chisq)`[2]
}

ps <- p.adjust(ps)
lme.ps <- p.adjust(  lme.ps)

colnames(stool.norm)[ps < 0.05]
colnames(stool.norm)[lme.ps < 0.05]

boxplot(log(meta.stool$Choline) ~ map.stool$pam, ylab = "Choline", xlab = "Enterotype")

boxplot(log(meta.stool$Ornithine) ~ map.stool$pam, ylab = "Ornithine", xlab = "Enterotype")

boxplot(log(meta.stool$Taurine) ~ map.stool$pam, ylab = "Taurine", xlab = "Enterotype")



colnames(stool.norm)[norm1 < 0.05]
colnames(stool.norm)[norm2 < 0.05]








