
Occ <- sim_func(n.pat = 100,
                n.vis = 10, 
                # Within and among host covariate effects for phi and gamma:
                bpat1g = 0, 
                bpat2g = 0,
                btime1g = 0,
                btime2g = 0,
                bpat1p = 0,
                bpat2p = 0,
                btime1p = 0,
                btime2p = 0,
                rag = 0, #rho.across.gamma
                rap = 0.9, #rho.across.phi
                rwg = 0, #rho.within.gamma
                rwp = 0, #rho.within.phi
                # sd across patients for each strain (phi and gamma):
                sap1 = 1, 
                sap2 = 1,
                sag1 = 1,
                sag2 = 1,
                # sd within patients for each strain (phi and gamma):
                swp1 = .4, 
                swp2 = .4,
                swg1 = .4,
                swg2 = .4,
                #Global probabilities:
                globphi = .5, 
                globgam = .5, 
                globpsi = .5 
)

t1 <- subset(Occ, Strain == 1)
t2 <- subset(Occ, Strain == 2)

plot(t1$lphi ~ t2$lphi)
plot(t1$lgam ~ t2$lgam)
plot(t1$psi ~ t2$psi)

comb <- cbind(t2$psi, t1$psi, t1[,6])
head(comb)
colnames(comb) <- c("psi2", "psi1", "Pat")
comb <- data.frame(comb)

comb <- subset(comb, Pat < 15)
library(ggplot2)
ggplot(comb, aes(x=psi2, y=psi1, color=factor(Pat)))+
  geom_point()+
  theme_classic()

# Weird GLM
m1 <- glm(t1$Y ~ t2$Y, family="binomial")
summary(m1)