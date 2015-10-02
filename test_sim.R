source("sim_func.R")

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
                # Correlations:
                manual = 0,
                # If 'manual' == T, then specify correlations
                ## Among-patients:
                raG1G2 = 0,
                raP1P2 = 0, 
                raP1G1 = 0,
                raP2G2 = 0,
                raP1G2 = 0,
                raP2G1 = 0,
                ## Within-patients:
                rwG1G2 = 0,
                rwP1P2 = 0, 
                rwP1G1 = 0,
                rwP2G2 = 0,
                rwP1G2 = 0,
                rwP2G1 = 0,
                ## sd across patients for each strain (phi and gamma):
                ## Assume all sd equal
                sa = 1,
                ## sd within patients for each strain (phi and gamma):
                ## Assume all sd equal
                sw = .4,
                # If 'manual' == FALSE, generate pos. definite covariance matrix automatically
                ## eta controls the degree of correlation:
                etaA = 2, # Among-patient eta
                etaW = 1.5, # Within-patient eta
                #Global probabilities:
                globphi = .5, 
                globgam = .5, 
                globpsi = .5
)

library(ggplot2)
n1 <- subset(Occ, Strain == 1)[,2:4]
n2 <- subset(Occ, Strain == 2)[,2:4]

df <- cbind(n1,n2)
colnames(df) <- c("psi1", "phi1", "gam1", "psi2", "phi2", "gam2")
df <- data.frame(df)


p1 <- qplot(x=phi1, y=phi2, data=df, geom="point", xlim=c(0,1), ylim=c(0,1),
            xlab = expression(paste(phi, " Strain 1")),
            ylab = expression(paste(phi, " Strain 2")))

p2 <- qplot(x=gam1, y=gam2, data=df, geom="point", xlim=c(0,1), ylim=c(0,1),
            xlab = expression(paste(gamma, " Strain 1")),
            ylab = expression(paste(gamma, " Strain 2")))

p3 <- qplot(x=phi1, y=gam2, data=df, geom="point", xlim=c(0,1), ylim=c(0,1),
            xlab = expression(paste(phi, " Strain 1")),
            ylab = expression(paste(gamma, " Strain 2")))

p4 <- qplot(x=phi2, y=gam1, data=df, geom="point", xlim=c(0,1), ylim=c(0,1),
            xlab = expression(paste(phi, " Strain 2")),
            ylab = expression(paste(gamma, " Strain 1")))

p5 <- qplot(x=psi1, y=psi2, data=df, geom="point", xlim=c(0,1), ylim=c(0,1),
            xlab = expression(paste(psi, " Strain 1")),
            ylab = expression(paste(psi, " Strain 2")))

library(gridExtra)
quartz(height=3, width=15)
grid.arrange(p1,p2,p3,p4,p5, nrow=1)

# Within and among correlations...

Among <- sim_func(raP1G2 = 0.9)

a1 <- subset(Among, Strain == 1)
a2 <- subset(Among, Strain == 2)

among <- cbind(a2$gam, a1$phi, a1[,6])
colnames(among) <- c("gam2", "phi1", "Pat")
among <- data.frame(among)
among <- subset(among, Pat < 16)

p1 <- ggplot(among, aes(x=phi1, y=gam2, color=factor(Pat)))+
  geom_point()+
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,1))+
  labs(x = expression(paste(phi, " Strain 1")), y = expression(paste(gamma, " Strain 2")))+
  ggtitle("Only Among-patient Corr")

Both <- sim_func(raP1G2 = 0.9, rwP1G2 = -0.9)
b1 <- subset(Both, Strain == 1)
b2 <- subset(Both, Strain == 2)

both <- cbind(b2$gam, b1$phi, b1[,6])
colnames(both) <- c("gam2", "phi1", "Pat")
both <- data.frame(both)
both <- subset(both, Pat < 16)

p2 <- ggplot(both, aes(x=phi1, y=gam2, color=factor(Pat)))+
  geom_point()+
  theme(legend.position = "none")+
  scale_y_continuous(limits=c(0,1))+
  scale_x_continuous(limits=c(0,1))+
  labs(x = expression(paste(phi, " Strain 1")), y = expression(paste(gamma, " Strain 2")))+
  ggtitle("Both Corr")

quartz(height=5, width=10)
grid.arrange(p1,p2, ncol=2)  
