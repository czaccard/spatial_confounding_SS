library(tidyverse)
library(geoR)
library(mvtnorm)

outpath = "results/"
inpath = "func/"
plot_path = "plots/"

setting <- readxl::read_excel(paste0(inpath, "setting.xlsx"), sheet = "Foglio2")
source(paste0(inpath, "data_generate.r"), local = T)
source(paste0(inpath, "sim_functions.r"), local = T)
source(paste0(inpath, "misc.r"), local = T)

beta.real = c(1,2) %>% as.matrix(ncol=1) # vector giving intercept coefficient and beta_x
ind = 2 # 2: indicates the exposure coefficient in a design matrix like cbind(1,X)

M = 2^6 # Length of one side of grid 
griddf = expand.grid(lngcoords = seq(0, 1, l=M+1)[-(M+1)],
                     latcoords = seq(0, 1, l=M+1)[-(M+1)])

n0 = 500 # number of sampled locations
n_new = 50 # to make predictions
set.seed(231)
yind <- sample(1:(M^2), size=n0+n_new)

coords = griddf[yind,]
D = dist(coords, diag = T, upper = T) %>% as.matrix() %>% unname()
rownames(D) <- NULL
coords_new = tail(coords, n_new)

n_sim = 100
mu1 = rep(0, n0+n_new)
mu2 = rep(0, n0+n_new)

type_dgm = "RelBias_deltafixed"
options_model = 1


configs = c(7, 3, 5)
out = out_m = data.frame()
for (j in configs) {
  DWA_ = readRDS(paste0(outpath, "DWA/sim", j,"_", type_dgm,".RDS"))
  ols = DWA_$beta_results[,1]
  subfolder = 'SS_fv'
  SS_ = readRDS(paste0(outpath, subfolder, options_model, "/sim", j,"_", type_dgm,".RDS"))
  subfolder = 'SS_nmig'
  SS_nmig = readRDS(paste0(outpath, subfolder, options_model, "/sim", j,"_", type_dgm,".RDS"))
  subfolder = 'SS_mom'
  SS_mom = readRDS(paste0(outpath, subfolder, options_model, "/sim", j,"_", type_dgm,".RDS"))
  
  appo1 = 
    data.frame(replicate = 1:n_sim, SS_fv = SS_$d_pmean[,2], 
               SS_nmig = SS_nmig$d_pmean[,2], SS_mom = SS_mom$d_pmean[,2]) %>% 
    gather(key = "Prior", value = "diff", -replicate) %>% group_by(Prior) %>%
    mutate(analytic = factor(1))
  appo2 = data.frame(replicate = 1:n_sim, SS_fv = SS_$pmeanbeta[,2]-ols, 
                     SS_nmig = SS_nmig$pmeanbeta[,2]-ols, 
                     SS_mom = SS_mom$pmeanbeta[,2]-ols) %>% 
    gather(key = "Prior", value = "diff", -replicate) %>% group_by(Prior) %>%
    mutate(analytic = factor(0))
  appo = bind_rows(appo1, appo2)
  out_m = rbind(out_m, appo)
  
  appo3 =
    data.frame(replicate = 1:n_sim, SS_fv = SS_$d_bias,
               SS_nmig = SS_nmig$d_bias, SS_mom = SS_mom$d_bias) %>% 
    gather(key = "Prior", value = "diff", -replicate) %>% group_by(Prior) %>%
    mutate(analytic = factor(1))
  out = rbind(out, appo3)
}

out$Configuration = c(rep('I',n_sim*3), rep('II',n_sim*3), rep('III',n_sim*3))
out$Configuration = factor(out$Configuration, levels = c('I', 'II', 'III'),
                           labels = c(expression('Configuration I ('*phi[x]* '= 0.05, '*phi[w]*'= 0.5)'), 
                                      expression('Configuration II ('*phi[x]* '= 0.5, '*phi[w]*'= 0.05)'),
                                      expression('Configuration III ('*phi[x]* '= 0.2, '*phi[w]*'= 0.2)')))
out$Difference = factor("italic(d[x])")
out_m$Configuration = rep(out$Configuration, each=2)
out_m$Difference = factor("italic(bar(d[x]^{phantom() ~ symbol('*') ~ phantom()}))")

# THIS IS FIGURE C.1 ####
g3 = ggplot(rbind(out, out_m), aes(x=Prior, y=diff, linetype = analytic)) + 
  geom_hline(aes(yintercept=0), linetype='dotted', color='grey') +
  geom_boxplot(fill = "#D0DDD7") + 
  scale_y_continuous(labels=function(x) sprintf("%.2f", x)) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  facet_grid(Difference ~ Configuration, scales='free', 
             labeller = label_parsed) +
  theme_bw() + xlab("Prior Structure") +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.title.y = element_blank())
g3


ggsave(paste0(plot_path, "differences d and dstar.pdf"), g3,
       width = 7.5, height = 4, device = "pdf")
####



threshold = 0.5
psplines = principal.spline.2D(coords, 1:n0, ~lngcoords+latcoords)
B = psplines$F_n[,-1]

configs = c("I" = 7,"II" = 3,"III" = 5)
out_all = out2_all = data.frame()
for(settingplot in configs) {
  DG = generate_data(setting, settingplot, type_dgm, coords, D, 0.15)
  x = DG$X[1,1:n0]
  corxB = drop(cor(x, B))
  corgB = apply(DG$G, 1, function(g) cor(B, g[1:n0]))
  
  out = out2 = data.frame()
  # out2 is for linear combination -> cor(x or w, lin.comb. of bases)
  for(subfolder in c('SS_fv', 'SS_nmig', 'SS_mom')) {
    SS_ = readRDS(paste0(outpath, subfolder, options_model, "/sim", settingplot,"_", type_dgm,".RDS"))
    cor_g_b = sapply(1:n_sim, function(i) (SS_$pip[i,]>threshold) * corgB[,i])
    cor_g_b[cor_g_b==0] = NA
    cor_x_b = sapply(1:n_sim, function(i) (SS_$pip[i,]>threshold) * corxB)
    cor_x_b[cor_x_b==0] = NA
    appo = data.frame(k = rep(1:ncol(B), n_sim), 
                      Exposure = as.vector(cor_x_b),
                      Confounder = as.vector(cor_g_b))
    out = rbind(out, appo)
    appo = data.frame(irep = 1:n_sim,
                      Exposure = sapply(1:n_sim, function(i) cor(x, B%*%SS_$pmeanbasis[i,])),
                      Confounder = sapply(1:n_sim, function(i) cor(DG$G[i,1:n0], B%*%SS_$pmeanbasis[i,])))
    out2 = rbind(out2, appo)
  }
  out$Prior = rep(c('SS_fv', 'SS_nmig', 'SS_mom'), each=ncol(B)*n_sim)
  out2$Prior = rep(c('SS_fv', 'SS_nmig', 'SS_mom'), each=n_sim)
  out = out %>% 
    pivot_longer(Exposure:Confounder, names_to = "Correlation with", values_to = "corr") %>% 
    mutate(acorr = abs(corr), Configuration = names(configs)[match(settingplot,configs)])
  out2 = out2 %>% 
    pivot_longer(Exposure:Confounder, names_to = "Correlation with", values_to = "corr") %>% 
    mutate(acorr = abs(corr), Configuration = names(configs)[match(settingplot,configs)])
  
  out_all = rbind(out_all, out)
  out2_all = rbind(out2_all, out2)
}

out2_all$Configuration = factor(out2_all$Configuration, levels = c('I', 'II', 'III'),
                           labels = c(expression('Configuration I ('*phi[x]* '= 0.05, '*phi[w]*'= 0.5)'), 
                                      expression('Configuration II ('*phi[x]* '= 0.5, '*phi[w]*'= 0.05)'),
                                      expression('Configuration III ('*phi[x]* '= 0.2, '*phi[w]*'= 0.2)')))

# THIS IS FIGURE C.2 ####
g8 = ggplot(out2_all %>% filter(`Correlation with`=="Confounder"), 
            aes(Prior, acorr)) +
  geom_boxplot(fill = "#D0DDD7") +
  ylab("Correlation (absolute value)") +
  theme_bw() + facet_grid(~Configuration, labeller = label_parsed) +
  theme(panel.grid = element_blank())
g8
ggsave(paste0(plot_path, "correlation g with lin comb.pdf"), g8,
       width = 7.5, height = 4, device = "pdf")




