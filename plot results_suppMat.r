# FIGURES IN WEB APPENDICES

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
compare_scenarios = tibble()
for (settingplot in configs) {
  DWA_ = readRDS(paste0(outpath, "DWA/sim", settingplot,"_", type_dgm,".RDS"))
  KS_ = readRDS(paste0(outpath, "KS/sim", settingplot,"_", type_dgm,".RDS"))
  subfolder = 'SS_fv'
  SS_ = readRDS(paste0(outpath, subfolder, options_model, "/sim", settingplot,"_", type_dgm,".RDS"))
  subfolder = 'SS_nmig'
  SS_nmig = readRDS(paste0(outpath, subfolder, options_model, "/sim", settingplot,"_", type_dgm,".RDS"))
  subfolder = 'SS_mom'
  SS_mom = readRDS(paste0(outpath, subfolder, options_model, "/sim", settingplot,"_", type_dgm,".RDS"))
  SA_ = readRDS(paste0(outpath, "SpectralAdj/sim", settingplot,"_", type_dgm,".RDS"))
  subfolder = 'SRE'
  SRE_ = readRDS(paste0(outpath, subfolder, "/sim", settingplot,"_", type_dgm,".RDS"))
  
  DWA_$input
  compare_names<-c("OLS","SpatialTP_fx","RSR_fx","gSEM_fx","Spatial+_fx","SpatialTP",
                   "RSR","gSEM","Spatial+", "KS", "SS_fv", "SS_nmig", "SS_mom", "SA",
                   "SRE")
  compare_beta = data.frame(DWA_$beta_results, KS_$beta_results, SS_$pmeanbeta[,ind],
                            SS_nmig$pmeanbeta[,ind], SS_mom$pmeanbeta[,ind],
                            SA_$pmeanbeta[,ind], SRE_$beta_results)
  colnames(compare_beta) = compare_names
  compare_beta = compare_beta %>% as_tibble() %>%
    relocate(SRE, .after = OLS)
  compare_names = colnames(compare_beta)
  compare_beta = compare_beta %>% select(!c(SpatialTP_fx, gSEM_fx, RSR_fx, RSR))
  # knitr::kable(t(colMeans(compare_beta)), 'latex', digits = 2)
  compare_scenarios = bind_rows(compare_scenarios, compare_beta)
}

compare_scenarios = bind_cols(compare_scenarios, Configuration = c(rep('I',n_sim), rep('II',n_sim), rep('III',n_sim))) %>%
  pivot_longer(cols = -Configuration, names_to = "Model", values_to = "Estimates") %>%
  mutate(Configuration = factor(Configuration, levels = c('I', 'II', 'III'),
                                labels = c(expression('Configuration I ('*phi[x]* '= 0.05, '*phi[w]*'= 0.5)'), 
                                           expression('Configuration II ('*phi[x]* '= 0.5, '*phi[w]*'= 0.05)'),
                                           expression('Configuration III ('*phi[x]* '= 0.2, '*phi[w]*'= 0.2)'))),
         Model = factor(Model, levels = unique(Model)))
model_names2 = c("OLS","SRE", "Spatial+_fx", "SpatialTP", "gSEM", "Spatial+", "KS",
                 "SA", "SS_fv", "SS_nmig", "SS_mom")
model_face = rep('plain', length(model_names2))
model_face[9:11] = 'bold'

#### THIS IS FIGURE 1 ####
fig2 = ggplot(compare_scenarios, aes(x = Model, y = Estimates)) +
  geom_boxplot(fill = "#D0DDD7") +
  geom_hline(yintercept=beta.real[ind], linetype="dashed",
             linewidth = 1.1, color = "red") +
  xlab("Models") + ylab("Estimates") +
  scale_x_discrete(limits=model_names2) +
  theme_bw() +
  facet_wrap(~Configuration, nrow = 3, scales = 'free_y', labeller = label_parsed) +
  theme(axis.text.x = element_text(face = model_face, angle = 90),
        panel.grid = element_blank())
ggsave(paste0(plot_path, "boxplot_3configurations.pdf"), fig2, width = 7.5, height = 5.5, device = "pdf")
#####################à


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

# THIS IS FIGURE 2 ####
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




n0 = 500
vec = c(7,3,5)
label_vec = c('I','II','III')# as.character(vec)
lr = ''
#vec = c(1:9, 37:45, 145:153)
coefs_mom = matrix(ncol = n0-1, nrow = length(vec))
rownames(coefs_mom) = paste0("config", vec)
coefs_george = coefs_nmig = coefs_mom
#coefs_mom = as_tibble(coefs_mom)
threshold = 0.5
cc = 1
for (settingplot in vec) {
  SS = readRDS(paste0(outpath, lr, "SS1/sim", settingplot,"_", type_dgm,".RDS"))
  coefs_george[cc,] = colMeans(SS$pip > threshold)
  SS = readRDS(paste0(outpath, lr, "SS_nmig1/sim", settingplot,"_", type_dgm,".RDS"))
  coefs_nmig[cc,] = colMeans(SS$pip > threshold)
  SS = readRDS(paste0(outpath, lr, "SS_mom1/sim", settingplot,"_", type_dgm,".RDS"))
  coefs_mom[cc,] = colMeans(SS$pip > threshold)
  
  cc = cc+1
}
coefs_george[coefs_george==0] = NA
coefs_nmig[coefs_nmig==0] = NA
coefs_mom[coefs_mom==0] = NA

df=data.frame(Coefficient = 1:(n0-1), Configuration='I', Prior = 'SS_fv',
              Value = coefs_george[1,])
df = rbind(df,
           data.frame(Coefficient = 1:(n0-1), Configuration='II', Prior = 'SS_fv',
                      Value = coefs_george[2,]),
           data.frame(Coefficient = 1:(n0-1), Configuration='III', Prior = 'SS_fv',
                      Value = coefs_george[3,]),
           data.frame(Coefficient = 1:(n0-1), Configuration='I', Prior = 'SS_nmig',
                      Value = coefs_nmig[1,]),
           data.frame(Coefficient = 1:(n0-1), Configuration='II', Prior = 'SS_nmig',
                      Value = coefs_nmig[2,]),
           data.frame(Coefficient = 1:(n0-1), Configuration='III', Prior = 'SS_nmig',
                      Value = coefs_nmig[3,]),
           data.frame(Coefficient = 1:(n0-1), Configuration='I', Prior = 'SS_mom',
                      Value = coefs_mom[1,]),
           data.frame(Coefficient = 1:(n0-1), Configuration='II', Prior = 'SS_mom',
                      Value = coefs_mom[2,]),
           data.frame(Coefficient = 1:(n0-1), Configuration='III', Prior = 'SS_mom',
                      Value = coefs_mom[3,]))

# THIS IS FIGURE 3 ####
gg1 = ggplot(df, aes(Coefficient, Value, color=Prior)) + 
  geom_line(linewidth=1.5) +
  # scale_fill_gradientn(colours = brewer.pal(n = 11, name = 'Spectral'),
  #                      na.value = "transparent", limits = c(0, 1),
  #                      breaks = seq(0, 1, 0.2)) +
  scale_color_discrete(na.value = "grey", breaks = c('SS_fv', 'SS_nmig', 'SS_mom')) +
  theme_bw() + facet_wrap(~Configuration, nrow=3, labeller = label_both) +
  theme(panel.grid = element_blank()) +
  labs(x = "Basis #", y = "Proportion")
gg1
ggsave(file.path(plot_path, 'proportion selected bases.pdf'), plot = gg1,
       width = 7.5, height = 5.5, device = pdf)

###########,












organize.data <- function(configs) {
  compare_scenarios = tibble()
  for (settingplot in configs) {
    DWA_ = readRDS(paste0(outpath, "DWA/sim", settingplot,"_", type_dgm,".RDS"))
    KS_ = readRDS(paste0(outpath, "KS/sim", settingplot,"_", type_dgm,".RDS"))
    subfolder = 'SS'
    SS_ = readRDS(paste0(outpath, lr, subfolder, options_model, "/sim", settingplot,"_", type_dgm,".RDS"))
    subfolder = 'SS_nmig'
    SS_nmig = readRDS(paste0(outpath, lr, subfolder, options_model, "/sim", settingplot,"_", type_dgm,".RDS"))
    subfolder = 'SS_mom'
    SS_mom = readRDS(paste0(outpath, lr, subfolder, options_model, "/sim", settingplot,"_", type_dgm,".RDS"))
    subfolder = 'SRE'
    SRE_ = readRDS(paste0(outpath, subfolder, "/sim", settingplot,"_", type_dgm,".RDS"))
    
    # compare_names<-c("OLS","SpatialTP_fx","RSR_fx","gSEM_fx","Spatial+_fx","SpatialTP",
    #                  "RSR","gSEM","Spatial+", "KS", "SS_fv", "SS_nmig", "SS_mom", "SRE")
    compare_names<-c("SS_fv", "SS_nmig", "SS_mom")
    # compare_beta = data.frame(DWA_$beta_results[,,'tprs150'], KS_$beta_results, SS_$pmeanbeta[,ind],
    #                           SS_nmig$pmeanbeta[,ind], SS_mom$pmeanbeta[,ind],
    #                           SRE_$beta_results)
    compare_beta = tibble(SS_$pmeanbeta[,ind], SS_nmig$pmeanbeta[,ind], SS_mom$pmeanbeta[,ind])
    compare_low = tibble(SS_$cilow[,ind], SS_nmig$cilow[,ind], SS_mom$cilow[,ind])
    compare_high = tibble(SS_$cihigh[,ind], SS_nmig$cihigh[,ind], SS_mom$cihigh[,ind])
    colnames(compare_beta) = compare_names
    high_minus_low = compare_high - compare_low
    colnames(high_minus_low) = paste0(compare_names, 'WIDTH')
    # compare_beta = compare_beta %>% relocate(SRE, .after = OLS)
    # compare_names = colnames(compare_beta)
    # compare_beta = compare_beta %>% select(!c(SpatialTP_fx, gSEM_fx, RSR_fx, RSR))
    compare_scenarios = bind_rows(compare_scenarios, bind_cols(compare_beta, high_minus_low))
  }
  compare_scenarios = bind_cols(compare_scenarios, Configuration = c(rep('I',n_sim), rep('II',n_sim), rep('III',n_sim))) %>%
    pivot_longer(cols = -Configuration, names_to = "Model", values_to = "Estimates") %>%
    mutate(Configuration = factor(Configuration, levels = c('I', 'II', 'III'),
                                  labels = c(expression('Configuration I ('*phi[x]* '= 0.05, '*phi[w]*'= 0.5)'), 
                                             expression('Configuration II ('*phi[x]* '= 0.5, '*phi[w]*'= 0.05)'),
                                             expression('Configuration III ('*phi[x]* '= 0.2, '*phi[w]*'= 0.2)'))),
           Model = factor(Model, levels = unique(Model)))
  return(compare_scenarios)
}

data_to_plot = bind_rows(
  bind_cols(organize.data(5 + c(2,-2,0)), 'Var.eps' = factor(0.25, levels = c(0.25,0.4,0.1))),
  bind_cols(organize.data(214 + c(2,-2,0)), 'Var.eps' = factor(0.4, levels = c(0.25,0.4,0.1))),
  bind_cols(organize.data(223 + c(2,-2,0)), 'Var.eps' = factor(0.1, levels = c(0.25,0.4,0.1)))
)

model_names2 = c("SS_fv", "SS_nmig", "SS_mom")
model_face = rep('bold', length(model_names2))

high_minus_low_plot = data_to_plot %>% filter(!Model %in% model_names2) %>% 
  group_by(Configuration, Var.eps, Model) %>% 
  summarise(mean_width = round(mean(Estimates),2))

# THIS IS FIGURE 4 ####
g1plus = ggplot(data_to_plot %>% filter(Model %in% model_names2), aes(x = Model, y = Estimates, fill = Var.eps)) +
  geom_boxplot(position=position_dodge(.8)) +
  geom_hline(yintercept=beta.real[ind], linetype="dashed",
             linewidth = 1.1, color = "red") +
  geom_label(aes(x=rep(1:3, 9), y = 1.9, label=mean_width), high_minus_low_plot,
             position=position_dodge(.8), show.legend = F, alpha=0.4) +
  xlab("Models") + ylab("Estimates") +
  scale_x_discrete(limits=model_names2) +
  ylim(1.8, NA) +
  theme_bw() +
  labs(fill=expression("Var["*epsilon[y]*"]")) + guides(col="none") +
  scale_fill_brewer(palette="Dark2") +
  # facet_grid(Configuration ~ ., scales='free', labeller = label_parsed, switch = 'y') +
  facet_wrap(~Configuration, nrow = 3, scales = 'free_y', labeller = label_parsed) +
  theme(axis.text.x = element_text(face = model_face, angle = 90),
        panel.grid = element_blank())
g1plus
ggsave(paste0(plot_path, "boxplot_3configurations_varyingSigmaeps.pdf"), g1plus, width = 7.5, height = 5.5, device = pdf)
###############.











