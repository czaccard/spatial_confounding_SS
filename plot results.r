library(tidyverse)
library(geoR)
library(mvtnorm)
beta.real = c(1,2) %>% as.matrix(ncol=1) # vector giving intercept coefficient and beta_x
ind = 2 # 2: indicates the exposure coefficient in a design matrix like cbind(1,X)

outpath = "results/"
inpath = "func/"
plot_path = "plots/"


type_dgm = "RelBias_deltafixed"
options_model = 1 # choose 1, 2, or 3

n_sim = 100


# Figures in the paper #####
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

#### THIS IS FIGURE 2 ####
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







source(paste0(inpath, "data_generate.r"), local = T)
source(paste0(inpath, "misc.r"), local = T)
setting <- readxl::read_excel(paste0(inpath, "setting.xlsx"), sheet = "Foglio2")
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
mu1 = rep(0, n0+n_new)
mu2 = rep(0, n0+n_new)
n_sim = 100
type_dgm = "conditional"
analysis1 = function(k, B, options_model) {
  nk = length(k)
  df = matrix(0, nk, 2, dimnames = list(paste0("d", k), c("part1", "part2")))
  df = data.frame(df)
  
  switch (as.character(options_model),
          '1' = {
            for (ii in seq(nk)) {
              row = paste0("d", k[ii])
              ind1k = if(k[ii]>0) 1:k[ii] else numeric(0)
              bb = B[, ind1k]
              if (NCOL(bb)!=0) {
                s = crossprod(bb) - t(bb) %*% temp2 %*% bb
                invs = solve(s)
                t = temp1 %*% bb %*% invs
                temp4 = DG$delta * sqrt(DG$sigma2w/DG$sigma2z) * t %*% t(bb)
                p1 = temp4 %*% temp2 %*% temp3
                p2 = temp4 %*% temp3
              } else {p1 = p2 = c(0,0)}
              
              df[row,] = c(p1[2], p2[2])
              cat('.')
            }
          },
          '2' = {
            temp4 = DG$delta * sqrt(DG$sigma2w/DG$sigma2z)
            temp5 = C %*% solve(crossprod(C)) %*% t(C)
            invCXmat = C %*% solve(crossprod(C, diag(n0) - temp2) %*% C) %*% t(C)
            invXCmat = solve(crossprod(xa, diag(n0) - temp5) %*% xa)
            p1 = temp4 * temp1 %*% invCXmat %*% temp2 %*% temp3
            p2 = temp4 * invXCmat %*% t(xa) %*% temp5 %*% temp3
            df$part1 = rep(p1[2], nk) %>% as.numeric()
            df$part2 = rep(p2[2], nk) %>% as.numeric()
          },
          '3' = {
            foo = 1
          }
  )
  df = df %>% mutate(difference = part1 - part2, n_bases = k, options_model = as.factor(options_model))
  return(as_tibble(df))
}
K = c(0, 1, 5, 10, 20, 50, 100, 150, 200, 250, 300, 350, 400, 450, 490, 495, 497)
C = unname(as.matrix(coords[1:n0,]))

# this takes few minutes...
out = data.frame()
for (rz in c(0.5, 0.2, 0.05)) {
  for (rg in c(0.5, 0.2, 0.05)) {
    j = (setting %>% filter(range1==rz, range2==rg, delta==0.5, sigma_1==1, sigma_2==1))$sim
    DG = generate_data(setting, j, type_dgm, coords, D, 0.15)
    x = DG$X[1,1:n0]
    xa = cbind(1,x)
    temp1 = solve(crossprod(xa)) %*% t(xa)
    temp2 = xa %*% solve(crossprod(xa)) %*% t(xa)
    temp_gls = DG$Lw[1:n0,1:n0] %*% solve(DG$Lz[1:n0,1:n0])
    temp3 = temp_gls %*% x
    p_formula = ~lngcoords+latcoords
    psplines = principal.spline.2D(coords, 1:n0, p_formula, x=DG$X[1,])
    E = psplines$F_n[,-(1:psplines$q)]
    df1 = analysis1(K, cbind(C, E), 1)
    p_formula = ~x+lngcoords+latcoords
    psplines = principal.spline.2D(coords, 1:n0, p_formula, x=DG$X[1,])
    E = psplines$F_n[,-(1:psplines$q)]
    df2 = analysis1(K, NA, 2)
    p_formula = ~x
    psplines = principal.spline.2D(coords, 1:n0, p_formula, x=DG$X[1,])
    E = psplines$F_n[,-(1:psplines$q)]
    df3 = analysis1(K, NA, 3)
    biasOlsGls = c(
      bias_ols(setting$sigma_2[j], setting$delta[j], DG$X[1,1:n0], setting$sigma_1[j], temp_gls),
      bias_gls(setting$sigma_2[j], setting$delta[j], DG$R_phiw[1:n0,1:n0], 0.25, DG$X[1,1:n0], setting$sigma_1[j], temp_gls)
    )
    df = rbind(df1, df2, df3) %>% mutate(rz = rz, rg = rg, diff_gls = diff(biasOlsGls))
    out = rbind(out, df)
    cat('\n')
  }
}

out_d = out %>% dplyr::select(n_bases, difference, diff_gls, rz, rg, options_model)

#### THIS IS FIGURE 1 ####
g5 = ggplot(out_d, aes(n_bases)) + 
  geom_line(aes(y=difference, linetype=options_model), linewidth=.6) +
  # geom_line(aes(y=diff_gls), color='darkgray', linetype='dotted') +
  facet_grid(rz ~ rg, scales='free',
             labeller=label_bquote(phi[x] == .(rz), phi[w] == .(rg))) +
  theme_bw() +
  xlab(latex2exp::TeX("Number of bases (\\it{$k$})")) +
  ylab(latex2exp::TeX("Difference (\\it{$d_x$})")) +
  # scale_color_manual(name = "Null-space\ntype", 
  #                    values=c('#AAAEB6', '#484B54', '#000000')) +
  scale_linetype_manual(name = "Null-space\ntype", 
                        values=c("solid", "dotdash", "longdash")) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
g5
ggsave(paste0(plot_path, "difference part1 minus part2_basis from 1 to k.pdf"),
       width = 7.5, height = 4, device = "pdf")



