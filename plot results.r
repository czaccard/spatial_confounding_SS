library(tidyverse)
library(geoR)
library(mvtnorm)
library(metR)
library(paletteer)
beta.real = c(1,2) %>% as.matrix(ncol=1) # vector giving intercept coefficient and beta_x
ind = 2 # 2: indicates the exposure coefficient in a design matrix like cbind(1,X)

outpath = "results/"
inpath = "func/"
plot_path = "plots/"


type_dgm = "RelBias_deltafixed"
options_model = 1 # choose 1, 2, or 3

n_sim = 100


# FIGURE 1 #####
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

#### Plot FIGURE 1 ####
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
ggsave(paste0(plot_path, "difference_part1_minus_part2_basis_from_1_to_k.pdf"),
       width = 7.5, height = 4, device = "pdf")


# FIGURE 2 ######
all.res = readRDS(paste0(outpath, '/ratio.models.OLS.RDS'))
model_names = c('SRE', 'Spatial+_fx', 'SpatialTP', 'gSEM', 
                'Spatial+', 'KS', 'SA', 'SS_mom')

type_plot = 'RMSE' # 'RMSE' or 'MAE' or 'Estimate'

ggdf2 = with(all.res, {
  a = select(setting2, range1, range2) %>% 
    rename(range_x = range1, range_w = range2)
  ggdf = as_tibble(lapply(all.res[-1], function(x) {
    if (type_plot=='RMSE') sqrt(rowMeans((x-2)^2))/sqrt(rowMeans((betaOLS-2)^2))
    else if (type_plot=='Estimate') rowMeans(x)/rowMeans(betaOLS)
    else if (type_plot=='MAE') rowMeans(abs(x-2))/rowMeans(abs(betaOLS-2))
  }))
  
  ggdf = a %>% bind_cols(ggdf) %>% select(-betaOLS)
  ggdf
})



ggdf3 = ggdf2 %>%
  pivot_longer(-c(range_x, range_w), names_to = 'Model', values_to = 'Model.over.OLS') %>%
  mutate(smaller1=factor(Model.over.OLS<1, levels = c('FALSE', 'TRUE')),
         Model=factor(Model, 
                      levels = c('betaSRE', 'betaSPPLUS_fx','betaSPATIALTP',
                                 'betaGSEM','betaSPPLUS','betaKS', 'betaSA',
                                 'betaPMOM'),
                      labels = model_names))


if (type_plot=='Estimate') {breaks = seq(0.94, 1.28, by=0.02); ttl = '(a)'}
if (type_plot=='RMSE') {breaks = seq(0.6, 3, by=0.2); ttl = '(b)'}
if (type_plot=='MAE') {breaks = seq(0.6, 3.2, by=0.2); ttl = '(a)'}
gg2 = ggplot(ggdf3, 
             aes(range_x, range_w, z = Model.over.OLS)) +
  geom_contour_filled(breaks = breaks) +
  geom_contour(breaks = breaks, color='black') +
  geom_contour(breaks = 1, color='black', linewidth=1) +
  geom_line(aes(range_x, range_w), filter(ggdf3, range_x==range_w),
            color='grey50', linewidth=1, linetype = 'dotted') +
  geom_text_contour(breaks = breaks, stroke = 0.15, skip = 0, size = 2.8,
                    label.placer = label_placer_fraction(), check_overlap = T) +
  facet_wrap(~Model, nrow = 2) +
  scale_fill_manual(
    values = c(
      paletteer_c("ggthemes::Blue", sum(breaks<1), dir=-1),
      paletteer_c("ggthemes::Red", sum(breaks>=1)))) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab(expression(paste('Range Exposure (',phi[x], ')'))) +
  ylab(expression(paste('Range Confounder (',phi[w], ')'))) +
  ggtitle(ttl) +
  theme(legend.position = 'none', panel.grid = element_blank(),
        title = element_text(size = 8),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 7),
        strip.text = element_text(size = 7),
        plot.margin=unit(c(0,0.4,0,0), "lines"),
        panel.spacing = unit(0.1,'lines'))

if (type_plot=='RMSE') {ggr = gg2; ggname = "ratios_RMSE.pdf"}
if (type_plot=='Estimate') {gge = gg2; ggname = "ratios_Estimate.pdf"}
if (type_plot=='MAE') {ggb = gg2; ggname = "ratios_MAE.pdf"}


library(gridExtra)
gg1 = arrangeGrob(ggb, ggr, nrow=2, ncol=1)
plot(gg1)


ggsave(paste0(plot_path, "ratios_maps.pdf"), gg1, 
       width = 6, height = 8, device = pdf)


ggdf3 %>% group_by(Model) %>% 
  summarise(m = mean(as.logical(smaller1)))

val1 = 0.98*(type_plot=='Estimate') + 0.8*(type_plot=='RMSE') + 0.8*(type_plot=='MAE')
val2 = 1.1*(type_plot=='Estimate') + 1.8*(type_plot=='RMSE') + 1.8*(type_plot=='MAE')
ggdf3 %>% group_by(Model) %>% filter(range_x<range_w) %>% 
  summarise(m = mean(Model.over.OLS<val1))

ggdf3 %>% group_by(Model) %>% filter(range_x<range_w & 0.2<range_x) %>% 
  summarise(m = mean(Model.over.OLS<val1))

ggdf3 %>% group_by(Model) %>%
  summarise(m = mean(Model.over.OLS>val2))



ggdf4 = with(all.res, {
  a = select(setting2, range1, range2) %>% 
    rename(range_x = range1, range_w = range2)
  ggdf = as_tibble(lapply(all.res[-1], function(x) {
    rowMeans(abs(x) < abs(betaOLS))
  }))
  
  ggdf = a %>% bind_cols(ggdf) %>% select(-betaOLS)
  ggdf
})
colMeans(ggdf4)
colMeans(ggdf4 %>% filter(range_x<range_w))
colMeans(ggdf4 %>% filter(range_x<range_w & 0.2<range_x))
