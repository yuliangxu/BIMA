

# generating a test example -----------------------------------------------
num_region = 1
side = 20 # > 50*50
p = side*side
region_idx = vector("list",num_region)
grids = GP.generate.grids(d=2L,num_grids=side)


idx_matr = matrix(1:(side*side),ncol = side)
side_per_region = side/sqrt(num_region)
for(r in 1:num_region){
  idx = rep(NA,(side_per_region)^2)
  colr = r - floor(r/sqrt(num_region))*sqrt(num_region);if(colr==0){
    colr = sqrt(num_region)
  }
  rowr = ceiling(r/sqrt(num_region));
  col_range = (max(colr-1,0):colr)*side_per_region;col_range[1] = col_range[1]+1;
  row_range = (max(rowr-1,0):rowr)*side_per_region;row_range[1] = row_range[1]+1;
  region_idx[[r]] = c(idx_matr[row_range[1]:row_range[2],col_range[1]:col_range[2]])
}



grids_df = as.data.frame(grids)
Soft_threshold = function(x,lambda){
  return( (x-sign(x)*lambda)*(abs(x)>lambda))
}


center = apply(grids,2,mean) + c(0.3,0.3)
rad = apply(grids,1,function(x){sum((x-center)^2)})
inv_rad = 2-rad
inv_rad_ST = Soft_threshold(inv_rad,1)
beta = 2*log(inv_rad_ST^2+1)
plot_img(beta,grids_df,"true beta")

center = apply(grids,2,mean) + c(-0.1,-0.1)
rad = apply(grids,1,function(x){sum((x-center)^2)})
inv_rad = 2-rad
inv_rad_ST = Soft_threshold(inv_rad,1.2)
alpha = 2*log(inv_rad_ST^2+1)
plot_img(alpha,grids_df,"true alpha")

plot_img(alpha*beta,grids_df,"true TIE")


alpha_true = alpha; beta_true = beta


lambda = 0.5
basis_sq = generate_sq_basis(grids, region_idx,a = 0.01,b = 10,poly_degree = 10)
L = length(unlist(basis_sq$Phi_D))

datsim = STGP_generate_theta_block(alpha_true,beta_true,
                                   2, 0.1, grids, basis_sq$Phi_Q,lambda=lambda,region_idx,basis_sq$L_all,
                                   n.sample=300)
datsim$lambda = lambda
save(datsim,file="./data/one_region.rda")
one_region = datsim
usethis::use_data(one_region)

# run BIMA ----------------------------------------------------------------

data(one_region)
BIMA_mcmc = BIMA(one_region$Y, one_region$X, one_region$M, one_region$C)



# show BIMA result --------------------------------------------------------


library(ggplot2)
library(viridis)
library(gtable)
library(grid)
# plot_img: Visualize 2D images
plot_img = function(img, grids_df,title="img",col_bar = NULL){
  ggplot(grids_df, aes(x=x1,y=x2)) +
    geom_tile(aes(fill = img)) +
    scale_fill_viridis_c(limits = col_bar, oob = scales::squish)+

    ggtitle(title)+
    theme(plot.title = element_text(size=20),legend.text=element_text(size=10))
}

e1=plot_img(apply(BIMA_mcmc$beta_sample,1,mean),as.data.frame(BIMA_mcmc$grids),"est beta_mean")
e2=plot_img(apply(BIMA_mcmc$alpha_sample,1,mean),as.data.frame(BIMA_mcmc$grids),"est alpha_mean")
e3=plot_img(apply(BIMA_mcmc$TIE_sample,1,mean),as.data.frame(BIMA_mcmc$grids),"est TIE_mean")


t1=plot_img(one_region$beta_test_ST,as.data.frame(BIMA_mcmc$grids),"true beta_mean")
t2=plot_img(one_region$alpha_test_ST,as.data.frame(BIMA_mcmc$grids),"true alpha_mean")
t3=plot_img(one_region$beta_test_ST*one_region$alpha_test_ST ,as.data.frame(BIMA_mcmc$grids),"true TIE_mean")

ge1 <- ggplotGrob(e1); ge2 <- ggplotGrob(e2); ge3 <- ggplotGrob(e3)
gt1 <- ggplotGrob(t1); gt2 <- ggplotGrob(t2); ge3 <- ggplotGrob(e3)
est_plot <- cbind(ggplotGrob(e1), ggplotGrob(e2),ggplotGrob(e3), size = "first")
true_plot <- cbind(ggplotGrob(t1), ggplotGrob(t2),ggplotGrob(t3), size = "first")
g = rbind( est_plot, true_plot)
# grid.newpage()
grid.draw(g)
