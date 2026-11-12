# From the repository root: Rscript Data_Analysis/plot_revision.R
# Optional arguments: RESULTS_DIR FIGURE_DIR SHAPEFILE
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) stop("Usage: Rscript Data_Analysis/plot_revision.R [RESULTS_DIR [FIGURE_DIR [SHAPEFILE]]]")
app <- if (length(args) >= 1L) args[[1L]] else "Data_Analysis/Results"
fig <- if (length(args) >= 2L) args[[2L]] else "Data_Analysis/Figures"
shapefile <- if (length(args) >= 3L) args[[3L]] else {
  "Data_Analysis/Data/county_boundaries/cb_2020_us_county_500k.shp"
}
dir.create(fig, recursive = TRUE, showWarnings = FALSE)
effects <- read.csv(file.path(app,'effects_all.csv'),check.names=FALSE)
fits <- readRDS(file.path(app,'application_fits.rds'))
regions <- c('USA','Midwest','Northeast','Southeast','West')
scenarios <- c('Main','Alternative contrast','Binary B (distance decay retained)','Without mobility')
fn <- c('plot_poli_main.pdf','plot_poli_supp_int.pdf','plot_poli_supp_adj.pdf','plot_poli_sens_mobi.pdf')
stopifnot(nrow(effects)==120L, identical(names(fits),regions),
          all(table(effects$Region,effects$Scenario)==6L))
for (region in regions) {
  run <- fits[[region]]
  stopifnot(all(effects$N[effects$Region==region]==nrow(run$data)),
            !any(run$data$STATEFP=='31'))
  for (k in names(run$effects)) {
    s <- run$effects[[k]]$results
    f <- run$fits[[if (k=='binary') 'Binary_B' else 'State_weighted']][[s$Mobility_Adjustment[1]]]
    shift <- s$Neighbor_Target[1] - s$Neighbor_Reference[1]
    b <- f$fit_y$coef; g <- f$fit_m$coef
    exact <- c(shift*f$h*b['WA'],shift*f$h*b['M']*g['WA'])
    stopifnot(max(abs(s$Estimate[s$Effect %in% c('tau4','tau5')]-exact))<1e-12)
    shown <- effects[effects$Region==region & effects$Scenario==s$Scenario[1],]
    cols <- c('Estimate','SE','Lower_95','Upper_95')
    stopifnot(max(abs(as.matrix(shown[cols])-as.matrix(s[cols])))<1e-12)
  }
}
effects$Region <- factor(effects$Region,levels=regions)
effects$Effect <- factor(effects$Effect,levels=paste0('tau',1:6),labels=paste('Effect',1:6))
# Use the same y scale within a pathway across the four figures.
axis_anchors <- do.call(rbind,lapply(split(effects,effects$Effect),function(z) {
  data.frame(Effect=z$Effect[1],Region=factor('USA',levels=regions),
             Estimate=range(c(0,z$Lower_95,z$Upper_95)))
}))
palette <- setNames(c('#0072B2','#D55E00','#009E73','#A84E91','#B16C00'),regions)
for (i in seq_along(scenarios)) {
  z <- effects[effects$Scenario==scenarios[i],]
  p <- ggplot(z,aes(x=Region,y=Estimate,colour=Region)) +
    geom_blank(data=axis_anchors,aes(x=Region,y=Estimate),inherit.aes=FALSE) +
    geom_hline(yintercept=0,linetype='dashed',linewidth=.45,colour='grey55') +
    geom_errorbar(aes(ymin=Lower_95,ymax=Upper_95),width=.16,linewidth=.62) +
    geom_point(size=1.85) +
    facet_grid(rows=vars(Effect),scales='free_y') +
    scale_colour_manual(values=palette,guide='none') +
    labs(x='Region',y='Effect on square-root deaths per 100 residents') +
    theme_bw(base_size=10.5) + theme(panel.grid.minor=element_blank(),
      panel.grid.major.x=element_blank(),panel.spacing=grid::unit(7,'pt'),
      strip.background=element_rect(fill='grey94',colour='grey74'),
      strip.text.y=element_text(face='bold'),axis.text.x=element_text(colour='black'))
  ggsave(file.path(fig,fn[i]),p,width=7.4,height=8.4,device=grDevices::pdf,useDingbats=FALSE)
}
# Map the original-scale variables for the same counties used in the USA fit.
geo <- sf::st_read(shapefile, quiet = TRUE)
conus <- c('01','04','05','06','08','09','10','11','12','13','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','42','44','45','46','47','48','49','50','51','53','54','55','56')
geo <- geo[geo$STATEFP %in% conus,]
fips <- paste0(geo$STATEFP,geo$COUNTYFP)
dat <- fits$USA$data
stopifnot(all(dat$FIPS %in% fips),!anyDuplicated(dat$FIPS),length(dat$FIPS)==2906L)
idx <- match(fips,dat$FIPS)
for (v in c('A','M_raw','Y_raw','C1','C2','C3')) geo[[v]] <- dat[[v]][idx]
geo$Deaths_per_100 <- 100*geo$Y_raw
map <- function(var,title) {
  ggplot(geo) + geom_sf(aes(fill=.data[[var]]),colour=NA) +
    scale_fill_viridis_c(option='plasma',na.value='grey89',guide=guide_colourbar(
      title=NULL,barwidth=grid::unit(.15,'in'),barheight=grid::unit(.95,'in'))) +
    scale_x_continuous(breaks=seq(-120,-70,by=10)) +
    scale_y_continuous(breaks=seq(25,50,by=5)) +
    coord_sf(crs=sf::st_crs(4326),datum=sf::st_crs(4326),
      xlim=c(-125,-66),ylim=c(24,50),expand=FALSE) +
    labs(title=title,x=NULL,y=NULL,fill=NULL) +
    theme_minimal(base_size=11) + theme(
      axis.text=element_text(colour='grey25',size=10),
      panel.grid.major=element_line(colour='grey90',linewidth=.38),
      panel.grid.minor=element_blank(),
      legend.position='right',legend.text=element_text(colour='black',size=10),
      legend.margin=margin(0,0,0,0),legend.box.spacing=grid::unit(7,'pt'),
      plot.title=element_text(hjust=.5,size=13,margin=margin(b=5)),
      plot.margin=margin(5,8,4,6))
}
ps <- list(map('A','(A) Democratic Votes'),map('M_raw','(B) Vaccine Hesitancy'),
  map('Deaths_per_100','(C) COVID-19 Mortality'),map('C1','(D) Racial Distribution'),
  map('C3','(E) Social Vulnerability Index'),map('C2','(F) Mobility'))
ggsave(file.path(fig,'map_plot.pdf'),patchwork::wrap_plots(ps,ncol=2),
       width=9.6,height=7.4,device=grDevices::pdf,useDingbats=FALSE)
q <- rbind(data.frame(Variable='Outcome: square-root mortality',Value=dat$Y),
           data.frame(Variable='Mediator: square-root hesitancy',Value=dat$M))
qq <- ggplot(q,aes(sample=Value)) + stat_qq(size=.68,alpha=.55,colour='#0072B2') +
  stat_qq_line(linetype='dashed',colour='grey25') + facet_wrap(~Variable,scales='free_y',nrow=1) +
  labs(x='Theoretical normal quantiles',y='Observed quantiles') +
  theme_bw(base_size=10.5) + theme(panel.grid.minor=element_blank(),strip.text=element_text(face='bold'))
ggsave(file.path(fig,'qqplot_sqrt.pdf'),qq,width=7.4,height=3.25,device=grDevices::pdf,useDingbats=FALSE)
cat('Verified all 120 effects, isolate factors, and matched map/QQ cohort. Wrote six vector figures.\n')
