library(FrF2)
setwd("C:/Users/stella/Documents/GitHub/monica-germany-q4b/factorial_design/")

##FUNCTIONS
main_effect_plot<-function(outfile, response){
  
  png(filename=outfile, width=1400, height=1200)
  
  par(mar=c(15, 7, 4.1, 2.1))
  #TODO: set correct labels!
  #x_labels<-c("GW_Level", "Imp_Layer", "Pheno", "Slope", "Sowing", "Harvest", "LandCover", "Nresponse", "Wresponse", "Yield_Calib")
  mar=c(10,1,1,1)
  plot(response, cex = 2.5,
       yaxt="n", ylab="",
       xaxt="n", xlab="",
       main = "")
  
  #mgp=c(axis.title.position, axis.label.position, axis.line.position))
  axis(1, cex.axis=3, mgp=c(0,1,0), at=1:10, labels=x_labels, las=2)
  axis(2, cex.axis=3, mgp=c(0,1,0), las=1)
  
  dev.off()
}

fancy_main_effect_plot<-function(outfile, response){
  
  png(filename=outfile, width=1400, height=900)
  MEPlot(response, abbrev = 5, cex.xax = 1.6, cex.main = 2, main = "")
  dev.off()
}

daniel_plot<-function(outfile, response){
  png(filename=outfile, width=1400, height=900)
  DanielPlot(response, code = TRUE, half = TRUE, alpha = 0.1,
             cex.main = 2, cex.pch = 1.5, cex.lab = 1.5, cex.fac = 1.5,
             cex.axis = 1.5, main = "Half Normal Plot r, alpha=0.1")
  dev.off()
}
###END FUNCTIONS

design <- FrF2(resolution = 5, randomize = FALSE, factor.names = list(
  #GroundWaterLevel = c("false", "true"), 
  #ImpenetrableLayer = c("false", "true"),
  SowingDate = c("de", "lk"),
  #HarvestDate = c("de", "auto"),
  Slope = c("false", "true"),
  LandCover = c("false", "true"), 
  WaterDeficitResponse = c("false", "true"),
  Nresponse_and_Fertil = c("false", "true"), 
  Phenology_cal = c("de", "lk"),
  Yield_cal = c("de", "lk")
  ))

design
summary(design)

export.design(design, filename = "design_res_V_new", type = "csv", OutDec = ".") 

res_r <- read.csv("C:/Users/stella/Documents/GitHub/monica-germany/calculate-indices/report_best_exps/r.csv",header=F)
res_rrmse <- read.csv("C:/Users/stella/Documents/GitHub/monica-germany/calculate-indices/report_best_exps/RRMSE.csv",header=F)
res_pbias <- read.csv("C:/Users/stella/Documents/GitHub/monica-germany/calculate-indices/report_best_exps/pBIAS.csv",header=F)
res_ai <- read.csv("C:/Users/stella/Documents/GitHub/monica-germany/calculate-indices/report_best_exps/AI.csv",header=F)

resp_r <- add.response(design, res_r)
resp_rrmse <- add.response(design, res_rrmse)
resp_pbias <- add.response(design, res_pbias)
resp_ai <- add.response(design, res_ai)

#removing leaf Ext modifier... (for the presentation)
resp_r$LeafExtensionModifier = NA
resp_rrmse$LeafExtensionModifier = NA
resp_pbias$LeafExtensionModifier = NA
resp_ai$LeafExtensionModifier = NA

#original plot:
#plot(resp_rrmse, cex = 1.5, cex.lab = 1.5, cex.axis = 1.2,
#     main = "Main effects plot rrmse", cex.main = 1.5)
#####

main_effect_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/r.png", response=resp_r)
main_effect_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/rrmse.png", response=resp_rrmse)
main_effect_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/pbias.png", response=resp_pbias)
main_effect_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/ai.png", response=resp_ai)

#fancy_main_effect_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/r_f.png", resp_r)
daniel_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/r_DanielPlot.png", resp_r)
daniel_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/rrmse_DanielPlot.png", resp_rrmse)
daniel_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/pbias_DanielPlot.png", resp_pbias)
daniel_plot("C:/Users/stella/Desktop/ZALF/ESA_2018/ai_DanielPlot.png", resp_ai)

#IAPlot(resp_rrmse, abbrev = 5, show.alias = TRUE, lwd = 2, cex = 2,
#       cex.xax = 1.2, cex.lab = 1.5)



#Definition of high elevation and shallow groundwater landkreise
topography <- read.csv("C:/Users/stella/Documents/GitHub/monica-germany/avg_elevation_latitude_gw_per_landkreis.csv",header=T)
elevations <- as.vector(topography['elevation'])
gwater <- as.vector(topography['groundwaterlevel'])

my_var = unlist(gwater)

hist(my_var, breaks=5, col="red")

x <- my_var 
h<-hist(x, breaks=10, col="red",  
        main="Histogram with Normal Curve") 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)

#example of RFD below
plan.annotated <- FrF2(16, 6, factor.names = list(
  DieOrif = c(2.093, 2.1448), PistDiam = c(9.462, 9.5),
  Temp = c(188.1, 191.1), DieClean = c("Dirty", "Clean"),
  SMass = c(4, 8), BarClean = c("Dirty", "Clean")),
  seed = 6285)

summary(plan.annotated)

MI <- c(35.77, 35.03, 38.5, 39.33, 35.7, 35.1, 39.27, 37, 41.07, 32.03,
           42, 37.63, 40.2, 37, 40.1, 35.03)

plan.resp <- add.response(plan.annotated, MI)

plot(plan.resp, cex = 1.2, cex.lab = 1.2, cex.axis = 1.2,
     main = "Main effects plot for MI", cex.main = 2)

#MEPlot(plan.resp, abbrev = 5, cex.xax = 1.6, cex.main = 2)

IAPlot(plan.resp, abbrev = 5, show.alias = TRUE, lwd = 2, cex = 2,
       cex.xax = 1.2, cex.lab = 1.5)

IAPlot(plan.resp, abbrev = 5, select = 3:6, lwd = 2, cex = 2,
       cex.xax = 1.2, cex.lab = 1.5)

summary(lm(plan.resp))

DanielPlot(plan.resp, code = TRUE, half = TRUE, alpha = 0.1,
           cex.main = 1.8, cex.pch = 1.2, cex.lab = 1.4, cex.fac = 1.4,
           cex.axis = 1.2)

summary(plan.resp)
