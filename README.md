#https://eurekastatistics.com/generating-random-survival-times-from-any-hazard-function/
require(survival)
require(rms)
require(ggplot2)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inverse = function(fn, min_x, max_x){
  # Returns the inverse of a function for a given range.
  # E.g. inverse(sin, 0, pi/2)(sin(pi/4)) equals pi/4 because 0 <= pi/4 <= pi/2
  fn_inv = function(y){
    uniroot((function(x){fn(x) - y}), lower=min_x, upper=max_x)[1]$root
  }
  return(Vectorize(fn_inv))
}

integrate_from_0 = function(fn, t){
  int_fn = function(t) integrate(fn, 0, t)
  result = sapply(t, int_fn)
  value  = unlist(result["value",])
  msg    = unlist(result["message",])
  value[which(msg != "OK")] = NA
  return(value)
}

random_survival_times = function(hazard_fn, n, max_time=10000){
  # Given a hazard function, returns n random time-to-event observations.
  cumulative_density_fn = function(t) 1 - exp(-integrate_from_0(hazard_fn, t))
  inverse_cumulative_density_fn = inverse(cumulative_density_fn, 0, max_time)
  return(inverse_cumulative_density_fn(runif(n)))
} 

log(4/7)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot some curves
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1st curve
median.time<-4
lambda<-log(2)/median.time
n =100
hazard_fn = function(t) rep(lambda, length(t))
survival_times = random_survival_times(hazard_fn, n) 
survfit(Surv(survival_times) ~ 1)

require("survival")
plot(survfit(Surv(survival_times, rep(1, length(survival_times)))~1), xlab="Time", ylab="Survival Probability",
     main="Sampled and Expected Survival Curves for h(t) = lambda")
neg_exp_fn = function(x){exp(-lambda * x)}
curve(expr=neg_exp_fn, from=0, to=max(survival_times), add=TRUE, col="red", lwd=2)

##2nd curve
median.time1<-7
lambda1<-log(2)/median.time1
hazard_fn = function(t) rep(lambda1, length(t))
survival_times2 = random_survival_times(hazard_fn, n) 
survfit(Surv(survival_times2) ~ 1)

require("survival")
plot(survfit(Surv(survival_times2, rep(1, length(survival_times2)))~1), xlab="Time", ylab="Survival Probability",
     main="Sampled and Expected Survival Curves for h(t) = lambda1")
neg_exp_fn = function(x){exp(-lambda1 * x)}
curve(expr=neg_exp_fn, from=0, to=max(survival_times2), add=TRUE, col="red", lwd=2)

##dataframe
d <- as.data.frame(cbind(grp = gl(2, n, labels = c("Control", "Treat"))  , 
                         time=c(survival_times,survival_times2 )))

 
#analyse, cox score test pvalue
logR <- survdiff(Surv(time ) ~ grp,data=d)
(coxphobject<-coxph(Surv(time)~grp,d))
summcph <- summary(coxphobject)
summcph$sctest


#expected HR
median.time/ median.time1
log(median.time/ median.time1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot together
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


require(rms)
f2 <- npsurv(Surv(time ) ~ grp, d)

x <- c("#e41a1c","#377eb8","#4daf4a","#984ea3")

survplot(f2,  n.risk=TRUE, levels.only=T, conf.int=T,
         aehaz=TRUE,
         conf=c("bands"), 
         col.fill=gray(seq(.95, .75, length=4)),
         #col.fill= c(rgb(0,0,0,0.1)), 
         type=c("kaplan-meier"),
         lty=1,  col=x, xlab="Months", abbrev.label=T, 
         label.curves = list(keys = "lines"), #bty='n',
         y.n.risk= 0, cex.n.risk=.6, time.inc=2)


#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### make a function to test this statement...
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# With at least 100 patients per group, 
# there was approximately 94% power to detect a 3-day difference 
# in median time to otorrhea cessation between Moxidex plus TT and TT only. 
# This was based on the assumption that the median time to cessation of otorrhea 
# for patients receiving Moxidex plus TT would be similar to the 4-day median time 
# to cessation of otorrhea observed in all Alcon-conducted studies of CIPRODEX.
# All power and sample size calculations were based on the assumption of an 
# exponential survival distribution for the 2 groups and the application of a 2-sided
# test at the a = 0.05 level of significance. 


expsim <- function(n=100, lambda1 =log(2)/4, lambda2=log(2)/7) {
  
  # 1st curve
  # lambda1 
  hazard_fn = function(t) rep(lambda1, length(t))
  survival_times = random_survival_times(hazard_fn, n) 
   
  # 2nd curve
  # lambda2 
  hazard_fn = function(t) rep(lambda2, length(t))
  survival_times2 = random_survival_times(hazard_fn, n) 
   
  # dataframe
  d <- as.data.frame(cbind(grp = gl(2, n, labels = c("Control", "Treat"))  , 
                           time=c(survival_times,survival_times2 )))

  # analyse, score test for Cox
  # logR <- survdiff(Surv(time ) ~ grp,data=d)
  coxphobject <- coxph(Surv(time)~grp,d)
  summcph <- summary(coxphobject)
  return((summcph$sctest[3][[1]]))
  
}


x <-  replicate(299 , expsim(n=25)) # n is in each group
mean(x<0.05)

## good approx to this?
# m expected total number of events over both groups.
powerSurvEpi::powerCT.default0(k = 1, m = 50, RR = 4/7, alpha = 0.05)

# this seems to agree, put power from the above
coef<- log(4/7)
d = (qnorm(.975) + qnorm(.49))^2 /(.5*.5*coef^2)
d


######################################################################################
######################################################################################
######################################################################################

expsim2 <- function(n=25, lambda1 =log(2)/4, lambda2=log(2)/7, lambda3=log(2)/1) {
  
  # 1st curve
  #lambda1 
  hazard_fn = function(t) rep(lambda1, length(t))
  survival_times = random_survival_times(hazard_fn, n) 
  
  ##2nd curve
  #lambda2 
  hazard_fn = function(t) rep(lambda2, length(t))
  survival_times2 = random_survival_times(hazard_fn, n) 
  
  ##dataframe
  d <- as.data.frame(cbind(grp = gl(2, n, labels = c("Control", "Treat"))  , 
                           time=c(survival_times,survival_times2 )))
  
  #analyse
  #logR <- survdiff(Surv(time ) ~ grp,data=d)
  coxphobject <- coxph(Surv(time)~grp,d)
  summcph <- summary(coxphobject)
  pv1 <- summcph$sctest[3][[1]]
  
  
  ##3nd curve
  hazard_fn = function(t) rep(lambda3, length(t))
  survival_times3 = random_survival_times(hazard_fn, n) 
  
  ##dataframe
  d <- as.data.frame(cbind(grp = gl(2, n, labels = c("Control", "Treat"))  , 
                           time=c(survival_times,survival_times3 )))
 
  #analyse
  #logR <- survdiff(Surv(time ) ~ grp,data=d)
  coxphobject <- coxph(Surv(time)~grp,d)
  summcph2 <- summary(coxphobject)
  pv2 <- summcph2$sctest[3][[1]]
  pv2
  return(min(pv1, pv2))
  
}

x <-  replicate(50  , expsim2())
mean(x<0.05)
##################


#/***************************************************************************#
#Filename           :      xxxx.Rmd
#Author             :      obrieea1
#Date               :      xx-xxx-xxxx
#R                  :      3.2.3 (2015-12-10) 
#Platform           :      x86_64-pc-linux-gnu (64-bit) 
#Project/Study      :      CIGE025E3401\pub_2  [example]
#Description        :      xxxxxx
#Assumptions:   
#Input              :      xxxxxxx
#Output             :      
#Macros used        :      none
#---------------------------------------------------------------------------
#MODIFICATION HISTORY: 
#    <DD-MON-YYYY>,
#    <Description> 
#    
#***************************************************************************/

        rm(list=ls())
 
        set.seed(1234)
        startTime<-proc.time()
        library(knitr)
        options(width=120)
        opts_chunk$set(comment = "", warning = FALSE, message = FALSE,
                       echo = FALSE, tidy = FALSE, size="tiny",  cache=FALSE,
                       progress=TRUE,
                       #out.width='500px', dpi=200,
                       fig.width=6, fig.height=4,
                       cache.path = 'program_Cache/',
                       fig.path='figure/')
        
        #knitr::opts_knit$set(global.par = TRUE)
        #opts_chunk$set(out.width='750px', dpi=200) 
        
        knitr::knit_hooks$set(inline = function(x) {
          knitr:::format_sci(x, 'md')
        })
        
                
        # create an R file of the code!
        # https://stackoverflow.com/questions/26048722/knitr-and-tangle-code-without-execution
        
         knit_hooks$set(purl = function(before, options) {
           if (before) return()
           input  = current_input()  # filename of input document
           output = paste(tools::file_path_sans_ext(input), 'R', sep = '.')
           if (knitr:::isFALSE(knitr:::.knitEnv$tangle.start)) {
           assign('tangle.start', TRUE, knitr:::.knitEnv)
           unlink(output)
         }
           cat(options$code, file = output, sep = '\n', append = TRUE)
         })

         

          # where<-"blah"
          # 
          # x  <- "data/R Scripts/PROJECTS/COPD Irish CPO Audit Study 29Sept2017"
          # x  <- "BusUnits/BSC/SSP/Data Science/Statistics/Clinical Trials/Scientific Analytics (Non-study data)/IRISH CPO COPD AUDIT FINDINGS" 
          # drive <- "G"
          # #x  <- "data/R Scripts/PROJECTS/COPD Irish CPO Audit Study 29Sept2017"
          # #drive <- "H"
          # 
          # path <- paste0(x,"/CODE")   # CHANGE TO CODE LATER TO FINALISE
          # path2 <- paste0(x,"/DATA")
          # path3 <- paste0(x,"/OUTPUT")  
          # 
          # work <-    paste0(drive,":/", path)
          # nonwork<- paste0("~/",drive,"/", path)
          # if (where=="home") {wd<- nonwork} else {wd<-work}
          # 
          # work2 <-   paste0(drive,":/", path2)
          # nonwork2 <- paste0("~/",drive,"/", path2)
          # if (where=="home") {wd2<- nonwork2} else {wd2<-work2}
          # 
          # work3 <-   paste0(drive,":/", path3)
          # nonwork3 <- paste0("~/",drive,"/", path3)
          # if (where=="home") {wd3<- nonwork3} else {wd3<-work3}
          # 
          # setwd(wd2)
          # 
          # opts_knit$set(root.dir = wd2)   
          # 
          
 
          p3 <- function(x) {formatC(x, format="f", digits=3)}
          p4 <- function(x) {formatC(x, format="f", digits=4)}
          p2 <- function(x) {formatC(x, format="f", digits=2)}
          p1 <- function(x) {print(formatC(x, format="f", digits=1),quote=FALSE)}
          p1 <- function(x) {formatC(x, format="f", digits=1)}
    
          # Function to construct 95% CI for variance
          var.interval = function(var, d.f, conf.level) {
             chilower = qchisq((1 - conf.level)/2, d.f)
             chiupper = qchisq((1 - conf.level)/2, d.f, lower.tail = FALSE)
             c(d.f * var/chiupper, d.f * var/chilower)
                  }
 
          options(mc.cores=parallel::detectCores())
         
          # function to calculate mode
          Mode <- function(x) {
            ux <- unique(x)
            ux[which.max(tabulate(match(x, ux)))]
          }
          




idx<-
  c("ID	flag trt	time	Censoring
CUKG489B2205AL_23861201	Y	Moxidex	1	0
CUKG489B2205AL_23861202	Y	Tubes	1	0
CUKG489B2205AL_23861205	Y	Tubes	5	0
CUKG489B2205AL_23861206	Y	Moxidex	2	1
CUKG489B2205AL_23861207	Y	Tubes	25	1
CUKG489B2205AL_23861208	Y	Moxidex	1	0
CUKG489B2205AL_23861211	Y	Tubes	4	1
CUKG489B2205AL_23861212	Y	Moxidex	1	0
CUKG489B2205AL_23861213	Y	Tubes	1	0
CUKG489B2205AL_23861214	Y	Moxidex	2	0
CUKG489B2205AL_23861216	Y	Tubes	4	1
CUKG489B2205AL_23861218	Y	Moxidex	1	0
CUKG489B2205AL_23861219	Y	Tubes	1	0
CUKG489B2205AL_23861221	Y	Moxidex	8	0
CUKG489B2205AL_23861222	Y	Moxidex	1	0
CUKG489B2205AL_23861224	Y	Tubes	1	0
CUKG489B2205AL_23861225	Y	Moxidex	15	1
CUKG489B2205AL_23861227	Y	Tubes	1	0
CUKG489B2205AL_23861228	Y	Tubes	1	0
CUKG489B2205AL_23861229	Y	Moxidex	1	0
CUKG489B2205AL_23861231	Y	Moxidex	1	0
CUKG489B2205AL_25061302	Y	Moxidex	1	0
CUKG489B2205AL_25061303	Y	Tubes	1	0
CUKG489B2205AL_25061305	Y	Tubes	1	0
CUKG489B2205AL_25061306	Y	Moxidex	7	1
CUKG489B2205AL_25061307	Y	Moxidex	7	1
CUKG489B2205AL_25061309	Y	Tubes	7	1
CUKG489B2205AL_25061310	Y	Tubes	1	0
CUKG489B2205AL_25061311	Y	Moxidex	1	0
CUKG489B2205AL_25061314	Y	Tubes	3	1
CUKG489B2205AL_25061315	Y	Moxidex	1	0
CUKG489B2205AL_25061316	Y	Moxidex	7	1
CUKG489B2205AL_25061317	Y	Tubes	1	0
CUKG489B2205AL_25061319	Y	Moxidex	1	0
CUKG489B2205AL_25061321	Y	Tubes	1	0
CUKG489B2205AL_25061322	Y	Tubes	5	1
CUKG489B2205AL_25061324	Y	Moxidex	1	0
CUKG489B2205AL_25061325	Y	Tubes	1	0
CUKG489B2205AL_25061327	Y	Moxidex	1	0
CUKG489B2205AL_25061329	Y	Tubes	1	0
CUKG489B2205AL_25061330	Y	Moxidex	8	1
CUKG489B2205AL_25061331	Y	Moxidex	1	0
CUKG489B2205AL_25061333	Y	Tubes	1	0
CUKG489B2205AL_25061334	Y	Moxidex	3	0
CUKG489B2205AL_25061336	Y	Tubes	2	0
CUKG489B2205AL_28421401	Y	Tubes	1	0
CUKG489B2205AL_28421402	Y	Moxidex	1	0
CUKG489B2205AL_28901501	Y	Tubes	2	1
CUKG489B2205AL_28901503	Y	Moxidex	8	1
CUKG489B2205AL_28901505	Y	Tubes	7	1
CUKG489B2205AL_28901506	Y	Moxidex	1	0
CUKG489B2205AL_28901507	Y	Moxidex	1	0
CUKG489B2205AL_28901509	Y	Tubes	8	1
CUKG489B2205AL_30671701	Y	Tubes	1	0
CUKG489B2205AL_30671703	Y	Moxidex	3	0
CUKG489B2205AL_30671704	Y	Tubes	13	1
CUKG489B2205AL_30671706	Y	Moxidex	1	0
CUKG489B2205AL_30671707	Y	Moxidex	1	0
CUKG489B2205AL_30671708	Y	Tubes	7	1
CUKG489B2205AL_30671711	Y	Tubes	1	0
CUKG489B2205AL_30671712	Y	Moxidex	1	0
CUKG489B2205AL_30671714	Y	Tubes	8	0
CUKG489B2205AL_30671715	Y	Moxidex	7	0
CUKG489B2205AL_30671717	Y	Moxidex	1	0
CUKG489B2205AL_30671718	Y	Tubes	2	1
CUKG489B2205AL_30671720	Y	Tubes	4	0
CUKG489B2205AL_30671721	Y	Moxidex	1	0
CUKG489B2205AL_30671722	Y	Moxidex	7	1
CUKG489B2205AL_30671723	Y	Tubes	2	1
CUKG489B2205AL_30671725	Y	Moxidex	5	0
CUKG489B2205AL_46541901	Y	Tubes	7	1
CUKG489B2205AL_46541902	Y	Moxidex	1	0
CUKG489B2205AL_46541904	Y	Tubes	7	1
CUKG489B2205AL_46541905	Y	Moxidex	1	0
CUKG489B2205AL_46541908	Y	Moxidex	1	0
CUKG489B2205AL_46541909	Y	Tubes	16	1
CUKG489B2205AL_46541911	Y	Tubes	3	0
CUKG489B2205AL_46541912	Y	Moxidex	1	0
CUKG489B2205AL_46541914	Y	Tubes	14	1
CUKG489B2205AL_46541915	Y	Moxidex	3	0
CUKG489B2205AL_46541917	Y	Tubes	4	0
CUKG489B2205AL_46541918	Y	Moxidex	1	0
CUKG489B2205AL_46541919	Y	Moxidex	1	0
CUKG489B2205AL_46541920	Y	Tubes	9	1
CUKG489B2205AL_46541922	Y	Moxidex	1	0
CUKG489B2205AL_46541923	Y	Tubes	9	1
CUKG489B2205AL_46541925	Y	Tubes	1	1
CUKG489B2205AL_46592301	Y	Tubes	4	1
CUKG489B2205AL_46592303	Y	Moxidex	1	0
CUKG489B2205AL_46592305	Y	Tubes	1	0
CUKG489B2205AL_46592306	Y	Moxidex	1	0
CUKG489B2205AL_46592307	Y	Moxidex	1	0
CUKG489B2205AL_46592309	Y	Tubes	7	1
CUKG489B2205AL_46592310	Y	Moxidex	1	0
CUKG489B2205AL_46592311	Y	Tubes	1	0
CUKG489B2205AL_46592313	Y	Tubes	1	0
CUKG489B2205AL_46592315	Y	Moxidex	1	0
CUKG489B2205AL_46592317	Y	Moxidex	1	0
CUKG489B2205AL_46592318	Y	Tubes	1	0
CUKG489B2205AL_46592319	Y	Tubes	1	0
CUKG489B2205AL_46592321	Y	Moxidex	1	0
CUKG489B2205AL_46592322	Y	Tubes	1	0
CUKG489B2205AL_46672002	Y	Tubes	7	1
CUKG489B2205AL_46672003	Y	Moxidex	9	1
CUKG489B2205AL_46672004	Y	Moxidex	2	0
CUKG489B2205AL_46672005	Y	Tubes	6	1
CUKG489B2205AL_46672008	Y	Moxidex	1	0
CUKG489B2205AL_46672009	Y	Tubes	1	0
CUKG489B2205AL_46672011	Y	Tubes	6	1
CUKG489B2205AL_46672012	Y	Moxidex	11	0
CUKG489B2205AL_46672014	Y	Tubes	2	1
CUKG489B2205AL_47382101	Y	Tubes	2	0
CUKG489B2205AL_47382102	Y	Moxidex	2	0
CUKG489B2205AL_47382104	Y	Moxidex	1	0
CUKG489B2205AL_47382106	Y	Tubes	1	0
CUKG489B2205AL_47382107	Y	Moxidex	8	0
CUKG489B2205AL_47382108	Y	Tubes	1	0
CUKG489B2205AL_47382110	Y	Tubes	1	0
CUKG489B2205AL_47382111	Y	Moxidex	1	0
CUKG489B2205AL_49942501	Y	Moxidex	2	0
CUKG489B2205AL_49942502	Y	Tubes	2	0
CUKG489B2205AL_50182601	Y	Tubes	7	1
CUKG489B2205AL_50182602	Y	Moxidex	14	1
CUKG489B2205AL_50182604	Y	Moxidex	1	0
CUKG489B2205AL_50182605	Y	Tubes	1	0
CUKG489B2205AL_50182607	Y	Moxidex	1	0
CUKG489B2205AL_50182609	Y	Tubes	7	1
CUKG489B2205AL_50182610	Y	Tubes	3	0
CUKG489B2205AL_50182611	Y	Moxidex	1	0
CUKG489B2205AL_50182613	Y	Tubes	1	0
CUKG489B2205AL_50182615	Y	Moxidex	1	0
CUKG489B2205AL_50182617	Y	Tubes	2	0
CUKG489B2205AL_50182618	Y	Moxidex	2	0
CUKG489B2205AL_50182619	Y	Tubes	14	1
CUKG489B2205AL_50182621	Y	Moxidex	8	0
CUKG489B2205AL_50182622	Y	Tubes	7	1
CUKG489B2205AL_50182623	Y	Moxidex	1	0
CUKG489B2205AL_50182625	Y	Moxidex	1	0
CUKG489B2205AL_50182627	Y	Tubes	4	1
CUKG489B2205AL_50182629	Y	Moxidex	7	0
CUKG489B2205AL_50182630	Y	Tubes	1	0
CUKG489B2205AL_50182631	Y	Moxidex	1	0
CUKG489B2205AL_50182633	Y	Tubes	7	1
CUKG489B2205AL_50182634	Y	Tubes	1	0
CUKG489B2205AL_50182635	Y	Moxidex	3	0
CUKG489B2205AL_50233101	Y	Moxidex	7	0
CUKG489B2205AL_50233102	Y	Tubes	8	1
CUKG489B2205AL_50253201	Y	Tubes	7	1
CUKG489B2205AL_50253203	Y	Moxidex	13	0
CUKG489B2205AL_50253205	Y	Tubes	2	0
CUKG489B2205AL_50402902	Y	Tubes	3	1
CUKG489B2205AL_50402903	Y	Moxidex	1	0
CUKG489B2205AL_50402904	Y	Tubes	1	0
CUKG489B2205AL_50402906	Y	Moxidex	7	1
CUKG489B2205AL_50402907	Y	Moxidex	3	0
CUKG489B2205AL_50402908	Y	Tubes	1	0
CUKG489B2205AL_50402910	Y	Tubes	16	1
CUKG489B2205AL_50402912	Y	Moxidex	1	0
CUKG489B2205AL_50402913	Y	Tubes	13	1
CUKG489B2205AL_50402915	Y	Moxidex	1	0
CUKG489B2205AL_50402916	Y	Moxidex	3	0
CUKG489B2205AL_50402918	Y	Tubes	1	0
CUKG489B2205AL_50402920	Y	Tubes	1	0
CUKG489B2205AL_50402921	Y	Moxidex	1	0
CUKG489B2205AL_50402922	Y	Tubes	6	1
CUKG489B2205AL_51232802	Y	Tubes	5	0
CUKG489B2205AL_51232803	Y	Moxidex	2	0
CUKG489B2205AL_51232804	Y	Tubes	1	0
CUKG489B2205AL_51232806	Y	Moxidex	5	0
CUKG489B2205AL_51232808	Y	Tubes	7	1
CUKG489B2205AL_51232809	Y	Moxidex	14	1
CUKG489B2205AL_51232810	Y	Tubes	1	0
CUKG489B2205AL_51232811	Y	Moxidex	1	0
CUKG489B2205AL_51232813	Y	Tubes	1	0
CUKG489B2205AL_51232815	Y	Moxidex	1	0
CUKG489B2205AL_51232816	Y	Moxidex	8	0
CUKG489B2205AL_51232818	Y	Tubes	4	1
CUKG489B2205AL_51232820	Y	Tubes	1	0
CUKG489B2205AL_51232821	Y	Moxidex	1	0
CUKG489B2205AL_51782201	Y	Tubes	11	0
CUKG489B2205AL_51782202	Y	Moxidex	7	1
CUKG489B2205AL_51782204	Y	Tubes	1	0
CUKG489B2205AL_51782205	Y	Moxidex	5	0
CUKG489B2205AL_51782207	Y	Moxidex	1	0
CUKG489B2205AL_51782209	Y	Tubes	8	1
CUKG489B2205AL_51782210	Y	Tubes	1	0
CUKG489B2205AL_51782211	Y	Moxidex	1	0
CUKG489B2205AL_51782213	Y	Moxidex	2	0
CUKG489B2205AL_51782215	Y	Tubes	1	0
CUKG489B2205AL_51782216	Y	Moxidex	4	0
CUKG489B2205AL_51782217	Y	Tubes	7	1
CUKG489B2205AL_51782220	Y	Moxidex	3	0
CUKG489B2205AL_51782221	Y	Tubes	7	1
CUKG489B2205AL_51782223	Y	Tubes	1	0
CUKG489B2205AL_51782224	Y	Moxidex	3	0
CUKG489B2205AL_51782225	Y	Tubes	4	0
CUKG489B2205AL_51782227	Y	Moxidex	3	0
CUKG489B2205AL_51782228	Y	Tubes	1	0
CUKG489B2205AL_51782229	Y	Moxidex	12	0
CUKG489B2205AL_51782231	Y	Tubes	7	1
CUKG489B2205AL_51782232	Y	Moxidex	3	0
CUKG489B2205AL_51782234	Y	Tubes	7	1
CUKG489B2205AL_51782235	Y	Moxidex	1	0")

Thinh<- read.table(textConnection(idx), header=TRUE)

d <- Thinh

d<-plyr::arrange(d, trt,time )



require(rms)
dd <- datadist(d)
options(datadist='dd')

addmargins(table(d$trt, d$Censoring)) # so 0 are events, censoring =1!  + indicates censoring


cat("\nConvert so that 1 indicates event and 0 censoring\n")
d$Censoring <- ifelse(d$Censoring==0,1,0)
addmargins(table(d$trt, d$Censoring)) # so 0 are events, censoring =1!  + indicates censoring

S <- Surv(d$time ,d$Censoring)
f <- cph(S ~  trt, x=TRUE, y=TRUE, data=d)
f
summary(f)

f2 <- npsurv(S ~ trt, data=d)

x <- c("#e41a1c","#377eb8")

survplot(f2,  n.risk=TRUE, levels.only=T, conf.int=T,
         aehaz=TRUE,
         conf=c("bands"), 
         col.fill=gray(seq(.95, .75, length=4)),
         #col.fill= c(rgb(0,0,0,0.1)), 
         type=c("kaplan-meier"),
         lty=1,  col=x, xlab="Days", abbrev.label=T, 
         label.curves = list(keys = "lines"), #bty='n',
         y.n.risk= 0, cex.n.risk=.6, time.inc=2)
#lines(f2,  col= x) 



fit <- survfit(Surv(time ,Censoring ) ~ trt, d)
fit
#plot(survfit(Surv(time ,Censoring ) ~ trt, d), lty = 2:3) 



quantile(fit)


summary(fit)



cat("\nBootstrap Kaplan-Meier Estimates\n")

Tu<-bootkm(S[d$trt=='Tubes'], q=0.5, B=2000,   pr=FALSE)
Mo<-bootkm(S[d$trt!='Tubes'], q=0.5, B=2000,   pr=FALSE)

quantile(Tu, c(.025,.975), na.rm=TRUE)
quantile(Mo, c(.025,.975), na.rm=TRUE)

 # x<-S[d$trt=='Tubes']
 # x<- x[!grepl("\\+", x)]
 # x<-as.vector(x)
 # median(x)  
 # sort(x)[qbinom(c(0.025, 0.975), size=length(x), prob=0.5)]            
 # 
 # x<-S[d$trt!='Tubes']
 # x<- x[!grepl("\\+", x)]
 # x<-as.vector(x)
 # median(x)  
 # sort(x)[qbinom(c(0.025, 0.975), size=length(x), prob=0.5)]  

 
     # options(width=70)
     # opts_knit$set(root.dir = wd)   # THIS SETS YOUR WORKING DIRECTORY
     # sessionInfo()
     # print(wd)
 
 
stopTime<-proc.time()
 

      # move stangle R file to a folder in GPS
      # put this at bottom & give it the same name as RMD file , replace any blanks with underscore
      # https://amywhiteheadresearch.wordpress.com/2014/11/12/copying-files-with-r/
     # xxx <- "001 COPD STUDY MASTER.R"
     # rcode <-  gsub(' ','_', trimws(xxx))          # replace blank with underscore, this is needed
     # file.copy(rcode, wd,  overwrite=TRUE)         # make a copy of the r code in a folder of choice












