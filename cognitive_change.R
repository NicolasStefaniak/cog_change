#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button in R Studio
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(shinythemes)
library(huxtable)
library(ggplot2)
library(ggpubr)
library(jtools)
library(ciTools)
library(dplyr)

# function from ciTools that have been corrected for avoiding erro --------

# Copyright (C) 2017 Institute for Defense Analyses
#
# This file is part of ciTools.
#
# ciTools is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ciTools is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ciTools. If not, see <http://www.gnu.org/licenses/>.

#' Prediction Intervals for Generalized Linear Models
#'
#' This function is one of the methods for \code{add_pi}, and is
#' called automatically when \code{add_pi} is used on a \code{fit} of
#' class \code{glm}.
#'
#' Prediction intervals are generated through simulation with the aid
#' \code{arm::sim}, which simulates the uncertainty in the regression
#' coefficients. At the moment, only prediction intervals for Poisson,
#' Quasipoisson, Gaussian, and Gamma GLMs are supported. Note that if
#' the response is count data, prediction intervals are only
#' approximate. Simulation from the QuasiPoisson model is performed
#' with the negative binomial distribution, see Gelman and Hill
#' (2007).
#'
#' @param df A data frame of new data.
#' @param fit An object of class \code{glm}.
#' @param alpha A real number between 0 and 1. Controls the confidence
#'     level of the interval estimates.
#' @param names \code{NULL} or character vector of length two. If
#'     \code{NULL}, prediction bounds automatically will be named by
#'     \code{add_pi}, otherwise, the lower prediction bound will be
#'     named \code{names[1]} and the upper prediction bound will be
#'     named \code{names[2]}.
#' @param yhatName A string. Name of the predictions vector.
#' @param nSims A positive integer. Determines the number of
#'     simulations to run.
#' @param ... Additional arguments.
#'
#' @return A dataframe, \code{df}, with predicted values, upper and lower
#'     prediction bounds attached.
#'
#' @seealso \code{\link{add_ci.glm}} for confidence intervals for
#'     \code{glm} objects, \code{\link{add_probs.glm}} for conditional
#'     probabilities of \code{glm} objects, and
#'     \code{\link{add_quantile.glm}} for response quantiles of
#'     \code{glm} objects.
#'
#' @examples
#' # Fit a Poisson model
#' fit <- glm(dist ~ speed, data = cars, family = "poisson")
#' # Add prediction intervals and fitted values to the original data frame
#' add_pi(cars, fit)
#' # Try a different confidence level
#' add_pi(cars, fit, alpha = 0.5)
#' # Try custom names for the prediction bounds (may be useful for plotting)
#' add_pi(cars, fit, alpha = 0.5, names = c("lwr", "upr"))
#'
#' @export


add_pi.glm <- function(df, fit, alpha = 0.05, names = NULL, yhatName = "pred",
                       nSims = 2000,  ...){
  
  
  if (is.null(names)) {
    names[1] <- paste("LPB", alpha/2, sep = "")
    names[2] <- paste("UPB", 1 - alpha/2, sep = "")
  }
  if ((names[1] %in% colnames(df)))
    warning ("These PIs may have already been appended to your dataframe. Overwriting.")
  
  
  if(fit$family$family == "binomial")
    if(max(fit$prior.weights) == 1)
    {stop("Prediction intervals for Bernoulli response variables aren't useful")}else {
      warning("Treating weights as indicating the number of trials for a binomial regression where the response is the proportion of successes")
      warning("The response variable is not continuous so Prediction Intervals are approximate")
    }
  
  if(fit$family$family %in% c("poisson", "quasipoisson"))
    warning("The response is not continuous, so Prediction Intervals are approximate")
  
  if(!(fit$family$family %in% c("poisson", "quasipoisson", "Gamma", "binomial", "gaussian")))
    stop("Unsupported family")
  
  if(fit$family$family == "gaussian")
    pi_gaussian(df, fit, alpha, names, yhatName)  else   sim_pi_glm(df, fit, alpha, names, yhatName, nSims)
}

pi_gaussian <- function(df, fit, alpha, names, yhatName){
  sigma_sq <- summary(fit)$dispersion
  inverselink <- fit$family$linkinv
  out <- predict(fit, newdata = df, se.fit = TRUE)
  se_terms <- out$se.fit
  t_quant <- qt(p = alpha/2, df = fit$df.residual, lower.tail = FALSE)
  se_global <- sqrt(sigma_sq + se_terms^2)
  lwr <- inverselink(out$fit) - t_quant * se_global
  upr <- inverselink(out$fit) + t_quant * se_global
  
  if(is.null(df[[yhatName]]))
    df[[yhatName]] <- inverselink(out$fit)
  df[[names[1]]] <- lwr
  df[[names[2]]] <- upr
  data.frame(df)
}

get_sim_response <- function(df, fit, nSims){
  
  nPreds <- NROW(df)
  modmat <- model.matrix(fit, data = df)
  response_distr <- fit$family$family
  inverselink <- fit$family$linkinv
  overdisp <- summary(fit)$dispersion
  sims <- arm::sim(fit, n.sims = nSims)
  sim_response <- matrix(NA, ncol = nSims, nrow = nPreds)
  
  for (i in 1:nSims){
    yhat <- inverselink(modmat %*% sims@coef[i,])
    if(response_distr == "poisson"){
      sim_response[,i] <- rpois(n = nPreds,
                                lambda = yhat)
    }
    if(response_distr == "quasipoisson"){
      a <- yhat / (overdisp - 1)
      sim_response[,i] <- rnegbin(n = nPreds,
                                  mu = yhat,
                                  theta = a)
    }
    if(response_distr == "Gamma"){
      sim_response[,i] <- rgamma(n = nPreds,
                                 shape = 1/overdisp,
                                 rate = 1/(yhat *overdisp))
    }
    if(response_distr == "binomial"){
      yhat <- inverselink(modmat %*% sims@coef[i,]) * fit$prior.weights
      sim_response[,i] <- rbinom(n = nPreds,
                                 size = fit$prior.weights,
                                 prob = yhat / fit$prior.weights)
    }
    if(response_distr == "gaussian"){
      yhat <- inverselink(modmat %*% sims@coef[i,])
      sim_response[,i] <- rnorm(n = nPreds,
                                mean = yhat,
                                sd = sqrt(overdisp))
    }
    
  }
  sim_response
}

sim_pi_glm <- function(df, fit, alpha, names, yhatName, nSims){
  out <- predict(fit, newdata = df, type = "response")
  sim_response <- get_sim_response(df = df, fit = fit, nSims = nSims)
  lwr <- apply(sim_response, 1, FUN = quantile, probs = alpha/2, type = 1, na.rm=T)
  upr <- apply(sim_response, 1, FUN = quantile, probs = 1 - alpha / 2, type = 1, na.rm=T)
  
  if(fit$family$family == "binomial"){
    out <- out * fit$prior.weights
    warning("For binomial models, add_pi's column of fitted values reflect E(Y|X) rather than typical default for logistic regression, pHat")
  }
  
  if(is.null(df[[yhatName]]))
    df[[yhatName]] <- out
  df[[names[1]]] <- lwr
  df[[names[2]]] <- upr
  data.frame(df)
}
# Define UI for application that draws a histogram
ui <- navbarPage("Cognitive change assessment",
                 tabPanel("How to use this app?",
                          shinyUI(
                            fluidPage(
                          h4("Welcome to this app" ),
                          p("This app provides the supplemental material to XXX et al. (submitted)'s paper entitled 
                            ?Assessing cognitive changes in multiple sclerosis: Criteria for a reliable decision?", 
                            style = "color:white"),
                          p("This app allows users to obtain standard scores and estimate cognitive change from 
                            seven tests taken from the French version of the BRB-N (Rao et al., 1991), 
                            known as the BCcogSEP (battery evaluating cognitive functions in multiple sclerosis; 
                            for more information, see Dujardin, Sockeel, Cabaret, De Seze, & Vermersch, 2004).", style = "color:white"),
                          p("To use this app, you need to click on the tab 'Assessing cognitive change' and enter", 
                            strong("raw data"), "for SRT immediate recall, SRT delayed recall, and the digit symbol, 
                            forward/backward digit span, phonemic fluency, and category fluency subtests.", style = "color:white"),
                          p(strong("Table 1"), "provides the standard scores for each task at T0 and T1, and indices of
                            cognitive change yielded by five different methods." ,style = "color:white"),
                          p(strong("SD method"), "is the difference between T0 and T1, divided by the SD of the control group.",
                            style = "color:white"),
                          p( withMathJax("$$SD_{method}=\\frac{(T1-T0)}{sd_{control}}$$")),
                          p(strong("WARNING: this method should not be used and is  presented for comparison purpose only
                                  because it represents a classical index used in clinical practice.")),
                         
                          p(strong("RCI",tags$sub("Chelune")), "is the difference between T0 and T1 divided by the standard error of difference 
                            (Jacobson and Truax, 1991) with a correction for practice effect (Chelune et al., 1993).",
                            style = "color:white"),
                          p( withMathJax("$$RCI_{Chelune} =\\frac{(T1-T0)-(M_1-M_0)}{sed}$$")),
                          p(HTML(paste0("where M",tags$sub("0"), " and M", tags$sub("1") ," are the control group's means at T0 and T1 "))), 
                          p("where sed is ",withMathJax("$$sed=\\sqrt{2*(sd_{T0} \\sqrt{1-r_{xy}})^2 }$$")),
                          p(HTML(paste0("and r",tags$sub("xy")  ," is the test-retest coefficient"))), 
                          
                          p(strong("RCI",tags$sub("Iverson")), "does not assume that the SEM remains constant at each testing time 
                            (Iverson, 2001).",style = "color:white"),
                          p("thus sed is obtained by",withMathJax("$$sed=\\sqrt{(sd_{T0} \\sqrt{1-r_{xy}})^2 + (sd_{T1} \\sqrt{1-r_{xy}})^2}$$")),
                          p(strong("SRB",tags$sub("(1 predictor)")), "is the standardized-regression-based method in which the value at T1 is predicted by
                            the value at T0  (McSweeny et al., 1993).",style = "color:white"),
                          p( withMathJax("$$SRB_{(1 predictor)}=\\frac{(T1-predicted_{T1})}{SEE}$$")),
                          p("where SEE is ",withMathJax("$$SEE=sd_{T0}\\sqrt{1-r^2}$$")),
                          p(HTML(paste0("r", tags$sup("2"), " is the coefficient of determination of the regression model"))), 
                          p(strong("SRB",tags$sub("(4 predictors)") ), "is based on the SRB",tags$sub("(1 predictor)"), "but predictions take into account age, sex, 
                            education level (for the formula of the standard error of estimates for a new observation, see Crawford et al., 2012)",style = "color:white"),
                          p(strong("GSRB",tags$sub("(4 predictors)") ), "is based on the generalized linead model to predict the observed value. The confidence interval 
                            is base on simulation, which could explain small discrepancies between two sessions with the same values.",style = "color:white")
                            
                          )
                          )
                          ),
                 tabPanel("Assessing cognitive change", 
                          shinyUI(
                              fluidPage(theme = shinytheme("superhero"), 
                                    fluidRow(
                                      column(4,h3("Demographic information"),
                                             numericInput("Age", label="Age", value=30 ,min =1, max=100, step=1, width='250px'),
                                             numericInput("Educ", label="Education level", value=12 ,min =3, max=15, step=1, width='250px'),
                                             selectInput("Sex", "Patient's sex", c("M", "F","other")),
                                             selectInput("ES", "Model used to compute effect size", c("linear model", "generalized linear model"))
                                      ),
                                      
                                      column(4,h3("Time 0"),
                                             numericInput("SRT.av", label="SRT average recall", value=12 ,min =0, max=14, step=0.01, width='100px'),
                                             numericInput("SRT.d", label="SRT delayed recall", value=14 ,min =0, max=14, step=1, width='100px'),
                                             numericInput("DS", label="Digit symbol", value=59 ,min =0, max=100, step=1, width='100px'),
                                             numericInput("DDS", label="Forward digit span", value=2 ,min =0, max=16, step=1, width='100px'),
                                             numericInput("BDS", label="Backward digit span", value=4 ,min =0, max=16, step=1, width='100px'),
                                             numericInput("PF", label="Phonemic fluency", value=18 ,min =0, max=100, step=1, width='100px'),
                                             numericInput("SF", label="Semantic fluency", value=21 ,min =0, max=100, step=1, width='100px')
                                                            ),
                                     column(4,h3("Time 1"),
                                            numericInput("SRT.av.T1", label="SRT average recall", value=12 ,min =0, max=14, step=0.01, width='100px'),
                                            numericInput("SRT.d.T1", label=" SRT delayed recall", value=14 ,min =0, max=14, step=1, width='100px'),
                                            numericInput("DS.T1", label="Digit symbol", value=59 ,min =0, max=100, step=1, width='100px'),
                                            numericInput("DDS.T1", label="Forward digit span", value=2 ,min =0, max=16, step=1, width='100px'),
                                            numericInput("BDS.T1", label="Backward digit span", value=4 ,min =0, max=16, step=1, width='100px'),
                                            numericInput("PF.T1", label="Phonemic fluency", value=18 ,min =0, max=100, step=1, width='100px'),
                                            numericInput("SF.T1", label="Semantic fluency", value=21 ,min =0, max=100, step=1, width='100px')
                                                            )),
                                    fluidRow(
                                    column(10,h4("Table 1. Scaled scores at T0 and T1, and 7 indexes of cognitive change between T0 and T1"),
                                           DTOutput("table0"),
                                           p(strong("* cognitive worsening")),
                                           p(strong("** absence of cognitive change")),
                                           p(strong("*** cognitive improvement"))
                                    )),
                                  fluidRow(
                                    column(12,h4("Figure 1. Patient's evolution in comparison to the reference group and effect size of observation (zop)"),
                                    plotOutput("Figure1", width = "80%", height="1600px"))
                                  )
                          )
                          )))
                          





# Define server logic required to draw a histogram
server <- function(input, output) {
    
    
    output$table0 <-renderDT ({
             # name of test
        label.tasks<-c("SRT average recall","SRT delayed recall", "Digit symbol", 
                 "Forward digit span", "Backward digit span", "Phonemic fluency", "Semantic fluency")
  
        # scale values at T0 and T1 
        mean.control<-c(11.6277528089888, 13.685393258427, 57.8089887640449, 7.43820224719101, 
          6.93258426966292, 16.9213483146067, 22.6404494382022)
        sd.control<-c(1.32668824933724, 2.29423935512195, 12.330557390349, 2.33527814276711, 
          2.29524089055132, 4.49804179388887, 5.61905684218534)


        Demo<-list(Age=input$Age, LoE= input$Educ, Sex=ifelse(input$Sex=="M",1, 0))
        T0.subject<-c(input$SRT.av, input$SRT.d,input$DS, input$DDS,input$BDS,
                       input$PF,input$SF)
        T1.subject<-c(input$SRT.av.T1,input$SRT.d.T1, input$DS.T1,input$DDS.T1,input$BDS.T1,
                      input$PF.T1, input$SF.T1)
        
        z0<-(T0.subject-mean.control)/sd.control
        z1<-(T1.subject-mean.control)/sd.control
        
        z0<-round(z0, 3)
        z1<-round(z1, 3 ) 
        
        # SD method
        SD_method<-(T1.subject - T0.subject)/sd.control
        SD_method<-round(SD_method, 3)
        SD_method<-ifelse(SD_method>1, paste0(SD_method, " ***"),
                          ifelse(SD_method<(-1),paste0(SD_method, " *"),paste0(SD_method, " **")) ) 
        
        
        # RCI Chelune
        cor.HC.T0<-c(0.546155948288594, 0.595723868426037, 0.84545495183893, 0.77641980242785, 
          0.780995419172628, 0.660691956080756, 0.692473183704526)
        SEM<-sd.control *(1-cor.HC.T0)^.5
        sed<-(2*SEM^2)^0.5
        
        
        mean.control.T1<- c(11.78, 14, 59.0674157303371, 7.75280898876404, 7.39325842696629, 
          17.4606741573034, 22.8202247191011)
        RCI2<-((T1.subject - T0.subject)-(mean.control.T1-mean.control))/sed
        RCI2<-round(RCI2, 3)  
        RCI2<-ifelse(RCI2>1.645, paste0(RCI2, " ***"),
                    ifelse(RCI2<(-1.645),paste0(RCI2, " *"),paste0(RCI2, " **")) ) 
        
        # RCI Iverson
        
        sd.T1.C<-c(1.21280197739105, 1.77097815807074, 12.3231099904435, 2.1913097934403, 
          2.28442342799352, 4.42116855124628, 5.12028452190805)
        SEM.T1<-sd.T1.C *(1-cor.HC.T0)^.5
        sed3<-(SEM^2+SEM.T1^2)^0.5
        RCI3<-((T1.subject - T0.subject)-(mean.control.T1-mean.control))/sed3
        RCI3<-round(RCI3, 3)  
        RCI3<-ifelse(RCI3>1.645, paste0(RCI3, " ***"),
                     ifelse(RCI3<(-1.645),paste0(RCI3, " *"),paste0(RCI3, " **")) ) 
        
        # SRB
       coef.list<- structure(list(T1.SRT.Rappel.moyen = c(5.97458230291646, 0.499272541517712), 
                                  T1.SRT.rappel.diff = c(7.70672423830779, 0.459853483407724), 
                                  T1.Code = c(10.222039338326, 0.844944314652864), 
                                  T1.Score.Empan.dir = c(2.33367671848661, 0.728554036336393), 
                                  T1.Score.Empan.inv = c(2.00445952496364,0.777314590402328), 
                                  T1.Fluences.phon = c(6.47194244604317, 0.649400479616307), 
                                  T1.Fluences.cat = c(8.53396067711072, 0.631006203343497)), 
                             class = "data.frame", row.names = c(NA, -2L))
       a<-coef.list[1,]
       b<-coef.list[2,]
       predits<-a+b*T0.subject

       see<-sd.T1.C*((1-cor.HC.T0^2)* (88/87))^0.5
       SRB1<-(T1.subject-predits)/see
       critical.tvalue<-qt(0.95, 89-2)
 
       # SRB 1 predictor
       
       see_nf<-function(syx, nrows, x0, meanx, sd2){
           syx*(
               1+
                   (1/nrows)+ 
                   (
                       (x0-meanx)^2/
                           (sd2*(nrows-1))
                   ))^0.5
       }
       
       see_n<-see_nf(syx= see, nrows= 89, x0= T0.subject, meanx=mean.control, sd2=sd.control^2)
       SRB2<-(T1.subject-predits)/see_n
       SRB2<-round(SRB2,3)
       SRB2<-ifelse(SRB2>critical.tvalue, paste0(SRB2, " ***"),
                    ifelse(SRB2<(-critical.tvalue),paste0(SRB2, " *"),paste0(SRB2, " **")) ) 
       
       # SRB - 4 factors
       
       pred.man.SRTac<-c(`(Intercept)` = 6.38052814916193, GenreH = 0.214540161854636, 
                         Age = -0.00641871267052169, N.etudes = -0.0108080561171946, T0.SRT.Rappel.moyen = 0.493778478008299)
       
       
       pred.man.SRTdel<-c(`(Intercept)` = 9.48187113082139, GenreH = 0.0260622158545004, 
                          Age = -0.015671372955047, N.etudes = -0.077727991935889, T0.SRT.rappel.diff = 0.453590252084855)
       
       pred.man.DS<-c(`(Intercept)` = 10.5752518324084, GenreH = 2.58710105706402, 
                      Age = -0.0414519341405448, N.etudes = 0.0347739333178569, T0.Code = 0.845742505576553)
       
       pred.man.direct<-c(`(Intercept)` = 3.38375064105519, GenreH = 0.098855104611957, 
                          Age = -0.000466621888344133, N.etudes = -0.0933661910317817, 
                          T0.Score.Empan.dir = 0.748939196449153)
       
       pred.man.backward<-c(`(Intercept)` = 2.57067714127993, GenreH = 0.917929236860799, 
                            Age = -0.0172850846067988, N.etudes = 0.00211359423463607, T0.Score.Empan.inv = 0.751572853984314)
       
       pred.man.PF<-c(`(Intercept)` = 1.45610768667931, GenreH = 1.2068308121188, 
                      Age = 0.020817348148685, N.etudes = 0.285298577754853, T0.Fluences.phon = 0.645650098715606)
       
       pred.man.CF<-c(`(Intercept)` = 6.47747318037186, GenreH = -0.71136911342446, 
                      Age = 0.0508947252981979, N.etudes = -0.0545707088473257, T0.Fluences.cat = 0.666652746663211
       )
       NewO<-data.frame( "Genre"=Demo$Sex,
                         "Age"=Demo$Age,
                         "N.etudes"=Demo$LoE,
                         "T0.SRT.Rappel.moyen"=T0.subject[1],
                         "T0.SRT.rappel.diff"=T0.subject[2],
                         "T0.Code"=T0.subject[3],
                         "T0.Score.Empan.dir"=T0.subject[4],
                         "T0.Score.Empan.inv"=T0.subject[5],
                         "T0.Fluences.phon"=T0.subject[6],
                         "T0.Fluences.cat"=T0.subject[7],
                         "T1.SRT.Rappel.moyen"=T1.subject[1],
                         "T1.SRT.rappel.diff"=T1.subject[2],
                         "T1.Code"=T1.subject[3],
                         "T1.Score.Empan.dir"=T1.subject[4],
                         "T1.Score.Empan.inv"=T1.subject[5],
                         "T1.Fluences.phon"=T1.subject[6],
                         "T1.Fluences.cat"=T1.subject[7])
       
      
       
       list.R<-list(structure(c(1, -0.0567158929167967, -0.0704207368324297, 
                                -0.237094689891434, -0.0567158929167967, 1, -0.333158192156128, 
                                -0.471940351983711, -0.0704207368324297, -0.333158192156128, 
                                1, 0.284704637700897, -0.237094689891434, -0.471940351983711, 
                                0.284704637700897, 1), .Dim = c(4L, 4L), 
                              .Dimnames = list(c("Genre", "Age", "N.etudes", "T0.SRT.Rappel.moyen"), 
                                               c("Genre", "Age", "N.etudes", "T0.SRT.Rappel.moyen"))), 
                    structure(c(1, -0.0567158929167967,-0.0704207368324297, -0.261221427747892, 
                                -0.0567158929167967, 1, -0.333158192156128, -0.373855418594585, 
                                -0.0704207368324297, -0.333158192156128, 1, 0.305363148456841, 
                                -0.261221427747892, -0.373855418594585, 0.305363148456841, 1), .Dim = c(4L, 4L), .Dimnames = list(
                                  c("Genre", "Age", "N.etudes", "T0.SRT.rappel.diff"),
                                  c("Genre", "Age", "N.etudes", "T0.SRT.rappel.diff"))), 
                    structure(c(1,-0.0567158929167967, -0.0704207368324297, -0.242269871154987, 
                                -0.0567158929167967, 1, -0.333158192156128, -0.497814603832284, 
                                -0.0704207368324297, -0.333158192156128, 1, 0.283885524750917, 
                                -0.242269871154987, -0.497814603832284, 0.283885524750917, 1), .Dim = c(4L, 
                                                                                                        4L), .Dimnames = list(c("Genre", "Age", "N.etudes", "T0.Code"
                                                                                                        ), c("Genre", "Age", "N.etudes", "T0.Code"))), 
                    structure(c(1,  -0.0567158929167967, -0.0704207368324297, 0.125614061292371, 
                                -0.0567158929167967, 1, -0.333158192156128, -0.162863884660968, 
                                -0.0704207368324297, -0.333158192156128, 1, 0.235378782373347, 
                                0.125614061292371, -0.162863884660968, 0.235378782373347, 1), .Dim = c(4L, 
                                                                                                       4L), .Dimnames = list(c("Genre", "Age", "N.etudes", "T0.Score.Empan.dir"
                                                                                                       ), c("Genre", "Age", "N.etudes", "T0.Score.Empan.dir"))), 
                    structure(c(1, -0.0567158929167967, -0.0704207368324297, 0.0736368797508948, 
                                -0.0567158929167967, 1, -0.333158192156128, -0.109961250245402, 
                                -0.0704207368324297, -0.333158192156128, 1, 0.261437019695802, 
                                0.0736368797508948, -0.109961250245402, 0.261437019695802, 1), .Dim = c(4L, 
                                                                                                        4L), .Dimnames = list(c("Genre", "Age", "N.etudes", "T0.Score.Empan.inv"
                                                                                                        ), c("Genre", "Age", "N.etudes", "T0.Score.Empan.inv"))), 
                    structure(c(1, -0.0567158929167967, -0.0704207368324297, -0.116523721015182, 
                                -0.0567158929167967, 1, -0.333158192156128, 0.0631267984517857, 
                                -0.0704207368324297, -0.333158192156128, 1, 0.0958487510607872, 
                                -0.116523721015182, 0.0631267984517857, 0.0958487510607872, 1
                    ), .Dim = c(4L, 4L), .Dimnames = list(c("Genre", "Age", "N.etudes", 
                                                            "T0.Fluences.phon"), c("Genre", "Age", "N.etudes", "T0.Fluences.phon"))), 
                    structure(c(1, -0.0567158929167967, -0.0704207368324297,  0.0244185434351168, 
                                -0.0567158929167967, 1, -0.333158192156128,   -0.220033024907953,
                                -0.0704207368324297, -0.333158192156128,  1, 0.333899168617547, 
                                0.0244185434351168, -0.220033024907953, 0.333899168617547, 1), .Dim = c(4L, 4L), 
                              .Dimnames = list(c("Genre", "Age", "N.etudes", "T0.Fluences.cat"), 
                                               c("Genre", "Age", "N.etudes", "T0.Fluences.cat"))))
       
       pred.values<-list(pred.man.SRTac, 
                         pred.man.SRTdel,
                         pred.man.DS,
                         pred.man.direct,
                         pred.man.backward,
                         pred.man.PF,
                         pred.man.CF)
       # Mean and SD for demographic variables     
       M<-c(Genre = 0.37078651685393, Age = 43.7191011235955, N.etudes = 13.0449438202247)
       ET<-c(Genre = 0.485752052162186, Age = 13.0879791863694, N.etudes = 2.48591538088063)
       rsq<-c(0.310843341680445, 0.371973676629549, 0.727794681552961, 0.614199604579674, 
             0.65947559562284, 0.473101976369185, 0.503209782036205)
       
       pred.values.man<-c()
       Sn1_tot<-c()
       
       
       for(i in 1:7){
         # solve R
         Rinv<-solve( list.R[[i]]) 
         rii<-diag(Rinv)
         rij<-Rinv[lower.tri(Rinv, diag = FALSE)]
         
         # center and scale 
         M1<-c(M, mean.control[i])
         ET1<-c(ET,  sd.control[i])
         zi<-mapply(function(x, M, sd){(x-M)/sd} , NewO[,c(1:3,(i+3))],M=M1,sd=ET1)
         pred.values.man<-c( pred.values.man, (pred.values[[i]][1] + sum(pred.values[[i]][2:5]*NewO[,c(1:3,(i+3))])))
         # multiply z
         zj<-c()
         for(j in 1:4){
           for(k in (j+1):5){
             zj<-c(zj, zi[j]*zi[k])
           }
         }
         zj<-zj[which(!is.na(zj))]
         
         
         
         
         # Compute Syx (equation 4 from Crawford et al., 2012) 
         sxy<-sqrt(
           ((1-(rsq[i]))* 
              sd.T1.C[i]^2*(89-1))/
             (89-1-4)
         )
         
         # Compute Sn+1 (equation 5 from Crawford et al., 2012) 
         s_n1<-sxy*(1+
                      1/89+ 
                      (1/88* sum(rii*(zi^2))) + 
                      (2/88*sum(rij*zj)))^0.5     
         Sn1_tot<-c(Sn1_tot,s_n1 )
       }
         


       
     t.crit<- 1.663197
     SRB4<-(T1.subject-pred.values.man)/Sn1_tot
     SRB4<-round(SRB4,3)
     SRB4<-ifelse(SRB4>critical.tvalue, paste0(SRB4, " ***"),
                  ifelse(SRB4<(-critical.tvalue),paste0(SRB4, " *"),paste0(SRB4, " **")) ) 
      
    # Generalized SRB
###########################################################     
     data.lm<-structure(list(ID = structure(c(9L, 31L, 18L, 62L, 67L, 22L, 
                                              24L, 35L, 16L, 59L, 78L, 61L, 81L, 88L, 58L, 49L, 82L, 46L, 48L, 
                                              11L, 26L, 7L, 86L, 66L, 54L, 56L, 17L, 51L, 14L, 60L, 50L, 63L, 
                                              15L, 68L, 69L, 30L, 40L, 27L, 13L, 87L, 53L, 25L, 10L, 43L, 39L, 
                                              71L, 89L, 34L, 4L, 5L, 6L, 42L, 72L, 45L, 64L, 8L, 41L, 47L, 
                                              79L, 38L, 1L, 80L, 12L, 44L, 20L, 70L, 23L, 19L, 77L, 76L, 75L, 
                                              33L, 65L, 32L, 3L, 57L, 85L, 74L, 2L, 21L, 73L, 83L, 37L, 55L, 
                                              84L, 36L, 29L, 28L, 52L), 
                                            .Label = c("ARNI54", "ARQU93", "AUMA93", "BAAL93", "BACA65", "BAJE64", "BAJF67", "BANA73", "BEJE68", "BLAU87", 
                                                       "BLBR54", "BOAN88", "BRCL55", "CHFR70", "CHJE86", "CHJO58", "CHMA89", "CHNA61", "CHYA78", "COAU86", "COGO69", "COJE70", "COPA84", "COVI73", 
                                                       "CRPA66", "DAFR67", "DAJE92", "DEAR74", "DELA94", "DESA78", "DOAN92","DUJE90", "DULE93", "FOGA74", "GELA60", "GOJE67", "GOME90", "GRCH54", 
                                                       "GRKA77", "HEVI74", "JAMA56", "JAST91", "LECA94", "LECL90", "LEMA65", "LOAN88", "LOAU93", "LOCH60", "LOPA54", "MACH93\r\n", "MAMY67", 
                                                       "MANA69", "MANO90", "MARI66", "MASA80", "MEIS66", "MODO59", "MOMA77","NACH62", "NIER70", "OLMI60", "ORIS63", "OUSA94", "PEAG88", "PEMA89", 
                                                       "PIJA93", "POFL81", "POGH64", "POMA67", "PRCA53", "PREM88", "PRMI85",  "QUGO93", "RIJA56", "RIPA66", "RISE74", "RISY59", "ROGI62", "SIDA54", 
                                                       "SIFR54", "SISO69", "STKU74", "SULA86", "TAMA61", "THMA61", "THSO74","VICH75", "VISA76", "VOCA81"), class = "factor"), 
                             DDN = structure(c(40L, 80L, 20L, 25L, 63L, 45L, 51L, 17L, 12L, 21L, 22L, 16L, 42L, 59L, 
                                               61L, 5L, 54L, 8L, 13L, 7L, 37L, 35L, 53L, 83L, 33L, 31L, 75L, 
                                               36L, 47L, 46L, 82L, 89L, 70L, 27L, 38L, 62L, 55L, 79L, 10L, 58L, 
                                               77L, 34L, 71L, 87L, 60L, 74L, 64L, 57L, 86L, 30L, 26L, 78L, 69L, 
                                               28L, 72L, 50L, 11L, 85L, 6L, 9L, 24L, 3L, 73L, 68L, 1L, 29L, 
                                               39L, 67L, 48L, 52L, 32L, 4L, 76L, 23L, 2L, 49L, 18L, 66L, 81L, 
                                               44L, 84L, 15L, 65L, 14L, 19L, 43L, 88L, 56L, 41L), 
                                             .Label = c("15/05/1986", "16/05/1993", "17/09/1954", "19/04/1993", "19814", "19885", "19937", 
                                                        "20/07/1988", "20033", "20219", "20786", "21240", "21974", "22/12/1980", 
                                                        "22/12/1986", "22059", "22091", "22352", "22423", "22640", "22724", 
                                                        "22957", "23/05/1990", "23/12/1954", "23107", "23505", "23533", 
                                                        "23896", "23957", "23985", "24137", "24270", "24276", "24303", 
                                                        "24475", "24494", "24496", "24700", "25/05/1984", "25002", "25269", 
                                                        "25286", "25408", "25453", "25833", "25836", "25922", "26/11/1956", 
                                                        "26/11/1959", "26666", "26865", "27063", "27074", "27189", "27264", 
                                                        "27345", "27388", "27558", "27941", "28281", "28318", "28503", 
                                                        "29706", "29857", "30/03/1990", "30/05/1956", "30/11/1978", "31/05/1990", 
                                                        "31367", "31723", "32132", "32288", "32396", "32463", "32616", 
                                                        "32758", "33002", "33564", "33622", "33868", "34014", "34081", 
                                                        "34213", "34250", "34271", "34324", "34368", "34396", "34414"), class = "factor"), 
                             Genre = structure(c(2L, 2L, 1L, 1L, 2L,  2L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 
                                                 2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 
                                                 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 
                                                 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 
                                                 1L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 
                                                 2L, 1L, 2L, 1L), .Label = c("F", "H"), class = "factor"), 
                             Age = c(49, 25, 56, 55, 36, 47, 44, 57, 59, 55, 55, 57, 48, 41, 40, 65, 44, 
                                     30, 58, 64, 51, 52, 45, 25, 52, 53, 29, 52, 48, 48, 25, 25, 32, 
                                     54, 52, 41, 45, 27, 64, 43, 28, 52, 31, 25, 41, 30, 37, 44, 25, 
                                     53, 54, 27, 33, 53, 30, 46, 62, 25, 65, 65, 65, 65, 30, 28, 32, 
                                     53, 35, 40, 59, 45, 52, 25, 29, 28, 25, 59, 58, 62, 25, 49, 25, 
                                     32, 29, 38, 58, 52, 25, 44, 50), 
                             N.etudes = c(15, 10, 12, 12, 15, 10, 10, 14, 15, 15, 15, 10, 14, 12, 15, 10, 13, 12, 8, 16, 
                                          11, 11, 11, 14, 11, 15, 12, 15, 11, 8, 15, 17, 12, 11, 11, 18, 
                                          14, 17, 17, 12, 15, 15, 12, 15, 12, 15, 12, 16, 15, 14, 16, 11, 
                                          17, 15, 15, 12, 14, 14, 14, 9, 10, 9, 14, 13, 14, 8, 12, 11, 
                                          10, 12, 10, 15, 17, 16, 12, 9, 11, 9, 17, 12, 14, 14, 11, 12, 
                                          17, 14, 16, 16, 12), 
                             Date.T0 = structure(c(13L, 12L, 11L, 10L, 
                                                   9L, 8L, 8L, 7L, 7L, 6L, 6L, 5L, 5L, 4L, 4L, 27L, 22L, 27L, 27L, 
                                                   23L, 14L, 24L, 41L, 40L, 45L, 45L, 41L, 45L, 42L, 37L, 38L, 38L, 
                                                   39L, 39L, 39L, 43L, 42L, 42L, 44L, 37L, 36L, 16L, 36L, 47L, 28L, 
                                                   28L, 38L, 27L, 17L, 20L, 21L, 39L, 38L, 15L, 31L, 24L, 40L, 40L, 
                                                   17L, 18L, 18L, 18L, 19L, 20L, 23L, 24L, 24L, 24L, 26L, 25L, 25L, 
                                                   30L, 30L, 32L, 32L, 33L, 33L, 35L, 26L, 42L, 29L, 3L, 34L, 49L, 
                                                   46L, 50L, 2L, 48L, 1L), 
                                                 .Label = c("16/11/2019", "18/10/2019", 
                                                            "23/02/2019", "43081", "43086", "43088", "43089", "43103", "43105", 
                                                            "43117", "43118", "43121", "43125", "43459", "43460", "43461", 
                                                            "43475", "43481", "43482", "43483", "43484", "43486", "43488", 
                                                            "43491", "43492", "43493", "43497", "43498", "43499", "43501", 
                                                            "43504", "43505", "43506", "43507", "43511", "43512", "43513", 
                                                            "43516", "43517", "43518", "43526", "43527", "43532", "43533", 
                                                            "43534", "43535", "43536", "43657", "43672", "43705"), class = "factor"), 
                             T0.SRT.Rappel.moyen = c(12.73, 12.6, 12, 10.5, 13.25, 9.3, 
                                                     12.5, 10, 11.73, 9.7, 11, 7.4, 10.3, 10.9, 12.09, 7.09, 12.5, 
                                                     11.18, 11.87, 11.27, 9.64, 11.82, 12.11, 10.91, 11.09, 12.25, 
                                                     12.5, 12.33, 12.67, 12.22, 12.7, 12.7, 12.44, 12.45, 11.9, 
                                                     12.8, 13, 13.17, 11.9, 11.44, 12.86, 11.4, 11.86, 13.5, 12.13, 
                                                     12.78, 12.75, 12.56, 11.18, 12, 13.4, 11.71, 11.75, 10.55, 
                                                     11.27, 12.83, 8.54, 11.5, 7.55, 11.09, 9.91, 10.91, 11.71, 
                                                     11.44, 12.33, 12.27, 12.75, 8.72, 11.5, 11.88, 12.19, 12.37, 
                                                     12.83, 11.86, 12.38, 9.18, 11.64, 10.35, 12, 12.57, 11.6, 
                                                     12.66, 12.17, 12.67, 12.38, 10, 13, 11.73, 12.71), 
                             T0.SRT.rappel.diff = c(15, 
                                                    15, 15, 14, 15, 8, 14, 13, 15, 11, 14, 5, 12, 13, 12, 4, 
                                                    15, 14, 13, 15, 12, 15, 15, 13, 15, 14, 15, 15, 15, 15, 15, 
                                                    15, 15, 15, 13, 15, 15, 15, 15, 14, 15, 15, 14, 15, 15, 15, 
                                                    15, 15, 14, 15, 15, 15, 14, 12, 15, 15, 7, 15, 10, 15, 11, 
                                                    13, 14, 14, 15, 13, 14, 9, 14, 14, 12, 15, 15, 15, 12, 7, 
                                                    13, 15, 14, 15, 15, 15, 15, 15, 15, 10, 15, 15, 15), 
                             T0.Code = c(46, 
                                         72, 61, 51, 62, 44, 64, 61, 39, 69, 43, 38, 63, 65, 25, 50, 
                                         61, 52, 51, 54, 48, 33, 75, 61, 35, 48, 73, 50, 56, 56, 66, 
                                         58, 64, 58, 39, 68, 80, 63, 42, 68, 64, 60, 56, 75, 70, 61, 
                                         74, 64, 54, 50, 63, 58, 64, 48, 72, 68, 50, 67, 40, 25, 55, 
                                         49, 57, 55, 63, 48, 67, 53, 56, 60, 33, 62, 77, 73, 59, 43, 
                                         52, 57, 70, 70, 68, 43, 55, 91, 68, 61, 75, 65, 65), 
                             T0.Score.Empan.dir = c(7, 
                                                    10, 6, 8, 6, 6, 7, 6, 6, 12, 6, 4, 7, 9, 4, 8, 10, 8, 4, 
                                                    8, 3, 7, 7, 10, 13, 5, 3, 12, 9, 7, 7, 5, 8, 9, 5, 5, 8, 
                                                    9, 8, 9, 8, 14, 5, 6, 8, 7, 8, 6, 8, 8, 5, 3, 9, 6, 7, 6, 
                                                    7, 9, 4, 4, 5, 6, 11, 7, 8, 2, 8, 5, 9, 8, 9, 7, 13, 9, 10, 
                                                    6, 9, 8, 7, 9, 10, 12, 7, 6, 9, 9, 8, 9, 7), 
                             T0.Score.Empan.inv = c(8, 
                                                    11, 5, 12, 4, 4, 7, 4, 6, 8, 6, 5, 4, 8, 4, 6, 6, 6, 4, 6, 
                                                    4, 8, 7, 2, 8, 6, 4, 9, 8, 6, 9, 9, 6, 6, 3, 7, 7, 9, 8, 
                                                    6, 4, 13, 5, 7, 6, 6, 7, 8, 7, 7, 6, 3, 11, 8, 8, 6, 9, 10, 
                                                    4, 5, 5, 7, 7, 9, 8, 4, 7, 7, 11, 7, 8, 4, 8, 10, 7, 7, 7, 
                                                    8, 7, 7, 11, 12, 5, 4, 7, 6, 9, 13, 9), 
                             T0.Fluences.phon = c(19, 
                                                  24, 15, 23, 24, 11, 17, 22, 21, 17, 16, 14, 10, 12, 13, 20, 
                                                  24, 23, 23, 22, 11, 12, 19, 21, 21, 17, 11, 24, 16, 20, 16, 
                                                  17, 23, 17, 14, 19, 14, 22, 11, 14, 17, 21, 16, 14, 15, 19, 
                                                  19, 16, 15, 21, 14, 11, 18, 16, 20, 10, 17, 18, 17, 21, 12, 
                                                  30, 13, 11, 15, 10, 12, 6, 24, 16, 10, 13, 22, 16, 14, 12, 
                                                  19, 12, 18, 18, 14, 25, 16, 20, 18, 10, 16, 19, 21), 
                             T0.Fluences.cat = c(27, 
                                                 28, 18, 23, 28, 13, 16, 10, 28, 27, 27, 14, 17, 19, 23, 24, 
                                                 21, 23, 24, 29, 25, 21, 30, 16, 20, 17, 31, 26, 27, 20, 18, 
                                                 25, 20, 26, 28, 28, 27, 31, 17, 25, 20, 24, 32, 23, 17, 27, 
                                                 28, 27, 21, 20, 25, 15, 28, 21, 18, 18, 27, 27, 16, 17, 16, 
                                                 22, 21, 17, 28, 15, 25, 18, 26, 24, 23, 20, 29, 33, 18, 16, 
                                                 15, 14, 17, 19, 24, 42, 22, 24, 26, 14, 25, 34, 20), 
                             Date.T1 = structure(c(10L, 
                                                   9L, 8L, 8L, 5L, 6L, 6L, 4L, 4L, 7L, 7L, 10L, 10L, 3L, 2L, 
                                                   13L, 12L, 13L, 24L, 11L, 15L, 23L, 38L, 31L, 46L, 46L, 40L, 
                                                   40L, 43L, 48L, 36L, 51L, 49L, 42L, 52L, 45L, 43L, 44L, 50L, 
                                                   49L, 41L, 47L, 52L, 41L, 38L, 38L, 37L, 30L, 19L, 19L, 19L, 
                                                   35L, 38L, 15L, 32L, 23L, 34L, 39L, 17L, 18L, 18L, 17L, 19L, 
                                                   16L, 20L, 33L, 21L, 21L, 27L, 27L, 22L, 25L, 25L, 26L, 27L, 
                                                   28L, 19L, 29L, 14L, 50L, 53L, 54L, 55L, 57L, 56L, 58L, 59L, 
                                                   60L, 1L), 
                                                 .Label = c("25/01/2020", "43123", "43131", "43145", 
                                                            "43146", "43149", "43150", "43152", "43153", "43156", "43504", 
                                                            "43516", "43517", "43525", "43534", "43542", "43544", "43545", 
                                                            "43548", "43550", "43552", "43556", "43561", "43562", "43565", 
                                                            "43567", "43568", "43570", "43571", "43575", "43577", "43579", 
                                                            "43583", "43600", "43603", "43605", "43609", "43610", "43612", 
                                                            "43613", "43614", "43615", "43619", "43621", "43624", "43628", 
                                                            "43629", "43632", "43634", "43635", "43636", "43637", "43646", 
                                                            "43683", "43705", "43708", "43737", "43780", "43868", "43871"
                                                 ), class = "factor"), 
                             T1.SRT.Rappel.moyen = c(12.36, 11.45, 
                                                     11.3, 12.5, 12.9, 8.5, 11.2, 11.78, 11.9, 11.2, 11.3, 9.2, 
                                                     10.5, 10.6, 11.18, 11.12, 7.73, 12.7, 11.28, 6.82, 12.13, 
                                                     12, 12.22, 11.8, 11.2, 12.43, 12.7, 12.22, 12.53, 12.4, 12.6, 
                                                     12.8, 12.42, 12.3, 11.8, 12.91, 12.96, 13.2, 12.32, 11.39, 
                                                     12.83, 11.5, 11.78, 13.6, 12.14, 12.5, 12.43, 12.44, 12.43, 
                                                     12.5, 13.33, 11.17, 11.71, 10.36, 11.36, 12.6, 8.54, 12.16, 
                                                     10.09, 10.64, 12.33, 12.67, 13.18, 10.85, 11.89, 10.82, 12.67, 
                                                     11.27, 11.83, 11.91, 13.33, 12.71, 11.5, 11.6, 12.25, 10.45, 
                                                     11.36, 10.91, 11, 11.66, 10.5, 12.83, 13, 12.43, 12.2, 11.7, 
                                                     12.83, 12.78, 14),
                             T1.SRT.rappel.diff = c(14, 14, 15, 15, 
                                                    14, 8, 13, 15, 15, 10, 13, 7, 11, 13, 12, 13, 13, 15, 14, 
                                                    7, 15, 15, 15, 11, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 
                                                    14, 14, 15, 15, 15, 14, 15, 15, 15, 15, 15, 15, 15, 15, 14, 
                                                    15, 15, 15, 15, 14, 15, 15, 10, 15, 12, 12, 15, 13, 15, 15, 
                                                    13, 15, 15, 12, 15, 15, 15, 14, 15, 15, 15, 11, 13, 14, 15, 
                                                    15, 12, 15, 15, 15, 14, 15, 15, 15, 15),
                             T1.Code = c(49, 
                                         68, 55, 58, 65, 54, 71, 64, 52, 67, 39, 44, 64, 68, 28, 61, 
                                         37, 41, 51, 34, 55, 34, 73, 67, 41, 52, 70, 53, 58, 56, 69, 
                                         62, 62, 60, 42, 71, 78, 66, 45, 67, 66, 62, 50, 73, 73, 63, 
                                         75, 65, 56, 54, 59, 61, 68, 53, 68, 70, 45, 65, 43, 28, 53, 
                                         50, 60, 53, 59, 50, 68, 62, 67, 38, 64, 67, 82, 74, 58, 39, 
                                         50, 53, 66, 69, 76, 50, 58, 85, 71, 68, 82, 68, 69), 
                             T1.Score.Empan.dir = c(7, 
                                                    11, 8, 10, 6, 8, 8, 10, 6, 12, 8, 7, 6, 10, 5, 10, 13, 9, 
                                                    4, 6, 4, 9, 8, 8, 13, 5, 6, 12, 8, 8, 7, 6, 7, 8, 4, 5, 9, 
                                                    8, 7, 8, 8, 13, 5, 5, 8, 7, 8, 7, 10, 5, 7, 5, 10, 6, 8, 
                                                    7, 8, 10, 6, 5, 5, 7, 12, 8, 10, 7, 7, 4, 9, 8, 8, 6, 10, 
                                                    7, 10, 7, 7, 7, 6, 10, 10, 11, 5, 6, 6, 7, 8, 12, 8), 
                             T1.Score.Empan.inv = c(7, 
                                                    14, 4, 10, 5, 6, 7, 7, 6, 6, 9, 6, 5, 6, 7, 9, 10, 8, 4, 
                                                    8, 5, 9, 8, 2, 7, 7, 4, 8, 9, 6, 6, 8, 6, 6, 5, 6, 8, 9, 
                                                    7, 5, 5, 13, 4, 8, 8, 7, 7, 7, 10, 8, 6, 4, 11, 7, 7, 7, 
                                                    8, 10, 5, 6, 3, 6, 10, 9, 9, 3, 10, 9, 9, 7, 12, 6, 9, 11, 
                                                    8, 7, 8, 7, 7, 6, 11, 10, 7, 5, 6, 9, 10, 13, 8), 
                             T1.Fluences.phon = c(25, 
                                                  21, 17, 22, 19, 16, 20, 21, 22, 21, 22, 12, 9, 18, 11, 15, 
                                                  12, 11, 20, 24, 17, 15, 16, 25, 19, 17, 20, 23, 22, 21, 20, 
                                                  20, 23, 17, 21, 18, 16, 24, 10, 14, 16, 20, 18, 14, 17, 21, 
                                                  19, 17, 12, 24, 17, 12, 19, 18, 19, 13, 19, 22, 15, 21, 12, 
                                                  29, 15, 14, 20, 6, 7, 9, 20, 15, 13, 16, 25, 17, 16, 20, 
                                                  15, 11, 18, 13, 15, 22, 11, 18, 20, 15, 14, 20, 19), 
                             T1.Fluences.cat = c(21, 
                                                 24, 21, 25, 23, 18, 23, 29, 23, 25, 17, 21, 21, 15, 20, 29, 
                                                 21, 20, 23, 30, 21, 23, 28, 13, 22, 20, 28, 27, 26, 23, 16, 
                                                 23, 21, 24, 25, 29, 28, 28, 18, 26, 20, 20, 33, 19, 23, 26, 
                                                 29, 28, 19, 17, 31, 18, 28, 26, 18, 22, 24, 21, 16, 15, 23, 
                                                 25, 24, 17, 29, 20, 28, 20, 31, 26, 20, 16, 30, 30, 15, 16, 
                                                 13, 14, 18, 25, 22, 35, 17, 20, 19, 15, 29, 33, 30) 
     ), class = "data.frame", row.names = c(NA, 89L))
     
     
     # add_pi require that the formula has been written 
     SRB.models4<-list()
     SRB.models4[[1]]<-glm(formula = T1.SRT.Rappel.moyen ~ Genre + Age + N.etudes  + T0.SRT.Rappel.moyen, 
                           family = Gamma(link = "identity"), data = data.lm )
     SRB.models4[[2]]<-glm(formula = T1.SRT.rappel.diff ~ Genre + Age + N.etudes +  
                             T0.SRT.rappel.diff, family = poisson(link = "identity"), data = data.lm)
     SRB.models4[[3]]<-glm(formula = T1.Code ~ Genre + Age + N.etudes + T0.Code, 
                           family = poisson(link = "identity"), data = data.lm)
     SRB.models4[[4]]<-glm(formula = T1.Score.Empan.dir ~ Genre + Age + N.etudes +  
                             T0.Score.Empan.dir, family = poisson(link = "identity"), data = data.lm)
     SRB.models4[[5]]<-glm(formula = T1.Score.Empan.inv ~ Genre + Age + N.etudes + 
                             T0.Score.Empan.inv, family = poisson(link = "identity"), data = data.lm)
     SRB.models4[[6]]<-glm(formula = T1.Fluences.phon   ~ Genre + Age + N.etudes + T0.Fluences.phon,
                           family = poisson(link = "identity"),  data = data.lm)
     SRB.models4[[7]]<-glm(formula = T1.Fluences.cat ~ Genre + Age + N.etudes + T0.Fluences.cat,  
                           family = poisson(link = "identity"), data = data.lm)
     

     NewO<-rbind(NewO, NewO)
     if(NewO$Genre[1]=="1") { 
       NewO[1,1]<-"H"
       NewO[2,1]<-"F"
       } else {
         NewO[1,1]<-"F"
         NewO[2,1]<-"H" }
     
     thresh<-qt(.95, 84)
     NewO$Genre<-factor(NewO$Genre)
     NewO<- NewO %>% mutate(Genre=forcats::fct_relevel(Genre, c("F","H"))) 

     fit2<-data.frame() # create data.frame for storing predicted values and CI
     for(i in 1:7){ # loop for the 7 tasks
       
       fit<- add_pi(df=NewO, SRB.models4[[i]], alpha=.1, names=c("lwr","upr"), yhatName="fit" )
       # select relevant variables 
       fit<-fit[1, c( i+10, 18:20)]
       
       
       # create data.frame with results
       names(fit)<-c( "T2","fit", "lwr", "upr")
       fit2<-rbind(fit2,fit)
     }
     
     fit2$se<-(fit2$fit-fit2$lwr)/thresh
     fit2$GSRB<-(fit2$T2-fit2$fit)/fit2$se
     
     GSRB<-round(fit2$GSRB,3)
     GSRB<-ifelse(GSRB>1.663197, paste0(GSRB, " ***"),
                  ifelse(GSRB<(-1.663197),GSRB(SRB4, " *"),paste0(GSRB, " **")) ) 
     
     

     
##########################################################          

          DT<-data.frame("Subtest" = label.tasks, 
                       "scaled score at T0"=z0, 
                       "scaled score at T1"=z1, 
                       "SD method" = SD_method,
                       "RCI (Chelune)"=RCI2,
                       "RCI (Iverson) "=RCI3,
                       "SRB (1 factor)" = c(SRB2),
                       "SRB (4 factors)" = c(SRB4),
                       "generalized SRB" = c(GSRB)
                       )
        DT<-as.matrix(DT)
        colnames(DT)<-c('<span style="color:white">Subtest</span>',
                        '<span style="color:white">scaled score at T0</span>',
                        '<span style="color:white">scaled score at T1</span>',
                        '<span style="color:white">SD method</span>',
                         '<span style="color:white">RCI (Chelune)</span>',
                         '<span style="color:white">RCI (Iverson)</span>',
                         '<span style="color:white">SRB (1 factor)</span>',
                         '<span style="color:white">SRB (4 factors)</span>',
                         '<span style="color:white">GSRB (4 factors)</span>'
                         )
        DT<-datatable(DT, escape = FALSE,  rownames=F, 
                      filter="none", autoHideNavigation=T, 
                      options = list(dom = 't'))%>% 
            formatStyle(1:11,backgroundColor = 'lightslategray', color="white", fontWeight="bold")
        
        
    })
    
    output$Figure1<-renderPlot({
      
      label.tasks<-c("SRT average recall","SRT delayed recall", "Digit symbol", 
                     "Forward digit span", "Backward digit span", "Phonemic fluency", "Semantic fluency")
      
      
      data.lm<-structure(list(ID = structure(c(9L, 31L, 18L, 62L, 67L, 22L, 
                                                24L, 35L, 16L, 59L, 78L, 61L, 81L, 88L, 58L, 49L, 82L, 46L, 48L, 
                                                11L, 26L, 7L, 86L, 66L, 54L, 56L, 17L, 51L, 14L, 60L, 50L, 63L, 
                                                15L, 68L, 69L, 30L, 40L, 27L, 13L, 87L, 53L, 25L, 10L, 43L, 39L, 
                                                71L, 89L, 34L, 4L, 5L, 6L, 42L, 72L, 45L, 64L, 8L, 41L, 47L, 
                                                79L, 38L, 1L, 80L, 12L, 44L, 20L, 70L, 23L, 19L, 77L, 76L, 75L, 
                                                33L, 65L, 32L, 3L, 57L, 85L, 74L, 2L, 21L, 73L, 83L, 37L, 55L, 
                                                84L, 36L, 29L, 28L, 52L), 
                                              .Label = c("ARNI54", "ARQU93", "AUMA93", "BAAL93", "BACA65", "BAJE64", "BAJF67", "BANA73", "BEJE68", "BLAU87", 
                                                         "BLBR54", "BOAN88", "BRCL55", "CHFR70", "CHJE86", "CHJO58", "CHMA89", "CHNA61", "CHYA78", "COAU86", "COGO69", "COJE70", "COPA84", "COVI73", 
                                                         "CRPA66", "DAFR67", "DAJE92", "DEAR74", "DELA94", "DESA78", "DOAN92","DUJE90", "DULE93", "FOGA74", "GELA60", "GOJE67", "GOME90", "GRCH54", 
                                                         "GRKA77", "HEVI74", "JAMA56", "JAST91", "LECA94", "LECL90", "LEMA65", "LOAN88", "LOAU93", "LOCH60", "LOPA54", "MACH93\r\n", "MAMY67", 
                                                         "MANA69", "MANO90", "MARI66", "MASA80", "MEIS66", "MODO59", "MOMA77","NACH62", "NIER70", "OLMI60", "ORIS63", "OUSA94", "PEAG88", "PEMA89", 
                                                         "PIJA93", "POFL81", "POGH64", "POMA67", "PRCA53", "PREM88", "PRMI85",  "QUGO93", "RIJA56", "RIPA66", "RISE74", "RISY59", "ROGI62", "SIDA54", 
                                                         "SIFR54", "SISO69", "STKU74", "SULA86", "TAMA61", "THMA61", "THSO74","VICH75", "VISA76", "VOCA81"), class = "factor"), 
                               DDN = structure(c(40L, 80L, 20L, 25L, 63L, 45L, 51L, 17L, 12L, 21L, 22L, 16L, 42L, 59L, 
                                                 61L, 5L, 54L, 8L, 13L, 7L, 37L, 35L, 53L, 83L, 33L, 31L, 75L, 
                                                 36L, 47L, 46L, 82L, 89L, 70L, 27L, 38L, 62L, 55L, 79L, 10L, 58L, 
                                                 77L, 34L, 71L, 87L, 60L, 74L, 64L, 57L, 86L, 30L, 26L, 78L, 69L, 
                                                 28L, 72L, 50L, 11L, 85L, 6L, 9L, 24L, 3L, 73L, 68L, 1L, 29L, 
                                                 39L, 67L, 48L, 52L, 32L, 4L, 76L, 23L, 2L, 49L, 18L, 66L, 81L, 
                                                 44L, 84L, 15L, 65L, 14L, 19L, 43L, 88L, 56L, 41L), 
                                               .Label = c("15/05/1986", "16/05/1993", "17/09/1954", "19/04/1993", "19814", "19885", "19937", 
                                                          "20/07/1988", "20033", "20219", "20786", "21240", "21974", "22/12/1980", 
                                                          "22/12/1986", "22059", "22091", "22352", "22423", "22640", "22724", 
                                                          "22957", "23/05/1990", "23/12/1954", "23107", "23505", "23533", 
                                                          "23896", "23957", "23985", "24137", "24270", "24276", "24303", 
                                                          "24475", "24494", "24496", "24700", "25/05/1984", "25002", "25269", 
                                                          "25286", "25408", "25453", "25833", "25836", "25922", "26/11/1956", 
                                                          "26/11/1959", "26666", "26865", "27063", "27074", "27189", "27264", 
                                                          "27345", "27388", "27558", "27941", "28281", "28318", "28503", 
                                                          "29706", "29857", "30/03/1990", "30/05/1956", "30/11/1978", "31/05/1990", 
                                                          "31367", "31723", "32132", "32288", "32396", "32463", "32616", 
                                                          "32758", "33002", "33564", "33622", "33868", "34014", "34081", 
                                                          "34213", "34250", "34271", "34324", "34368", "34396", "34414"), class = "factor"), 
                               Genre = structure(c(2L, 2L, 1L, 1L, 2L,  2L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 
                                                   2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 
                                                   2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L, 
                                                   1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 
                                                   1L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 
                                                   2L, 1L, 2L, 1L), .Label = c("F", "H"), class = "factor"), 
                               Age = c(49, 25, 56, 55, 36, 47, 44, 57, 59, 55, 55, 57, 48, 41, 40, 65, 44, 
                                       30, 58, 64, 51, 52, 45, 25, 52, 53, 29, 52, 48, 48, 25, 25, 32, 
                                       54, 52, 41, 45, 27, 64, 43, 28, 52, 31, 25, 41, 30, 37, 44, 25, 
                                       53, 54, 27, 33, 53, 30, 46, 62, 25, 65, 65, 65, 65, 30, 28, 32, 
                                       53, 35, 40, 59, 45, 52, 25, 29, 28, 25, 59, 58, 62, 25, 49, 25, 
                                       32, 29, 38, 58, 52, 25, 44, 50), 
                               N.etudes = c(15, 10, 12, 12, 15, 10, 10, 14, 15, 15, 15, 10, 14, 12, 15, 10, 13, 12, 8, 16, 
                                            11, 11, 11, 14, 11, 15, 12, 15, 11, 8, 15, 17, 12, 11, 11, 18, 
                                            14, 17, 17, 12, 15, 15, 12, 15, 12, 15, 12, 16, 15, 14, 16, 11, 
                                            17, 15, 15, 12, 14, 14, 14, 9, 10, 9, 14, 13, 14, 8, 12, 11, 
                                            10, 12, 10, 15, 17, 16, 12, 9, 11, 9, 17, 12, 14, 14, 11, 12, 
                                            17, 14, 16, 16, 12), 
                               Date.T0 = structure(c(13L, 12L, 11L, 10L, 
                                                     9L, 8L, 8L, 7L, 7L, 6L, 6L, 5L, 5L, 4L, 4L, 27L, 22L, 27L, 27L, 
                                                     23L, 14L, 24L, 41L, 40L, 45L, 45L, 41L, 45L, 42L, 37L, 38L, 38L, 
                                                     39L, 39L, 39L, 43L, 42L, 42L, 44L, 37L, 36L, 16L, 36L, 47L, 28L, 
                                                     28L, 38L, 27L, 17L, 20L, 21L, 39L, 38L, 15L, 31L, 24L, 40L, 40L, 
                                                     17L, 18L, 18L, 18L, 19L, 20L, 23L, 24L, 24L, 24L, 26L, 25L, 25L, 
                                                     30L, 30L, 32L, 32L, 33L, 33L, 35L, 26L, 42L, 29L, 3L, 34L, 49L, 
                                                     46L, 50L, 2L, 48L, 1L), 
                                                   .Label = c("16/11/2019", "18/10/2019", 
                                                              "23/02/2019", "43081", "43086", "43088", "43089", "43103", "43105", 
                                                              "43117", "43118", "43121", "43125", "43459", "43460", "43461", 
                                                              "43475", "43481", "43482", "43483", "43484", "43486", "43488", 
                                                              "43491", "43492", "43493", "43497", "43498", "43499", "43501", 
                                                              "43504", "43505", "43506", "43507", "43511", "43512", "43513", 
                                                              "43516", "43517", "43518", "43526", "43527", "43532", "43533", 
                                                              "43534", "43535", "43536", "43657", "43672", "43705"), class = "factor"), 
                               T0.SRT.Rappel.moyen = c(12.73, 12.6, 12, 10.5, 13.25, 9.3, 
                                                       12.5, 10, 11.73, 9.7, 11, 7.4, 10.3, 10.9, 12.09, 7.09, 12.5, 
                                                       11.18, 11.87, 11.27, 9.64, 11.82, 12.11, 10.91, 11.09, 12.25, 
                                                       12.5, 12.33, 12.67, 12.22, 12.7, 12.7, 12.44, 12.45, 11.9, 
                                                       12.8, 13, 13.17, 11.9, 11.44, 12.86, 11.4, 11.86, 13.5, 12.13, 
                                                       12.78, 12.75, 12.56, 11.18, 12, 13.4, 11.71, 11.75, 10.55, 
                                                       11.27, 12.83, 8.54, 11.5, 7.55, 11.09, 9.91, 10.91, 11.71, 
                                                       11.44, 12.33, 12.27, 12.75, 8.72, 11.5, 11.88, 12.19, 12.37, 
                                                       12.83, 11.86, 12.38, 9.18, 11.64, 10.35, 12, 12.57, 11.6, 
                                                       12.66, 12.17, 12.67, 12.38, 10, 13, 11.73, 12.71), 
                               T0.SRT.rappel.diff = c(15, 
                                                      15, 15, 14, 15, 8, 14, 13, 15, 11, 14, 5, 12, 13, 12, 4, 
                                                      15, 14, 13, 15, 12, 15, 15, 13, 15, 14, 15, 15, 15, 15, 15, 
                                                      15, 15, 15, 13, 15, 15, 15, 15, 14, 15, 15, 14, 15, 15, 15, 
                                                      15, 15, 14, 15, 15, 15, 14, 12, 15, 15, 7, 15, 10, 15, 11, 
                                                      13, 14, 14, 15, 13, 14, 9, 14, 14, 12, 15, 15, 15, 12, 7, 
                                                      13, 15, 14, 15, 15, 15, 15, 15, 15, 10, 15, 15, 15), 
                               T0.Code = c(46, 
                                           72, 61, 51, 62, 44, 64, 61, 39, 69, 43, 38, 63, 65, 25, 50, 
                                           61, 52, 51, 54, 48, 33, 75, 61, 35, 48, 73, 50, 56, 56, 66, 
                                           58, 64, 58, 39, 68, 80, 63, 42, 68, 64, 60, 56, 75, 70, 61, 
                                           74, 64, 54, 50, 63, 58, 64, 48, 72, 68, 50, 67, 40, 25, 55, 
                                           49, 57, 55, 63, 48, 67, 53, 56, 60, 33, 62, 77, 73, 59, 43, 
                                           52, 57, 70, 70, 68, 43, 55, 91, 68, 61, 75, 65, 65), 
                               T0.Score.Empan.dir = c(7, 
                                                      10, 6, 8, 6, 6, 7, 6, 6, 12, 6, 4, 7, 9, 4, 8, 10, 8, 4, 
                                                      8, 3, 7, 7, 10, 13, 5, 3, 12, 9, 7, 7, 5, 8, 9, 5, 5, 8, 
                                                      9, 8, 9, 8, 14, 5, 6, 8, 7, 8, 6, 8, 8, 5, 3, 9, 6, 7, 6, 
                                                      7, 9, 4, 4, 5, 6, 11, 7, 8, 2, 8, 5, 9, 8, 9, 7, 13, 9, 10, 
                                                      6, 9, 8, 7, 9, 10, 12, 7, 6, 9, 9, 8, 9, 7), 
                               T0.Score.Empan.inv = c(8, 
                                                      11, 5, 12, 4, 4, 7, 4, 6, 8, 6, 5, 4, 8, 4, 6, 6, 6, 4, 6, 
                                                      4, 8, 7, 2, 8, 6, 4, 9, 8, 6, 9, 9, 6, 6, 3, 7, 7, 9, 8, 
                                                      6, 4, 13, 5, 7, 6, 6, 7, 8, 7, 7, 6, 3, 11, 8, 8, 6, 9, 10, 
                                                      4, 5, 5, 7, 7, 9, 8, 4, 7, 7, 11, 7, 8, 4, 8, 10, 7, 7, 7, 
                                                      8, 7, 7, 11, 12, 5, 4, 7, 6, 9, 13, 9), 
                               T0.Fluences.phon = c(19, 
                                                    24, 15, 23, 24, 11, 17, 22, 21, 17, 16, 14, 10, 12, 13, 20, 
                                                    24, 23, 23, 22, 11, 12, 19, 21, 21, 17, 11, 24, 16, 20, 16, 
                                                    17, 23, 17, 14, 19, 14, 22, 11, 14, 17, 21, 16, 14, 15, 19, 
                                                    19, 16, 15, 21, 14, 11, 18, 16, 20, 10, 17, 18, 17, 21, 12, 
                                                    30, 13, 11, 15, 10, 12, 6, 24, 16, 10, 13, 22, 16, 14, 12, 
                                                    19, 12, 18, 18, 14, 25, 16, 20, 18, 10, 16, 19, 21), 
                               T0.Fluences.cat = c(27, 
                                                   28, 18, 23, 28, 13, 16, 10, 28, 27, 27, 14, 17, 19, 23, 24, 
                                                   21, 23, 24, 29, 25, 21, 30, 16, 20, 17, 31, 26, 27, 20, 18, 
                                                   25, 20, 26, 28, 28, 27, 31, 17, 25, 20, 24, 32, 23, 17, 27, 
                                                   28, 27, 21, 20, 25, 15, 28, 21, 18, 18, 27, 27, 16, 17, 16, 
                                                   22, 21, 17, 28, 15, 25, 18, 26, 24, 23, 20, 29, 33, 18, 16, 
                                                   15, 14, 17, 19, 24, 42, 22, 24, 26, 14, 25, 34, 20), 
                               Date.T1 = structure(c(10L, 
                                                     9L, 8L, 8L, 5L, 6L, 6L, 4L, 4L, 7L, 7L, 10L, 10L, 3L, 2L, 
                                                     13L, 12L, 13L, 24L, 11L, 15L, 23L, 38L, 31L, 46L, 46L, 40L, 
                                                     40L, 43L, 48L, 36L, 51L, 49L, 42L, 52L, 45L, 43L, 44L, 50L, 
                                                     49L, 41L, 47L, 52L, 41L, 38L, 38L, 37L, 30L, 19L, 19L, 19L, 
                                                     35L, 38L, 15L, 32L, 23L, 34L, 39L, 17L, 18L, 18L, 17L, 19L, 
                                                     16L, 20L, 33L, 21L, 21L, 27L, 27L, 22L, 25L, 25L, 26L, 27L, 
                                                     28L, 19L, 29L, 14L, 50L, 53L, 54L, 55L, 57L, 56L, 58L, 59L, 
                                                     60L, 1L), 
                                                   .Label = c("25/01/2020", "43123", "43131", "43145", 
                                                              "43146", "43149", "43150", "43152", "43153", "43156", "43504", 
                                                              "43516", "43517", "43525", "43534", "43542", "43544", "43545", 
                                                              "43548", "43550", "43552", "43556", "43561", "43562", "43565", 
                                                              "43567", "43568", "43570", "43571", "43575", "43577", "43579", 
                                                              "43583", "43600", "43603", "43605", "43609", "43610", "43612", 
                                                              "43613", "43614", "43615", "43619", "43621", "43624", "43628", 
                                                              "43629", "43632", "43634", "43635", "43636", "43637", "43646", 
                                                              "43683", "43705", "43708", "43737", "43780", "43868", "43871"
                                                   ), class = "factor"), 
                               T1.SRT.Rappel.moyen = c(12.36, 11.45, 
                                                       11.3, 12.5, 12.9, 8.5, 11.2, 11.78, 11.9, 11.2, 11.3, 9.2, 
                                                       10.5, 10.6, 11.18, 11.12, 7.73, 12.7, 11.28, 6.82, 12.13, 
                                                       12, 12.22, 11.8, 11.2, 12.43, 12.7, 12.22, 12.53, 12.4, 12.6, 
                                                       12.8, 12.42, 12.3, 11.8, 12.91, 12.96, 13.2, 12.32, 11.39, 
                                                       12.83, 11.5, 11.78, 13.6, 12.14, 12.5, 12.43, 12.44, 12.43, 
                                                       12.5, 13.33, 11.17, 11.71, 10.36, 11.36, 12.6, 8.54, 12.16, 
                                                       10.09, 10.64, 12.33, 12.67, 13.18, 10.85, 11.89, 10.82, 12.67, 
                                                       11.27, 11.83, 11.91, 13.33, 12.71, 11.5, 11.6, 12.25, 10.45, 
                                                       11.36, 10.91, 11, 11.66, 10.5, 12.83, 13, 12.43, 12.2, 11.7, 
                                                       12.83, 12.78, 14),
                               T1.SRT.rappel.diff = c(14, 14, 15, 15, 
                                                      14, 8, 13, 15, 15, 10, 13, 7, 11, 13, 12, 13, 13, 15, 14, 
                                                      7, 15, 15, 15, 11, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 
                                                      14, 14, 15, 15, 15, 14, 15, 15, 15, 15, 15, 15, 15, 15, 14, 
                                                      15, 15, 15, 15, 14, 15, 15, 10, 15, 12, 12, 15, 13, 15, 15, 
                                                      13, 15, 15, 12, 15, 15, 15, 14, 15, 15, 15, 11, 13, 14, 15, 
                                                      15, 12, 15, 15, 15, 14, 15, 15, 15, 15),
                               T1.Code = c(49, 
                                           68, 55, 58, 65, 54, 71, 64, 52, 67, 39, 44, 64, 68, 28, 61, 
                                           37, 41, 51, 34, 55, 34, 73, 67, 41, 52, 70, 53, 58, 56, 69, 
                                           62, 62, 60, 42, 71, 78, 66, 45, 67, 66, 62, 50, 73, 73, 63, 
                                           75, 65, 56, 54, 59, 61, 68, 53, 68, 70, 45, 65, 43, 28, 53, 
                                           50, 60, 53, 59, 50, 68, 62, 67, 38, 64, 67, 82, 74, 58, 39, 
                                           50, 53, 66, 69, 76, 50, 58, 85, 71, 68, 82, 68, 69), 
                               T1.Score.Empan.dir = c(7, 
                                                      11, 8, 10, 6, 8, 8, 10, 6, 12, 8, 7, 6, 10, 5, 10, 13, 9, 
                                                      4, 6, 4, 9, 8, 8, 13, 5, 6, 12, 8, 8, 7, 6, 7, 8, 4, 5, 9, 
                                                      8, 7, 8, 8, 13, 5, 5, 8, 7, 8, 7, 10, 5, 7, 5, 10, 6, 8, 
                                                      7, 8, 10, 6, 5, 5, 7, 12, 8, 10, 7, 7, 4, 9, 8, 8, 6, 10, 
                                                      7, 10, 7, 7, 7, 6, 10, 10, 11, 5, 6, 6, 7, 8, 12, 8), 
                               T1.Score.Empan.inv = c(7, 
                                                      14, 4, 10, 5, 6, 7, 7, 6, 6, 9, 6, 5, 6, 7, 9, 10, 8, 4, 
                                                      8, 5, 9, 8, 2, 7, 7, 4, 8, 9, 6, 6, 8, 6, 6, 5, 6, 8, 9, 
                                                      7, 5, 5, 13, 4, 8, 8, 7, 7, 7, 10, 8, 6, 4, 11, 7, 7, 7, 
                                                      8, 10, 5, 6, 3, 6, 10, 9, 9, 3, 10, 9, 9, 7, 12, 6, 9, 11, 
                                                      8, 7, 8, 7, 7, 6, 11, 10, 7, 5, 6, 9, 10, 13, 8), 
                               T1.Fluences.phon = c(25, 
                                                    21, 17, 22, 19, 16, 20, 21, 22, 21, 22, 12, 9, 18, 11, 15, 
                                                    12, 11, 20, 24, 17, 15, 16, 25, 19, 17, 20, 23, 22, 21, 20, 
                                                    20, 23, 17, 21, 18, 16, 24, 10, 14, 16, 20, 18, 14, 17, 21, 
                                                    19, 17, 12, 24, 17, 12, 19, 18, 19, 13, 19, 22, 15, 21, 12, 
                                                    29, 15, 14, 20, 6, 7, 9, 20, 15, 13, 16, 25, 17, 16, 20, 
                                                    15, 11, 18, 13, 15, 22, 11, 18, 20, 15, 14, 20, 19), 
                               T1.Fluences.cat = c(21, 
                                                   24, 21, 25, 23, 18, 23, 29, 23, 25, 17, 21, 21, 15, 20, 29, 
                                                   21, 20, 23, 30, 21, 23, 28, 13, 22, 20, 28, 27, 26, 23, 16, 
                                                   23, 21, 24, 25, 29, 28, 28, 18, 26, 20, 20, 33, 19, 23, 26, 
                                                   29, 28, 19, 17, 31, 18, 28, 26, 18, 22, 24, 21, 16, 15, 23, 
                                                   25, 24, 17, 29, 20, 28, 20, 31, 26, 20, 16, 30, 30, 15, 16, 
                                                   13, 14, 18, 25, 22, 35, 17, 20, 19, 15, 29, 33, 30) 
      ), class = "data.frame", row.names = c(NA, 89L))
      
      
      

      
      
      
      T0.subject<-c(input$SRT.av, input$SRT.d,input$DS, input$DDS,input$BDS,
                    input$PF,input$SF)
      T1.subject<-c(input$SRT.av.T1,input$SRT.d.T1, input$DS.T1,input$DDS.T1,input$BDS.T1,
                    input$PF.T1, input$SF.T1)
      Demo<-list(Age=input$Age, LoE= input$Educ, Sex=ifelse(input$Sex=="M","H", "F"))
      NewO<-data.frame( "Genre"=Demo$Sex,
                        "Age"=Demo$Age,
                        "N.etudes"=Demo$LoE,
                        "T0.SRT.Rappel.moyen"=T0.subject[1],
                        "T0.SRT.rappel.diff"=T0.subject[2],
                        "T0.Code"=T0.subject[3],
                        "T0.Score.Empan.dir"=T0.subject[4],
                        "T0.Score.Empan.inv"=T0.subject[5],
                        "T0.Fluences.phon"=T0.subject[6],
                        "T0.Fluences.cat"=T0.subject[7],
                        "T1.SRT.Rappel.moyen"=T1.subject[1],
                        "T1.SRT.rappel.diff"=T1.subject[2],
                        "T1.Code"=T1.subject[3],
                        "T1.Score.Empan.dir"=T1.subject[4],
                        "T1.Score.Empan.inv"=T1.subject[5],
                        "T1.Fluences.phon"=T1.subject[6],
                        "T1.Fluences.cat"=T1.subject[7])
 
      ES<-input$ES
      predicted<-c()
      lm.list<-list() # store linear models 
      if( ES== "linear model"){
        r.out<-c(0.310843341680445, 0.371973676629549, 0.727794681552961, 0.614199604579674, 
                 0.65947559562284, 0.473101976369185, 0.503209782036205)
        for(i in 1:7){
         
            f4<-as.formula(paste0(names(data.lm)[i+14], "~","Genre+Age+N.etudes+", 
                                  names(data.lm)[i+6]))
            
          # perform linear model with 1 and 4 factors
            lm.out<-lm(f4, data.lm)
            lm.list[[i]]<-lm.out
            predicted<-c(predicted, predict.lm( lm.out, NewO))
        }
        
      }else{
        r.out<-c(0.310574414104724, 0.371850374423101, 0.72758339262502, 0.614109387718531, 
                 0.657675334561725, 0.472487624745481, 0.502110199777427) 
        
        for(i in 1:7){

          f4<-as.formula(paste0(names(data.lm)[i+14], "~","Genre+Age+N.etudes+", 
                                names(data.lm)[i+6]))
          
          # perform linear model with 1 and 4 factors
          lm.out<-glm(f4, data.lm, family  =poisson(link="identity"))
          if(i == 1) lm.out<-glm(f4, data.lm, family =Gamma(link="identity"))
          lm.list[[i]]<-lm.out
          predicted<-c(predicted, predict.glm( lm.out, NewO))
        }
        
      }
      
      

      # compute the three effect size
      
    
      
      
      
      zop<-mapply(function(o,data, p,r2){
        s=sd(data)
        (o-p)/(s*(1-r2))}, o=NewO[,11:17] , data=data.lm[,c('T1.SRT.Rappel.moyen','T1.SRT.rappel.diff','T1.Code',
                                         'T1.Score.Empan.dir',
                                         'T1.Score.Empan.inv','T1.Fluences.phon', 'T1.Fluences.cat')],
        p= predicted ,r2= r.out )
      zop<-round(zop,3)
     
     
      
      p1<-effect_plot(lm.list[[1]], pred = T0.SRT.Rappel.moyen,  interval = TRUE)+xlab("T0") + 
        ylab("T1")+ geom_point(aes(x=T0.subject[1], y=T1.subject[1]), colour="blue", size=5)+
        annotate(geom="text", label= paste0("Zop =", zop[1]), x=8, y=12)
      
      p2<-effect_plot(lm.list[[2]], pred = T0.SRT.rappel.diff,  interval = TRUE)+xlab("T0") + 
        ylab("T1")+ geom_point(aes(x=T0.subject[2], y=T1.subject[2]), colour="blue", size=5)+
        annotate(geom="text", label= paste0("Zop =", zop[2]), x=5, y=12)
      
      p3<-effect_plot(lm.list[[3]], pred = T0.Code,  interval = TRUE)+xlab("T0") + 
        ylab("T1")+ geom_point(aes(x=T0.subject[3], y=T1.subject[3]), colour="blue", size=5)+
        annotate(geom="text", label= paste0("Zop =", zop[3]), x=25, y=75)
      
      p4<-effect_plot(lm.list[[4]], pred = T0.Score.Empan.dir,  interval = TRUE)+xlab("T0") + 
        ylab("T1")+ geom_point(aes(x=T0.subject[4], y=T1.subject[4]), colour="blue", size=5)+
        annotate(geom="text", label= paste0("Zop =", zop[4]), x=4, y=11)
      
      p5<-effect_plot(lm.list[[5]], pred = T0.Score.Empan.inv,  interval = TRUE)+xlab("T0") + 
        ylab("T1")+ geom_point(aes(x=T0.subject[5], y=T1.subject[5]), colour="blue", size=5)+
        annotate(geom="text", label= paste0("Zop =", zop[5]), x=4, y=11)
      
      p6<-effect_plot(lm.list[[6]], pred = T0.Fluences.phon,  interval = TRUE)+xlab("T0") + 
        ylab("T1")+ geom_point(aes(x=T0.subject[6], y=T1.subject[6]), colour="blue", size=5)+
        annotate(geom="text", label= paste0("Zop =", zop[6]), x=10, y=25)
      
      p7<-effect_plot(lm.list[[7]], pred = T0.Fluences.cat,  interval = TRUE)+xlab("T0") + 
        ylab("T1")+ geom_point(aes(x=T0.subject[7], y=T1.subject[7]), colour="blue", size=5)+
        annotate(geom="text", label= paste0("Zop =", zop[7]), x=10, y=25)
      
      ggarrange(p1, p2, p3, p4, p5, p6, p7, 
                labels =label.tasks, hjust=-2, vjust=3, common.legend = F,
                font.label= list(size=18, face="plain"),
                ncol=2, nrow=4)
     
      
      
     
       
      
    })
    
    
}


# Run the application 
shinyApp(ui = ui, server = server)



