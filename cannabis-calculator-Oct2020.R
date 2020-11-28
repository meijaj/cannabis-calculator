##################################################
# Thermal stability of dry cannabis: a kinetic study 
# Author: Juris Meija, NRC Canada
# Date: 2020, version October 2020
##################################################

# Brief description
# This calculator simulates the degradation profile of cannabinoids in dried cannabis samples stored in the dark. 
# A network of consecutive first-order reactions is assumed and the values for the rate constants are 
# obtained from the stability study at the NRC Canada involving dried cannabis reference material (2019-2020).
# This work is an interactive supplement to publication in Anal Bioanal Chem (2020) doi: 10.1007/s00216-020-03098-2

require(shiny)    # shiny
require(shinyjs)  # java-shiny for collapsible features
require(mkin)     # kinetic model integrator
require(mvtnorm)  # draw samples of model parameters from multivariate normal distribution

# Molar masses of cannabinoids
M = list(id=c('THCA','THC','CBNA','CBN','CBDA','CBD','CBG','CBGA','CBC','CBCA','THCV','THCVA','CBDV','CBDVA'),
         M=c(358,314,354,310,358,314,316.5,360.5,314.5,357.5,286.4,330.4,286.4,330.4))

# Transparent gray (50% gray60)
gray.transp = rgb(col2rgb("gray60")[1], col2rgb("gray60")[2], col2rgb("gray60")[3], max = 255, alpha = (100 - 50) * 255 / 100)

### server file
server <- function(input, output, session) {

  # Reactive values
  values <- reactiveValues(x = NULL, A = NULL, N = NULL, q=NULL, can=NULL, D=NULL, U=NULL,
                           unit=NULL, A.q = NULL, N.q = NULL, AN.q = NULL, times = NULL, temperature = NULL)

  observeEvent(input$cannabinoid, {
    if(input$cannabinoid!='CBN') updateTextInput(session, 'initialvalueA', label = paste0('Initial value of ',input$cannabinoid,'A (mg/g)'))
    if(input$cannabinoid!='CBN') updateTextInput(session, 'initialvalueN', label = paste0('Initial value of ',input$cannabinoid,' (mg/g)'))
    if(input$cannabinoid=='CBN') updateTextInput(session, 'initialvalueA', label = paste0('Initial value of THCA (mg/g)'))
    if(input$cannabinoid=='CBN') updateTextInput(session, 'initialvalueN', label = paste0('Initial value of THC (mg/g)'))
    # Set initial cannabinoid values to the NRC Certified Reference Material of cannabis
    updateTextInput(session, 'initialvalueA', value = c(125, 125, 24, 1.9, 4.3, 0.13, 0.49)[input$cannabinoid == c('THC','CBN','CBD','CBC','CBG','CBDV','THCV')])
    updateTextInput(session, 'initialvalueN', value = c(60, 60, 9.1, 1.2, 2.1, 0.046, 0.30)[input$cannabinoid == c('THC','CBN','CBD','CBC','CBG','CBDV','THCV')])
    updateTextInput(session, 'initialvalueCBNA', value = 2.0, label = paste0('Initial value of CBNA (mg/g)'))
    updateTextInput(session, 'initialvalueCBN', value = 2.9, label = paste0('Initial value of CBN (mg/g)')) 
  })

  observeEvent(input$go, {
    # validate inputs
    validate(
      need(input$initialvalueA>=0, ''),
      need(input$initialvalueN>=0, ''),
      need(input$initialvalueCBNA>=0, ''),
      need(input$initialvalueCBN>=0, ''),
      need(input$temperature>=0&input$temperature<100, ''),
      need(input$duration>=0, ''),
      need(input$out>=50, ''),
      need(input$mc>=100, '')
    )
  
  withProgress(message = 'Please wait...', value = 0, {
    
    values$can <- input$cannabinoid
    values$D <- input$duration
    values$U <- input$durationunit
    values$temperature <- input$temperature
    
    # Duration (time) in weeks
    time.max = values$D*c(1/(7*24), 1/7, 1, 4.345, 4.345*12)[match(values$U, c('hour','day','week','month','year'))]
    
    # Initial amount of cannabinoids (mol)
    MA = M$M[M$id == paste0(values$can,'A')]
    MN = M$M[M$id == values$can]
    
    c0.A = input$initialvalueA/MA
    c0.N = input$initialvalueN/MN
    if(values$can=='THC'|values$can=='CBN') c0.CBNA = input$initialvalueCBNA/354
    if(values$can=='THC'|values$can=='CBN') c0.CBN = input$initialvalueCBN/310
    
    Nmc  <- input$mc  # number of Monte-Carlo simulations
    Nout <- input$out # number of time points
    
    # Calculate the rate constants for the given cannabinoid and temperature
    if(values$can=='THC'|values$can=='CBN'){
      #fit-THC-CBN-20Oct2020.rds
      alpha = c(64.43896, 65.64605, 42.69356, 20.07184, 21.43670)
      beta = c(-20696.73217, -21616.79421, -13989.66665,  -8320.80883,  -8194.69217)
      cov.aabb = matrix(c(35.5, 102, -2.2, -3.7, 2.2, -11034.9, -31783.8, 655.2, 1082.3, -685.5, 102, 788.7, -17.3, -0.7, 6.1, -31821.2, -245949, 5504.4, 431.3, -1936.1, -2.2, -17.3, 80.4, -3.1, -20.2, 719.3, 5441.5, -25557.7, -51.9, 6396.4, -3.7, -0.7, -3.1, 16.6, -3.6, 1145.9, 206.1, 1160.7, -4784.1, 1086.8, 2.2, 6.1, -20.2, -3.6, 16.4, -697.9, -1921.1, 6416.7, 1356, -5132.8, -11034.9, -31821.2, 719.3, 1145.9, -697.9, 3433156.1, 9915862.4, -212548.7, -334953.4, 216267, -31783.8, -245949, 5441.5, 206.1, -1921.1, 9915862.4, 76696990.3, -1734148.3, -132519.1, 610432.6, 655.2, 5504.4, -25557.7, 1160.7, 6416.7, -212548.7, -1734148.3, 8130927.3, -30560, -2037500.9, 1082.3, 431.3, -51.9, -4784.1, 1356, -334953.4, -132519.1, -30560, 1391651.2, -414110.2, -685.5, -1936.1, 6396.4, 1086.8, -5132.8, 216267, 610432.6, -2037500.9, -414110.2, 1610307.1), ncol=10)
    }
    if(values$can=='CBD'){
      #fit-CBD-20Oct2020.rds
      alpha = c(46.97385, 48.21181)
      beta = c(-15382.51161, -16258.14084)
      cov.aabb = matrix(c(9.1, 33.3, -2823.8, -10341.4,
                          33.3, 302.4, -10350.5, -94190.2,
                          -2823.8, -10350.5, 876688.5, 3217891.5,
                          -10341.4, -94190.2, 3217891.5, 29342027.9), ncol=4)
    }
    if(values$can=='CBC'){
      #fit-CBC-20Oct2020.rds
      alpha = c(46.97385, 48.21181)
      beta = c(-15382.51161, -16258.14084)
      cov.aabb = matrix(c(9.1, 33.3, -2823.8, -10341.4,
                          33.3, 302.4, -10350.5, -94190.2,
                          -2823.8, -10350.5, 876688.5,  3217891.5,
                          -10341.4, -94190.2, 3217891.5, 29342027.9), ncol=4)
    }
    if(values$can=='CBG'){
      #fit-CBG-20Oct2020.rds
      alpha = c(47.48775, 42.85833)
      beta = c(-15590.44528, -14513.67998)
      cov.aabb = matrix(c(28.1, 66.7, -8727.8,  -20763.8,
                          66.7, 279.1, -20776.5, -86950.4,
                          -8727.8, -20776.5, 2713358.9,  6465087.4,
                          -20763.8, -86950.4, 6465087.4, 27087593.9), ncol=4)
    }
    if(values$can=='CBDV'){
      #fit-CBDV-20Oct2020.rds
      alpha = c(50.15791, 45.86154)
      beta = c(-16372.91, -15415.05)
      cov.aabb = matrix(c(11.78, 26.89, -3656.61, -8359.25, 26.89, 136.39, -8374.21, -42501.1, -3656.61, -8374.21, 1135633.84, 2603376.58, -8359.25, -42501.1, 2603376.58, 13245450.67), ncol=4)
    }
    if(values$can=='THCV'){
      #fit-THCV-20Oct2020.rds
      alpha = c(54.70056, 28.76249)
      beta = c(-17627.03, -10010.39)
      cov.aabb = cov.aabb = matrix(c(7.51, 6.93, -2334.13, -2155.34, 6.93, 52.24, -2155.82, -16249.81, -2334.13, -2155.82, 725370.83, 670590, -2155.34, -16249.81, 670590, 5055138.14), ncol=4)
    }  
    
    ab.mc = unlist(rmvnorm(Nmc, mean=c(alpha, beta), sigma=cov.aabb, method = 'svd'))
    
    tinv = c(1, 1/(input$temperature + 273.15))
    if(values$can=='THC'|values$can=='CBN') {lnk.mc = cbind(ab.mc[,c(1,6)]%*%tinv, ab.mc[,c(2,7)]%*%tinv, ab.mc[,c(3,8)]%*%tinv, ab.mc[,c(4,9)]%*%tinv, ab.mc[,c(5,10)]%*%tinv )}
    else { lnk.mc = cbind(ab.mc[,c(1,3)]%*%tinv, ab.mc[,c(2,4)]%*%tinv) }
    
    ### Kinetic model
    incProgress(0.10, detail = "Compiling the kinetic model")
    if(values$can=='THC'|values$can=='CBN'){   
      model = mkinmod(parent = mkinsub("SFO", c("m1", "m2"), sink = FALSE),
                      m1 = mkinsub("SFO", c("m3")),
                      m2 = mkinsub("SFO", c("m3"), sink = FALSE),
                      m3 = mkinsub("SFO", sink = FALSE), use_of_ff = "min", quiet = TRUE)
    }
    else{
      model = mkinmod(parent = mkinsub("SFO", "m1", sink=FALSE), m1 = mkinsub("SFO"), use_of_ff = "min", quiet = TRUE)
    }   
    ### Monte Carlo for uncertainties
    model.res.mc.A = array(dim=c(Nout,Nmc))
    model.res.mc.N = array(dim=c(Nout,Nmc))
    model.res.mc.CBNA = array(dim=c(Nout,Nmc))
    model.res.mc.CBN = array(dim=c(Nout,Nmc))
    
    for(i in 1:Nmc){
      incProgress(1/Nmc, detail = paste("Monte Carlo draw", i))
      
      if(input$cannabinoid=='THC'|input$cannabinoid=='CBN'){
        model.res.mc = data.frame(mkinpredict(model, odeparms = c(k_parent_m1 = exp(lnk.mc[i,1]), k_m1_sink = exp(lnk.mc[i,2]),
                                                                  k_m2_m3 = exp(lnk.mc[i,3]), 
                                                                  k_parent_m2 = exp(lnk.mc[i,4]), k_m1_m3 = exp(lnk.mc[i,5])), 
                                              odeini = c(parent = c0.A, m1 = c0.N, m2 = c0.CBNA, m3 = c0.CBN), 
                                              outtimes = seq(0, time.max, length.out = Nout), solution_type = "eigen")
        )      
        model.res.mc.A[,i] = model.res.mc[,2]
        model.res.mc.N[,i] = model.res.mc[,3]
        model.res.mc.CBNA[,i] = model.res.mc[,4]
        model.res.mc.CBN[,i] = model.res.mc[,5]
      }
      else{
        model.res.mc = data.frame(mkinpredict(model, odeparms = c(k_parent_m1 = exp(lnk.mc[i,1]), k_m1_sink = exp(lnk.mc[i,2])), 
                                              odeini = c(parent = c0.A, m1 = c0.N), 
                                              outtimes = seq(0, time.max, length.out = Nout), solution_type = "analytical"))
        model.res.mc.A[,i] = model.res.mc[,2]
        model.res.mc.N[,i] = model.res.mc[,3]
      }
    }
    
    ### Average rate constants (lnk)
    lnk = apply(lnk.mc,2,mean)
    
    # Calculate 20%, 50%, and 80% quantiles
    values$A.q <- MA*apply(model.res.mc.A, 1, quantile, c(0.2, 0.5, 0.8) )
    values$N.q <- MN*apply(model.res.mc.N, 1, quantile, c(0.2, 0.5, 0.8) )
    values$AN.q <- MN*apply(model.res.mc.A+model.res.mc.N, 1, quantile, c(0.2, 0.5, 0.8) )
    if(values$can=='CBN') {
      values$A.q <- 354*apply(model.res.mc.CBNA, 1, quantile, c(0.2, 0.5, 0.8) )
      values$N.q <- 310*apply(model.res.mc.CBN, 1, quantile, c(0.2, 0.5, 0.8) )
      values$AN.q <- 310*apply(model.res.mc.CBNA+model.res.mc.CBN, 1, quantile, c(0.2, 0.5, 0.8) )
    }
    ##### END Monte Carlo
    values$times <- seq(0, values$D, length.out = Nout) # Duration in user units
    
  })
  
   # Calculate cannabis model age
    delta_k = exp(lnk[1]) - exp(lnk[2]) + ifelse(values$can=='THC', exp(lnk[4]) - exp(lnk[5]), 0)
    R_NA = (input$initialvalueN/MN)/(input$initialvalueA/MA)
    age = (1/delta_k)*log(1 + R_NA*delta_k/exp(lnk[1]))
    
    output$age = renderText({
      paste0('<b>The model age of such cannabis sample stored in dark conditions at ',values$temperature,' &deg;C temperature is ', round(age),' weeks.</b><br>In other words, it takes that long to attain the specified ',values$can,'A:',values$can, ' ratio at ',values$temperature,' &deg;C temperature from the time cannabis is harvested.')
    })
    # Data
    output$data1 = renderText({ paste0('<b>The following Arrhenius parameters are used</b>') })
    if(values$can=='THC'|values$can=='CBN') {id.n = c('THCA->THC','THC->X','CBNA->CBN','THCA->CBNA','THC->CBN')}
    else {id.n = c(paste0(values$can,'A->',values$can), paste0(values$can,'->X'))}
    output$datatable = renderTable({
      cbind('Reaction'=id.n,'a'=round(alpha,0), 'u(a)'=round(sqrt(diag(cov.aabb)[1:length(alpha)]),0), 'b'=round(beta,0), 'u(b)'=round(sqrt(diag(cov.aabb)[-c(1:length(alpha))]),0)) 
    })
    output$data2 = renderText({ paste0('The rate constants given by k = exp(a + b/(T/K)) have units of 1/week.') })  
  })
  
  # Plot results
  observeEvent(input$plot_hover$x, {
      if(!is.null(values$A.q)){
        values$x <- input$plot_hover$x
        id <- which.min(abs(input$plot_hover$x - values$times))
        dig <- ifelse(values$AN.q[2,1]>10,1,(3-0.4*log(values$AN.q[2,1]))%/%1) # establish significant digits
        values$A <-  formatC(values$A.q[2,id],dig,format='f')
        values$N <-  formatC(values$N.q[2,id],dig,format='f')
        values$AN <- formatC(values$AN.q[2,id],dig,format='f')
        
        values$Arel <-  ifelse(values$A.q[2,1]==0,'n/a',formatC(100*values$A.q[2,id]/values$A.q[2,1], 0, format='f'))
        values$Nrel <-  ifelse(values$N.q[2,1]==0,'n/a',formatC(100*values$N.q[2,id]/values$N.q[2,1], 0, format='f'))
        values$ANrel <- formatC(100*values$AN.q[2,id]/values$AN.q[2,1], 0,format='f')    
      }
      
    })
    
  output$plot <- renderPlot({
      if( !is.null(values$A.q) ){
        par(mar=c(5,10,6,9))
        plot(x=values$times, y=values$A.q[2,], type='l', col='steelblue3', lwd=3,
             ylim=c(min(values$A.q[1,], values$N.q[1,]), max(values$AN.q[3,])),
             xlim=c(0, values$D), cex.lab=1.5, cex.axis=1.5, main='', 
             xlab = paste0('Storage time (', values$U,'s) at ', values$temperature, ' \u00B0C'), 
             ylab = "Mass fraction (mg/g)"  )
        abline(v=values$x, lty=2) # Vertical hover line
        # Uncertainty bound and mean for the acid cannabinoid
        if(input$uncertainty=='Yes') polygon(c(values$times, rev(values$times)), c(values$A.q[1,], rev(values$A.q[3,])), col=gray.transp, border=gray.transp)
        lines(x=values$times, y=values$A.q[2,], col='steelblue3', lwd=3)
        
        # Uncertainty bound and mean for neutral cannabinoid
        if(input$uncertainty=='Yes') polygon(c(values$times, rev(values$times)), c(values$N.q[1,], rev(values$N.q[3,])), col=gray.transp, border=gray.transp)
        lines(x=values$times, y=values$N.q[2,], col='tomato', lwd=3)
        
        # Uncertainty bound and mean for total cannabinoid equivalent
        if(input$uncertainty=='Yes') polygon(c(values$times, rev(values$times)), c(values$AN.q[1,], rev(values$AN.q[3,])), col=gray.transp, border=gray.transp)
        lines(x=values$times, y=values$AN.q[2,], col='gray40', lwd=3)
        
        # Hover info text
        if(!is.null(values$x)) if(values$x>=0) if(values$x <= values$D){
          mtext(text=paste0(values$can,'A: ', values$A,' mg/g (',values$Arel,' % initial value)'), side=3, line=2.0, at=values$x, adj=0.5, cex=1.3, font=2, col='steelblue3')
          mtext(text=paste0(values$can,': ', values$N,' mg/g (',values$Nrel,' % initial value)'), side=3, line=0.6, at=values$x, adj=0.5, cex=1.3, font=2, col='tomato')
          mtext(text=paste0('total ', values$can, ': ', values$AN,' mg/g (',values$ANrel,' % initial value)'), side=3, line=3.4, at=values$x, adj=0.5, cex=1.3, font=2, col='gray40')
          points(x=rep(values$x, 3), y=c(values$A, values$N, values$AN), pch=21, bg='white', col=c('steelblue3','tomato','gray40'), cex=2.3)
        }
        graphics::box(col='black')
      }
      else {return(NULL)}
    })
  
  shinyjs::onclick("toggleextra", shinyjs::toggle(id = "filterextra", anim = TRUE))
  shinyjs::onclick("togglenotes", shinyjs::toggle(id = "filternotes", anim = TRUE))
  shinyjs::onclick("toggleextra.mol", shinyjs::toggle(id = "filterextra", anim = TRUE))
  shinyjs::onclick("togglenotes.mol", shinyjs::toggle(id = "filternotes", anim = TRUE))
  
}

### user interface
ui <- fluidPage(
  tags$style("@import url(https://use.fontawesome.com/releases/v5.7.2/css/all.css);"),
  titlePanel( title="Cannabis stability calculator" ),
  shinyjs::useShinyjs(),
  sidebarLayout(
    
    sidebarPanel(
      
      fluidRow(column(12, selectInput("cannabinoid", label = tags$div(HTML(' <i class="fa fa-cannabis" style = "color:#808000;"></i> Select cannabinoid')), choices = c('THC','CBD','CBN','CBC','CBG','CBDV','THCV'), selected = 'THC')) ),
      fluidRow(column(12, h5(tags$div(HTML(' <i class="fas fa-exclamation" style = "color:#ff6633;"></i> Note that 10 mg/g = 1 %'))
                             )
                      ) 
               ),
      fluidRow(column(6, numericInput("initialvalueA", label = "Initial value of THCA (mg/g)", value = 125, min = 0, step = 1)),
               column(6, numericInput("initialvalueN", label = "Initial value of THC (mg/g)", value = 60, min = 0, step = 1))
               ),
      conditionalPanel('input.cannabinoid=="THC"|input.cannabinoid=="CBN"',
                       fluidRow(column(6, numericInput("initialvalueCBNA", label = "Initial value of CBNA (mg/g)", value = 2.0, min = 0, step = 1)),
                                column(6, numericInput("initialvalueCBN", label = "Initial value of CBN (mg/g)", value = 2.9, min = 0, step = 1))
                       )),
      br(),br(),
      fluidRow(column(6, numericInput("temperature", label = tags$div(HTML('<i class="fa fa-thermometer-quarter" style = "color:#0072B2;"></i> Temperature (&deg;C)')), value = 40, min = 4, max = 100, step = 1)) ),
      fluidRow(column(6, numericInput("duration", label = tags$div(HTML('<i class="fa fa-clock" style = "color:#0072B2;"></i> Duration')), value = 1, min = 0, step = 1)),
               column(6, selectInput("durationunit", label = "Unit", choices = c('year','month','week','day','hour'), selected = 'month'))
      ),
      h5("Additional settings:", a(id = "toggleextra", "show/hide")),
      shinyjs::hidden(div(id = "filterextra", 
                      fluidRow(column(6, radioButtons('uncertainty', label='Show 80 % uncertainty bounds', choices = c('Yes','No'), selected = 'No', inline=TRUE)) ),
                      fluidRow(column(6, numericInput("out", label = 'Number of points to draw curves', value = 100, min=50, max=500)),
                      column(6, numericInput("mc", label = 'Number of Monte Carlo draws', value = 500, min=500, max=5000))
                      )
                      )
                ),
      br(),
      actionButton("go", "Calculate!")
      
    ),
    
    mainPanel(
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      plotOutput("plot", hover = hoverOpts(id = "plot_hover", delay = 10, delayType = "debounce")),br(),
      htmlOutput("data1"), tableOutput("datatable"),
      htmlOutput("data2"), br(),
      htmlOutput("age"), br(),
      h5("Notes and explanations:", a(id = "togglenotes", "show/hide")),
      shinyjs::hidden(div(id = "filternotes",
                          fluidRow(column(
                            p("This calculator simulates the degradation profiles of cannabinoids in dried cannabis samples stored in the dark. A network of consecutive first-order reactions is assumed and the values for the rate constants are taken from the stability study at the NRC involving dried cannabis reference material (2019-2020). The following defintions of duration are adopted: 1 year = 12 months = 365 days. Uncertainty of the predicted profiles is evaluated from the model parameters using the Monte Carlo method."),
                            p("Created using R and Shiny. Juris Meija (2020) NRC Canada"), width=11)
                          ))),
      br(), 
      h5("NRC Cannabis stability calculator (October 2020)", a("See Anal Bioanal Chem (2020) publication", href = "https://doi.org/10.1007/s00216-020-03098-2"),"by Meija et al  for more details.")
      
    )
  )
)

shinyApp(ui = ui, server = server)
