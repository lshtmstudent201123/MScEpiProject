home<-"~/Desktop/LSHTM/Thesis/"
setwd(home)

?ode
#Basic model to test function#
library(deSolve)
library(tidyverse)
new.basic.model <- function(t, state, parms) { 
  with(as.list(c(state, parms)),{ 
    #Parameters
    t <- parms["t"]
    a1 <- parms["a1"]
    a2 <- parms["a2"]
    a3 <- parms["a3"]
    a4 <- parms["a4"]
    a5 <- parms["a5"]
    a6 <- parms["a6"]
    a7 <- parms["a7"]
    d <- parms["d"]
    D <- parms["D"]
    L <-parms["L"]
    D1 <- parms["D1"]
    D2 <- parms["D2"]
    D3 <- parms["D3"]
    NS <- parms["NS"]
    e <- parms["e"]
    e1 <- parms["e1"]
    e2 <- parms["e2"]
    e3 <- parms["e3"]
    e4 <- parms["e4"]
    e5 <- parms["e5"]
    e6 <- parms["e6"]
    e7 <- parms["e7"]
    e8 <- parms["e8"]
    gamma <- parms["gamma"]
    R1 <- parms["R1"]
    R2 <- parms["R2"]
    R3 <- parms["R3"]
    R4 <- parms["R4"]
    R5 <- parms["R5"]
    R6 <- parms["R6"]
    SA <- parms["SA"]
    SB <- parms["SB"]
    SC <- parms["SC"]
    
    #State Basics
    WA <- state["WA"]
    WB <- state["WB"]
    WC <- state["WC"]
    TA <- state["TA"]
    TB <- state["TB"]
    TC <- state["TC"]
    NE <- state["NE"]
    ND <- state["NE"]
    NR <- state["NE"]
    TA <- state["TA"]
    TB <- state["TB"]
    TC <- state["TC"]
    h <- state["h"]
    
    
    #Tree Diagram states
    AES <- state["AES"]
    AESD <- state["AESD"]
    AESO  <- state["AESO"]
    AESOR <- state["AESOR"]
    AESOD <- state["AESOD"]
    AESOEF <- state["AESOEF"]
    AESOES <- state["AESOES"]
    AESOEFE <- state["AESOEFE"]
    AESOEFO <- state["AESOEFO"]
    AESOEFD <- state["AESOEFD"]
    AESOESD <- state["AESOESD"]
    AESOESO <- state["AESOESO"]
    AESOEFED <- state["AESOEFED"]
    AESOEFEO <- state["AESOEFEO"]
    AESOEFOR <- state["AESOEFOR"]
    AESOEFOD <- state["AESOEFOD"]
    AESOEFOE <- state["AESOEFOE"]
    AESOESOR <- state["AESOESOR"]
    AESOESOD <- state["AESOESOD"]
    AESOEFEOR <- state["AESOEFEOR"]
    AESOEFEOD <- state["AESOEFEOD"]
    AESOEFOED <- state["AESOEFOED"]
    AESOEFOEO <- state["AESOEFOEO"]
    AESOEFOEOR <- state["AESOEFOEOR"]
    AESOEFOEOD <- state["AESOEFOEOD"]
    
    AESEF  <- state["AESEF"]
    AESEFE <- state["AESEFE"]
    AESEFO <- state["AESEFO"]
    AESEFD <- state["AESEFD"]
    AESEFED <- state["AESEFED"]
    AESEFEO <- state["AESEFEO"]
    AESEFOR <- state["AESEFOR"]
    AESEFOD <- state["AESEFOD"]
    AESEFOE <- state["AESEFOE"]
    AESEFOEO <- state["AESEFOEO"]
    AESEFEOR <- state["AESEFEOR"]
    AESEFEOD <- state["AESEFEOD"]
    AESEFOED <- state["AESEFOED"]
    AESEFOEOR <- state["AESEFOEOR"]
    AESEFOEOD <- state["AESEFOEOD"]
    AESES <- state["AESES"]
    AESESD <- state["AESESD"]
    AESESO <- state["AESESO"]
    AESESOR <- state["AESESOR"]
    AESESOD <- state["AESESOD"]
    
    BEF <- state["BEF"]
    BEFE <- state["BEFE"]
    BEFO <- state["BEFO"]
    BEFD <- state["BEFD"]
    BEFED <- state["BEFED"]
    BEFEO <- state["BEFEO"]
    BEFOR <- state["BEFOR"]
    BEFOD <- state["BEFOD"]
    BEFOE <- state["BEFOE"]
    BEFEOR <- state["BEFEOR"]
    BEFEOD <- state["BEFEOD"]
    BEFOED <- state["BEFOED"]
    BEFOEO <- state["BEFOEO"]
    BEFOEOR <- state["BEFOEOR"]
    BEFOEOD <- state["BEFOEOD"]
    
    CES <- state["CES"]
    CESD <- state["CESD"]
    CESO <- state["CESO"]
    CESOR <- state["CESOR"]
    CESOD <- state["CESOD"]
    
    
    #Equations
    dAES = NS*h*WA - AES*(e1*D1 + d + a1 + a2)
    dAESD = AES*e1*D1
    dAESO = AES*d - AESO*(e1*D1 + R1 + a4 + a5)
    dAESOR = AESO*R1
    dAESOD = AESO*e1*D1
    dAESOEF = AESO*a4 - AESOEF*(a6 + d + e2*D2)
    dAESOES = AESO*a5 - AESOES*(d + e2*D3)
    dAESOEFE = AESOEF*a6 - AESOEFE*(d + e3*D3)
    dAESOEFO = AESOEF*d - AESOEFO*(R2 + e2*D2 + a7)
    dAESOEFD = AESOEF*e2*D2
    dAESOESD = AESOES*e2*D3
    dAESOESO = AESOES*d - AESOESO*(R2 + e2*D3)
    dAESOEFED = AESOEFE*e3*D3
    dAESOEFEO = AESOEFE*d - AESOEFEO*(R2 + e3*D3)
    dAESOEFOR = AESOEFO*R2
    dAESOEFOD = AESOEFO*e2*D2
    dAESOEFOE = AESOEFO*a7 - AESOEFOE*(d + e3*D3)
    dAESOESOR = AESOESO*R2
    dAESOESOD = AESOESO*e2*D3
    dAESOEFEOR = AESOEFEO*R2
    dAESOEFEOD = AESOEFEO*e3*D3
    dAESOEFOED = AESOEFOE*e3*D3
    dAESOEFOEO = AESOEFOE*d - AESOEFOEO*(R2 + e3*D3)
    dAESOEFOEOR = AESOEFOEO*R2
    dAESOEFOEOD = AESOEFOEO*e3*D3
    
    dAESEF = AES*a1 - AESEF*(a6 + d + e4*D2)
    dAESEFE = AESEF*a6 - AESEFE*(e5*D3 + d)
    dAESEFO = AESEF*d - AESEFO*(a7 + R3 + e4*D2)
    dAESEFD = AESEF*e4*D2
    dAESEFED = AESEFE*e5*D3
    dAESEFEO = AESEFE*d - AESEFEO*(R2 + e5*D3)
    dAESEFOR = AESEFO*R3
    dAESEFOD = AESEFO*e4*D2
    dAESEFOE = AESEFO*a7 - AESEFOE*(e5*D3 + d)
    dAESEFOEO = AESEFOE*d - AESEFOEO*(R2 + e5*D3)
    dAESEFEOR = AESEFEO*R2
    dAESEFEOD = AESEFEO*e5*D3
    dAESEFOED = AESEFOE*e5*D3
    dAESEFOEOR = AESEFOEO*R2
    dAESEFOEOD = AESEFOEO*e5*D3
    
    dAESES = AES*a2 - AESES*(e4*D3 + d)
    dAESESD = AESES*e4*D3
    dAESESO = AESES*d - AESESO*(R4 + e4*D3)
    dAESESOR = AESESO*R4
    dAESESOD = AESESO*e4*D3
    
    dBEF = NS*h*WB - BEF*(a3 + d + e6*D2)
    dBEFE = BEF*a3 - BEFE*(d + e7*D3)
    dBEFO = BEF*d - BEFO*(R5 + a7 + e6*D2)
    dBEFD = BEF*e6*D2
    dBEFED = BEFE*e7*D3
    dBEFEO = BEFE*d - BEFEO*(R2 + e7*D3)
    dBEFOR = BEFO*R5
    dBEFOD = BEFO*e6*D2
    dBEFOE = BEFO*a7 - BEFOE*(d + e7*D3)
    dBEFEOR = BEFEO*R2
    dBEFEOD = BEFEO*e7*D3
    dBEFOED = BEFOE*e7*D3
    dBEFOEO = BEFOE*d - BEFOEO*(R2 + e7*D3)
    dBEFOEOR = BEFOEO*R2
    dBEFOEOD = BEFOEO*e7*D3
    
    dCES = NS*h*WC - CES*(d + e8*D3)
    dCESD = CES*e8*D3
    dCESO = CES*d - CESO*(R6 + e8*D3)
    dCESOR = CESO*R6
    dCESOD = CESO*e8*D3
    
    #Outcome equations - NEEDS WORK!
    dh = h*L
    
    dND = CESO*e8*D3 +  CES*e8*D3 + BEFOEO*e7*D3 + BEFOE*e7*D3 + BEFEO*e7*D3 + BEF*e6*D2 + BEFE*e7*D3 + BEFO*e6*D2 + AESESO*e4*D3 + AESEFE*e5*D3 + AESEFO*e4*D2 + AESEFEO*e5*D3 + AESEFOE*e5*D3 + AESEFOEO*e5*D3 + AESES*e4*D3 + AES*e1*D1 + AESO*e1*D1 + AESOEF*e2*D2 + AESOES*e2*D3 + AESOEFE*e3*D3 + AESOEFO*e2*D2 + AESOESO*e2*D3 + AESOEFEO*e3*D3 + AESOEFOE*e3*D3 + AESOEFOEO*e3*D3
    
    dNR = AESO*R1 + AESOEFO*R2 + AESOESO*R2 + AESOEFEO*R2 + AESOEFOEO*R2 + AESEFO*R3 + AESEFEO*R2 + AESEFOEO*R2 + AESESO*R4 + BEFO*R5 + BEFEO*R2 + BEFOEO*R2 + CESO*R6
    
    dTA = AES*d - AESO*(e1*D1 + R1 + a4 + a5) + AES*a1 - AESEF*(a6 + d + e4*D2) + AES*a2 - AESES*(e4*D3 + d) + AESO*a4 - AESOEF*(a6 + d + e2*D2) + AESO*a5 - AESOES*(d + e2*D3) + SA*(NS*h*WA - AES*(e1*D1 + d +a1 +a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3))
    
    dTB = AESOEF*a6 - AESOEFE*(d + e3*D3) + AESOEF*d - AESOEFO*(R2 + e2*D2 + a7) + AESOEFO*a7 - AESOEFOE*(d + e3*D3) + AESEF*a6 - AESEFE*(e5*D3 + d) + AESEF*d - AESEFO*(a7 + R3 + e4*D2) + AESEFO*a7 - AESEFOE*(e5*D3 + d) + BEF*a3 - BEFE*(d + e7*D3) + BEF*d - BEFO*(R5 +a7 + e6*D2) + BEFO*a7 - BEFOE*(d + e7*D3) + SB*(NS*h*WA - AES*(e1*D1 + d + a1 + a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3))
    
    dTC = AESOEFE*d - AESOEFEO*(R2 + e3*D3) + AESOEFOE*d - AESOEFOEO*(R2 + e3*D3) + AESOES*d - AESOESO*(R2 + e2*D3) + AESEFE*d - AESEFEO*(R2 + e5*D3) + AESEFOE*d - AESEFOEO*(R2 + e5*D3) + AESES*d - AESESO*(R4 + e4*D3) + BEFE*d - BEFEO*(R2 + e7*D3) + BEFOE*d - BEFOEO*(R2 + e7*D3) + CES*d - CESO*(R6 + e8*D3) + SC*(NS*h*WA - AES*(e1*D1 + d +a1 +a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3))
    
    dNE = NS*h*WA - AES*(e1*D1 + d +a1 +a2) + NS*h*WB - BEF*(a3 + d + e6*D2) + NS*h*WC - CES*(d + e8*D3)
    
    dWA = -1*gamma*TA - (1/3)*gamma*TB - (1/3)*gamma*TC
    dWB = 1*gamma*TA - (1/3)*gamma*TB - (2/3)*gamma*TC
    dWC = (2/3)*gamma*TB + 1*gamma*TC
    
    dTAA = 1*TA
    dTBB = 1*TB
    dTCC = 1*TC
    
    
    dxdt <- c(dh, dTAA, dTBB, dTCC, dND, dNR, dTA, dTB, dTC, dNE, dWA, dWB, dWC, dAES, dAESD, dAESO, dAESOR, dAESOD, dAESOEF, dAESOES, dAESOEFE, dAESOEFO, dAESOEFD, dAESOESD, dAESOESO, dAESOEFED, dAESOEFEO, dAESOEFO, dAESOEFOD, dAESOEFOE, dAESOESOR, dAESOESOD, dAESOEFEOR, dAESOEFEOD, dAESOEFOED, dAESOEFOEO, dAESOEFOEOR, dAESOEFOEOD, dAESEF, dAESEFE, dAESEFO, dAESEFD, dAESEFED, dAESEFEO, dAESEFOR, dAESEFOD, dAESEFOE, dAESEFOEO, dAESEFEOR, dAESEFEOD, dAESEFOED, dAESEFOEOR, dAESEFOEOD, dAESES, dAESESD, dAESESO, dAESESOR, dAESESOD, dBEF, dBEFE, dBEFO, dBEFD, dBEFED, dBEFEO, dBEFOR, dBEFOD,dBEFOE, dBEFEOR, dBEFEOD, dBEFOED, dBEFOEO, dBEFOEOR, dBEFOEOD, dCES, dCESD, dCESO, dCESOR, dCESOD)
    ## return result as a list
    list(dxdt)}
  )
}

#Scenario A, Baseline
parms <- c(L=0.000066, t=7, a1=0.001, a2=0.001/5, a3=0.001, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^2, e4=(1.3^((1/0.5)/2)), e5=(1.3^((1/0.5)/2))^2, e6=(1.3^((1/0.5)/2)), e7=(1.3^((1/0.5)/2))^2, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/7), R4=(1/7), R5=(1/7), R6=(1/7), SA=1, SB=0, SC=0)
state <- c(h=0.00000212, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
library(tidyverse)
library(deSolve)
outputBaselineA<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))
view(outputBaselineA)

#Add new column
outputBaselineA$cum_mort_rate <- outputBaselineA$ND/(outputBaselineA$NR + outputBaselineA$ND)*100
view(outputBaselineA)
outputBaselineA$Scenario <- "A"

#Graphs
options(scipen=999)  # turn off scientific notation like 1e+06
library(ggplot2)

#Scenario B, Baseline
parms <- c(L=0.000066, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^1, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^1, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^1, e8=(1.3^((1/0.5)/2)), gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/7), R5=(1/(7-(1/0.5))), R6=(1/7), SA=0, SB=1, SC=0)
state <- c(h=0.00000212, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
library(tidyverse)
library(deSolve)
outputBaselineB<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))
view(outputBaselineB)

#Add new column
outputBaselineB$cum_mort_rate <- outputBaselineB$ND/(outputBaselineB$NR + outputBaselineB$ND)*100
view(outputBaselineB)
outputBaselineB$Scenario <- "B"


#Scenario C, Baseline
parms <- c(L=0.000066, t=7, a1=0.001*1.375, a2=0.001*1.125, a3=0.001*1.5, a4=0.001, a5=0.001/5, a6=0.001, a7=0.001*1.5, d=0.5, D=0.01, D1=0.01, D2=0.01*1.5, D3=0.01*2, NS=100000, e=1.3^((1/0.5)/2), e1=(1.3^((1/0.5)/2))^0, e2=(1.3^((1/0.5)/2))^0, e3=(1.3^((1/0.5)/2))^0, e4=(1.3^((1/0.5)/2))^0, e5=(1.3^((1/0.5)/2))^0, e6=(1.3^((1/0.5)/2))^0, e7=(1.3^((1/0.5)/2))^0, e8=(1.3^((1/0.5)/2))^0, gamma=0.000018, R1=(1/(7-(1/0.5))), R2=1/7, R3=(1/(7-(1/0.5))), R4=(1/(7-(1/0.5))), R5=(1/(7-(1/0.5))), R6=(1/(7-(1/0.5))), SA=0, SB=0, SC=1)
state <- c(h=0.00000212, TAA=0, TBB=0, TCC=0, ND=0, NR=0, TA=0, TB=0, TC=0, NE=0, WA= 0.914, WB=0.036, WC=0.050, AES= 0, AESD= 0, AESO= 0, AESOR= 0, AESOD= 0, AESOEF= 0, AESOES= 0, AESOEFE= 0, AESOEFO= 0, AESOEFD= 0, AESOESD= 0, AESOESO= 0, AESOEFED= 0, AESOEFEO= 0, AESOEFOR= 0, AESOEFOD= 0, AESOEFOE= 0, AESOESOR= 0, AESOESOD= 0, AESOEFEOR= 0, AESOEFEOD= 0, AESOEFOED= 0, AESOEFOEO= 0, AESOEFOEOR= 0, AESOEFOEOD= 0, AESEF= 0, AESEFE= 0, AESEFO= 0, AESEFD= 0, AESEFED= 0, AESEFEO= 0, AESEFOR= 0, AESEFOD= 0, AESEFOE= 0, AESEFOEO= 0, AESEFEOR= 0, AESEFEOD= 0, AESEFOED= 0, AESEFOEOR= 0, AESEFOEOD= 0, AESES= 0, AESESD= 0, AESESO= 0, AESESOR= 0, AESESOD= 0, BEF= 0, BEFE= 0, BEFO= 0, BEFD= 0, BEFED= 0, BEFEO= 0, BEFOR= 0, BEFOD= 0, BEFOE= 0, BEFEOR= 0, BEFEOD= 0, BEFOED= 0, BEFOEO= 0, BEFOEOR= 0, BEFOEOD= 0, CES= 0, CESD= 0, CESO= 0, CESOR= 0, CESOD= 0)
times <- seq(from=0,to=365*5, by=25)
library(tidyverse)
library(deSolve)
outputBaselineC<-as.data.frame(ode(func=new.basic.model, y=state, times=times, parms=parms))
view(outputBaselineC)

#Add new column
outputBaselineC$cum_mort_rate <- outputBaselineC$ND/(outputBaselineC$NR + outputBaselineC$ND)*100
view(outputBaselineC)
outputBaselineC$Scenario <- "C"

#Combining data frames
DatasetBaseline <- rbind(outputBaselineA, outputBaselineB, outputBaselineC)
view(DatasetBaseline)

#GridPlot Look

#ND
Deathplot <- ggplot(data = DatasetBaseline, aes(time/365, ND)) +
  geom_line(color = "seagreen3", size = 0.1) +
  geom_point(color="seagreen3", size = 0.5) + 
  labs(title = "Total # of Deaths",
       y = "# Deaths", x = "Time (Years)") +
  theme(plot.title = element_text(size = 9)) +
  facet_wrap(~ Scenario)
Deathplot

#MRate
MRateplot <- ggplot(data = DatasetBaseline, aes(time/365, cum_mort_rate)) +
  geom_line(color = "seagreen3", size = 0.1) +
  geom_point(color="seagreen3", size = 0.5) + 
  labs(title = "Time-updated Case Fatality Rate (%)",
       y = "%", x = "") +
  coord_cartesian(xlim=c(1, 5), ylim=c(7, 8)) +
  theme(plot.title = element_text(size = 9)) +
  facet_wrap(~ Scenario)

#Res 1st line
Res1plot <- ggplot(data = DatasetBaseline, aes(time/365, WB*100)) +
  geom_line(color = "seagreen3", size = 0.5) +
  geom_point(color="seagreen3", size = 0.5) + 
  labs(title = "% Resistance, 1st Line Therapy",
       y = "%", x = "") +
  theme(plot.title = element_text(size = 9)) +
  facet_wrap(~ Scenario)

#Res 2nd line
Res2plot <- ggplot(data = DatasetBaseline, aes(time/365, WC*100)) +
  geom_line(color = "seagreen3", size = 0.5) +
  geom_point(color="seagreen3", size = 0.5) + 
  labs(title = "% Resistance, 2nd Line Therapy",
       y = "%", x = "") +
  theme(plot.title = element_text(size = 9)) +
  facet_wrap(~ Scenario)

#Days 1st Line AB
firstlineplot <- ggplot(data = DatasetBaseline, aes(time/365, TAA)) +
  geom_line(color = "seagreen3", size = 0.5) +
  geom_point(color="seagreen3", size = 0.5) + 
  labs(title = "# Days Treatment, 1st Line Antibiotics",
       y = "Days", x = "") +
  theme(plot.title = element_text(size = 9)) +
  facet_wrap(~ Scenario)

#Days 2nd Line AB
secondlineplot <- ggplot(data = DatasetBaseline, aes(time/365, TBB)) +
  geom_line(color = "seagreen3", size = 0.5) +
  geom_point(color="seagreen3", size = 0.5) + 
  labs(title = "# Days Treatment, 2nd Line Antibiotics",
       y = "Days", x = "") +
  theme(plot.title = element_text(size = 9)) +
  facet_wrap(~ Scenario)

#Days 3rd Line AB
thirdlineplot <- ggplot(data = DatasetBaseline, aes(time/365, TCC)) +
  geom_line(color = "seagreen3", size = 0.5) +
  geom_point(color="seagreen3", size = 0.5) + 
  labs(title = "# Days Treatment, 3rd Line Antibiotics",
       y = "Days", x = "") +
  theme(plot.title = element_text(size = 9)) +
  facet_wrap(~ Scenario)


#Combining into one
library("cowplot")
BaselineCombined <- ggdraw() +
  draw_plot(Deathplot, x = 0.02, y = 6/7, width = .95, height = 1/7) +
  draw_plot(MRateplot, x = 0.02, y = 5/7, width = .95, height = 1/7) +
  draw_plot(Res1plot, x = 0.02, y = 4/7, width = .95, height = 1/7) +
  draw_plot(Res2plot, x = 0.02, y = 3/7, width = .95, height = 1/7) +
  draw_plot(firstlineplot, x = 0.02, y = 2/7, width = .95, height = 1/7) +
  draw_plot(secondlineplot, x =0.02, y = 1/7, width = .95, height = 1/7) +
  draw_plot(thirdlineplot, x = 0.02, y = 0/7, width = .95, height = 1/7)

BaselineCombined





