# Lab 11: vector-borne infections

library(deSolve)

############################################################
## PART 1 Model Dynamics
############################################################

# CODE UP THE MODEL ON YOUR OWN
SIRMosVec = function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    dSH = vH - r*THM*SH*IM - muH*SH
    dIH = r*THM*SH*IM - gamma*IH - muH*IH
    
    dSM = vM - r*TMH*SM*IH - muM*SM
    dIM = r*TMH*SM*IH - muM*IM
    
    list(c(dSH, dIH, dSM, dIM))
  })
}

# initial states and parameters
NH=1e7; # human population size
IH = 1; # initial number of infection in humans
SH=NH-IH; # initial susceptible people
muH=1/50/365; # human mortality rate in days, human life span: 50 yr
vH=NH*muH; # number of newborns per day
TMH = 0.2; # prob infection from human to mosquito;
THM = 0.1; # prob infection from mosquito to human;
gamma=1/7; # infectious period: 7 days
NM=1e8; # mosquito population size
IM=1; # initial number of infection in mosquitoes
SM=NM-IM; # initial susceptible mosquitoes
muM = 1/7; # mosquito mortality rate, 1 week life span for mosquito
vM=NM*muM; # births in mosquito population
b=.5; # number of bite per mosquito per day;
r = b / NH; # bite rate per human per mosquito

parameters = c(muH = muH, muM = muM,
               vH = vH, vM = vM, 
               THM = THM, TMH = TMH, 
               gamma = gamma, r = b / NH)

state = c(SH = SH, IH = 1, SM = SM, IM = 1)

times=1:(365*500);

sim = ode(y = state, times = times, func = SIRMosVec, parms = parameters)

# [LQ1] Set up the code for the mosquito-borne disease model using the questions in Slide 5 and Run the model using parameters/initial conditions in Slide 6, for 500 years. (1pt)
# [LQ2] Based on the model-simulated results, what are the peak (i.e. maximum) prevalence values for human and mosquito? (0.5pt)
max(sim[ , 'IH'])
max(sim[ , 'IM'])

# [LQ3] Based on the model-simulated results, what is the prevalence of infection for human at equilibrium?
tail(sim[ , 'IH'])
tail(sim[ , 'IM'])

# [LQ4] Why the prevalence of human infection at equilibrium is much lower than that at the peak? (0.5pt)
# Hint: to understand this, plot the % immune to infection in humans (i.e. the ‘R’ class for humans, or NH – SH – IH) vs. time
RH = NH - sim[ , 'SH'] - sim[ , 'IH']

plot(x = sim[, 'time'], y = RH, )

############################################################
## PART 2 Factors shaping the model Dynamics
############################################################
## Calculating R0
## check the equatioin for R0 in the lecture slides
R0 = (b*THM)/(muM) * (b*TMH*NM)/((gamma+muH)*NH)

tail(sim[ , 'SH'], 1)/NH 
#####################
## Test NM/NH vs R0 vs epidemic dynamics
## use the R Shiny App: ShinyApp_VectorBorneDis.R
