library(tidyverse)
library(readxl)

## Expenditure Side ----
Exp <- read.csv("./DataRaw/GDP-Exp.csv")[,-1] %>% as.matrix()
Exp <- ts(Exp, start=c(1984,4), frequency=4)

nb_e <- 53 #Number of most disaggregate series
n_e <- ncol(Exp) #Total number of series in the hierarchy
na_e <- n_e - nb_e

C_EXP <- read.csv("./DataRaw/S_mat.csv")[1:na_e,-1] %>% as.matrix()

## Income Side ----
MDI <- read_excel("./DataRaw/Master Data File.xlsx", sheet=5, skip = 9) #Master Data for Income

#a_t = Gdpi, Tfi, TfiGos, TfiCoe, TfiGosCop, TfiGosCopNfn,
#b_t = TfiGosCopNfnPub, TfiGosCopNfnPvt, TfiGosCopFin,TfiGosGvt,TfiGosDwl,TfiGmi,TfiCoeWns,TfiCoeEsc,Tsi,Sdi

Inc <- tibble(Gdpi = MDI %>% pull("A2302467A") %>% ts(start=c(1959,3), frequency=4) %>%
                window(start=c(1984,4))) #GDP(I)

#Total factor income (GOS + GMI)
Inc %>% add_column(Tfi = MDI %>% pull("A2302411R") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Total factor income

#All sectors GOS (Corporations + Gen. Govn. + Dwellings)
Inc %>% add_column(TfiGos = MDI %>% pull("A2302409C") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #All sectors ;  Gross operating surplus

#Compensation of employees (Wages and salaries + Social contribution)
Inc %>% add_column(TfiCoe = MDI %>% pull("A2302401K") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Compensation of employees

#Total corporations (Non-financial(Pub + Pvt) + Financial)
Inc %>% add_column(TfiGosCop = MDI %>% pull("A2302406W") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Total corporations ;  Gross operating surplus

#Non-financial corporations (Public + Private)
Inc %>% add_column(TfiGosCopNfn = MDI %>% pull("A2302404T") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Non-financial corporations ;  Gross operating surplus

#Public non-financial corporations
Inc %>% add_column(TfiGosCopNfnPub = MDI %>% pull("A2302403R") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Public non-financial corporations ;  Gross operating surplus

#Private non-financial corporations
Inc %>% add_column(TfiGosCopNfnPvt = MDI %>% pull("A2323369L") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Private non-financial corporations ;  Gross operating surplus

#Financial corporations
Inc %>% add_column(TfiGosCopFin = MDI %>% pull("A2302405V") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Financial corporations ;  Gross operating surplus


#General government
Inc %>% add_column(TfiGosGvt = MDI %>% pull("A2298711F") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #General government ;  Gross operating surplus

#Dwellings
Inc %>% add_column(TfiGosDwl = MDI %>% pull("A2302408A") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Dwellings owned by persons ;  Gross operating surplus

#Gross mixed income
Inc %>% add_column(TfiGmi = MDI %>% pull("A2302410L") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Gross mixed income

#Compensation of employees - Wages and salaries
Inc %>% add_column(TfiCoeWns = MDI %>% pull("A2302399K") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc 

#Compensation of employees - Employers' social contributions
Inc %>% add_column(TfiCoeEsc = MDI %>% pull("A2302400J") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc #Compensation of employees - Employers' social contributions

#Taxes less subsidies (I)
Inc %>% add_column(Tsi = MDI %>% pull("A2302412T") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc 

#Statistical Discrepancy (I)
Inc %>% add_column(Sdi = MDI %>% pull("A2302413V") %>% ts(start=c(1959,3), frequency=4) %>% 
                     window(start=c(1984,4))) -> Inc 
Inc <- ts(Inc, start=c(1984,4), frequency=4)

nb_i <- 10L       # Number of most disaggregate series
n_i <- ncol(Inc) # Total number of series in the hierarchy
na_i <- n_i - nb_i
# S: Summing matrix
C_INC <- matrix(c(1,1,1,1,1,1,1,1,1,1,
                  1,1,1,1,1,1,1,1,0,0,
                  1,1,1,1,1,0,0,0,0,0,
                  0,0,0,0,0,0,1,1,0,0,
                  1,1,1,0,0,0,0,0,0,0,
                  1,1,0,0,0,0,0,0,0,0), nrow=nb_i) %>% t


## Complete mts ----
AusGDP <- bind_cols(Inc[,1:na_i],Exp[,2:na_e],Inc[,-(1:na_i)],Exp[,-(1:na_e)]) |>
  rename(Gdp = Gdpi)
AusGDP <- ts(AusGDP, start=c(1984,4), frequency=4)

# Number of aggregate level TS
na <- na_i + na_e - 1
# Number of bottom level TS
nb <- nb_i + nb_e
# Number of all TS
n <- na + nb

M1 <- matrix(c(1,rep(0,(na_i-1)),rep(0,(na_e-1)),-rep(1,nb_i), rep(0,nb_e),
               1,rep(0,(na_i-1)),rep(0,(na_e-1)), rep(0,nb_i),-rep(1,nb_e)), nrow=2, byrow=T)
A11 <- matrix(c(rep(0,(na_i-1))),nrow=(na_i-1))
A12 <- diag(na_i-1)
A13 <- matrix(c(rep(0,(na_i-1)*(na_e-1))),nrow=(na_i-1))
A15 <- matrix(c(rep(0,(na_i-1)*(nb_e))),nrow=(na_i-1))
A <- cbind(A11,A12,A13,-C_INC[-1,],A15)

B11 <- matrix(c(rep(0,(na_e-1))),nrow=(na_e-1))
B12 <- matrix(c(rep(0,(na_e-1)*(na_i-1))),nrow=(na_e-1))
B13 <- diag(na_e-1)
B14 <- matrix(c(rep(0,(na_e-1)*(nb_i))),nrow=(na_e-1))
B <- cbind(B11,B12,B13,B14,-C_EXP[-1,])

Ut <- rbind(M1,A,B)
colnames(Ut) <- colnames(AusGDP)

library(FoReco)
extract_C <- FoReco::lcmat(Gt = Ut)

C <- extract_C$Cbar
nb <- NCOL(C)
na <- NROW(C)
n <- na + nb
DATA <- cbind(as.matrix(AusGDP[, -c(1:na)]%*%t(C)), AusGDP[, -c(1:na)])
colnames(DATA) <- colnames(AusGDP)

K <- c(1,2,4)
H <- 4
m <- 4
fixed_length <- 40
save(DATA, C, K, H, m, fixed_length, file = "DATA.RData")

