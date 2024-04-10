#### R files for HW1 ####
#### Date: Apr 7, 2024 ####
#### Author: Ting Ye ####

# If any of the R packages is not previously installed, use install.packages("tidyverse") and install.packages("MASS") to install them first. 
# We will also use the `tidyverse` package for data manipulation, 
# and `MASS` for simulating negative binomial outcomes and fitting negative binomial GLMs
library(tidyverse)
library(MASS)

# Install the RobinCar package which is an integrated package for covariate adjustment
# It may require installing the `devtools` package using install.packages("devtools")
devtools::install_github("mbannick/RobinCar")
library(RobinCar)

# This is a function we will use for data generation. A=1 is treated and A=0 is control
Fun_datagen <- function(Fun.n = 500, 
                        Fun.y.type = c("continuous","binary","count"), 
                        Fun.p = 2/3){
  
  # generate the covariates Xc, Xb
  df <- tibble(
    Xc = runif(Fun.n), 
    Xb = rbinom(Fun.n, size = 1, prob = 0.5)
  )
  
  # generate treatment assignment A by simple randomization
  df <- df %>% mutate(A = rbinom(n=Fun.n, size=1, prob=Fun.p))
  
  # generate outcome y
  if(Fun.y.type == "continuous"){
    df <- df %>%
      mutate(y = (1-A)*(Xc+0.1*Xb) + A*(0.3*Xc+0.3*Xb))
  }else if(Fun.y.type == "binary"){
    df <- df %>% mutate(y = rbinom(n = Fun.n, size = 1, prob = (1-A)*(0.5*Xc+0.25*Xc^2+0.1*Xb) + A*exp(Xc+Xc^2+0.3*Xb)/(1+exp(Xc+Xc^2+0.3*Xb))))
  }else if(Fun.y.type == "count"){
    df <- df %>% mutate(y = MASS::rnegbin(n = Fun.n, mu = A*(2*Xc+5*Xc^2+0.1*Xb) + (1-A)*log(6*Xc^3+2+0.3*Xb), theta = 4))
  }
  
  df <- df %>% mutate(A = factor(A))
  return(df)
}
