setwd("C:/Users/Fran/Dropbox/Pesquisa/repairable_systems")
setwd("~/Dropbox/Pesquisa/repairable_systems")


# =============================================================================
# data set and some objects
db <- read.csv("data/Gilardoni2007.csv", sep = " ", header = FALSE)
cost          <- db[1, 1]    # CMR/CPM
db            <- db[-1, ]    # data set
nsystem       <- nrow(db)  
rownames(db)  <- 1:nsystem
nfailureTimes <- ncol(db)
failureTimes  <- matrix(0, ncol = nfailureTimes, nrow = 1)
censorTime    <- 0


# =============================================================================
# failures and censoring times
for (i in 1:nsystem){
    aux1 <- sum(is.na(db[i, ]) == FALSE)
    if(aux1 == 1){
        failureTimes   <- rbind(failureTimes, rep(0,nfailureTimes))     
        censorTime[i]  <- db[i, 1]
    } else {
        failureTimes   <- rbind(failureTimes, as.numeric(c(db[i, 1:(aux1-1)],
                                rep(0,(nfailureTimes-aux1+1)))))     
        censorTime[i]  <- db[i,aux1]
    }
}
failureTimes <- failureTimes[-1, ]
failureTimes # 0 == censuring
censorTime


# =============================================================================
# mean cumulative number of failures 
t1   <- sort(unique(c(0, failureTimes, censorTime)))
qq   <- 0
SCNF <- 0
MCNF <- 0

for(i in seq_along(t1)){
    qq[i]  <- sum(censorTime > t1[i])  # number of uncensored system until time t1[i]
    if (sum(censorTime <= t1[i]) > 0){     
        aux     <- failureTimes[- which(censorTime <= t1[i]), ]
        aux     <- c(aux)
        aux     <- aux[aux > 0]
        SCNF[i] <- sum(aux <= t1[i])   # mean cumulative number of failures until time t1[i]
    }
}
MCNF <- SCNF / qq

plot(t1, MCNF, type = "s", xlab = "time (hours)", 
     ylab = " mean cumulative number of failures")



# =============================================================================
# mean cumulative number of failures 
t1   <- sort(unique(c(0, failureTimes, censorTime)))
qq   <- 0
SCNF <- 0
MCNF <- 0

for(i in seq_along(t1)){
    qq[i]  <- sum(censorTime > t1[i])  # number of uncensored system until time t1[i]
    if (sum(censorTime <= t1[i]) > 0){     
        aux     <- failureTimes[- which(censorTime <= t1[i]), ]
        aux     <- c(aux)
        aux     <- aux[aux > 0]
        SCNF[i] <- sum(aux <= t1[i])   # mean cumulative number of failures until time t1[i]
    }
}
MCNF <- SCNF / qq

plot(t1, MCNF, type = "s", xlab = "time (hours)", 
     ylab = " mean cumulative number of failures")




# ====================================================================
# Returns the sorted time of each failure considering all systems.
require(lattice)
dotplot(as.matrix(db), main = "Failures and censorings times", 
        xlab = "time (hours)", ylab = "system", pch = 19)












# for (i in 1:ndata){
#     aux1 <- sum(is.na(db[i, ]) == FALSE)
#     if(aux1 == 1){
#         failureTimes    <- rbind(failureTimes, c(db[i, 1], rep(0,(nfailureTimes-aux1))))     
#         censorTime[i]   <- db[i, 1]
#     } else {
#           aux2          <- diff(as.numeric(db[i,1:aux1 ]), lag = 1)
#           failureTimes  <- rbind(failureTimes, c(db[i,1], aux2, rep(0,nfailureTimes-aux2+1)))
#           censorTime[i] <- db[i, aux1]
#     }
# }




