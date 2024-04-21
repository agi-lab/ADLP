
## Define the offsets values used for Gamma and Log-Normal models:
tau <- 5
tau_LN <- tau
tau_Ga <- tau

## Define the component density and CDF functions

dODP<-function(y,lambda,phi) {
    z <- ifelse(rep(phi, length(y)) < 1e-6, dpois(floor(y), lambda), dpois(floor(y/phi), floor(lambda/phi))/phi)
    return (z)
}

pODP<-function(y,lambda,phi) {
    lambda <- as.vector(lambda)
    phi <- as.vector(phi)
    return (
        sapply(
            1:length(y),
            function (x)
                tweedie::ptweedie(q = y[x],
                                  mu = lambda[(x-1) %% length(lambda) + 1],
                                  phi = phi[(x-1) %% length(phi) + 1],
                                  power = 1)
        )
    )
}

dens_ODP_GLM <- function(y, model, newdata) {
    pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
    pred_ODP_mu<-pred_ODP$fit
    pred_ODP_phi<-(pred_ODP$residual.scale)^2
    return(dODP(y,lambda=pred_ODP_mu,phi=pred_ODP_phi))
}

CDF_ODP_GLM<-function(y, model, newdata){
    pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
    pred_ODP_mu<-pred_ODP$fit
    pred_ODP_phi<-(pred_ODP$residual.scale)^2
    return(pODP(y,lambda=pred_ODP_mu,phi=pred_ODP_phi))
}

mu_ODP_GLM<-function(model, newdata){
    pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
    pred_ODP_mu<-pred_ODP$fit
    return(pred_ODP_mu)
}

sim_ODP_GLM<-function(model, newdata){
    pred_ODP<-predict.glm(model,newdata=newdata,type="response",se.fit=TRUE)
    pred_ODP_mu<-pred_ODP$fit
    pred_ODP_phi<-(pred_ODP$residual.scale)^2
    simy<-replicate(1,rtweedie(length(pred_ODP_mu),xi=1,mu=pred_ODP_mu,phi=pred_ODP_phi))
    return(simy)
}

dens_LN <- function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(dLNO(x=y,mu=pred_mu,sigma=pred_sigma))
}

CDF_LN<-function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(pLNO(q=y+tau,mu=pred_mu,sigma=pred_sigma))
}

mu_LN<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    mean<-exp(pred_mu+pred_sigma^2/2)
    return(mean - tau)
}

sim_LN<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    simy<-replicate(1, rLNO(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
    simy[simy<0]<-0
    return(simy)
}

dens_GA <- function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(dGA(x=y,mu=pred_mu,sigma=pred_sigma))
}

CDF_GA<-function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(pGA(q=y+tau,mu=pred_mu,sigma=pred_sigma))
}

mu_GA<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
    }))
    return(pred_mu - tau)
}

sim_GA<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    simy<-replicate(1, rGA(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
    simy[simy<0]<-0
    return(simy)
}

dens_NO <- function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(dNO(x=y,mu=pred_mu,sigma=pred_sigma))
}

CDF_Normal<-function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(pNO(q=y,mu=pred_mu,sigma=pred_sigma))
}

mu_Normal<-function(model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
    }))
    return(pred_mu)
}

sim_NO<-function(model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    simy<-replicate(1,rNO(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma))
    simy[simy<0]<-0
    return(simy)
}

dens_ZAGA<-function(y, model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    pred_nu<-predict(model,what="nu",newdata=newdata_nu,type="response")
    return(dZAGA(x=y,mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
}

CDF_ZAGA<-function(y, model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    pred_nu<-predict(model,what="nu",newdata=newdata_nu,type="response")
    return(gamlss.dist::pZAGA(q=y,mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
}

mu_ZAGA<-function(model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
    pred_nu<-predict(model,what="nu",newdata=newdata_nu,type="response")
    mean<-(1-pred_nu)*pred_mu
    return(mean)
}

sim_ZAGA<-function(model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
    pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    pred_nu<-predict(model,what="nu",newdata=newdata_nu,type="response")
    pred_nu<-mean(pred_nu)
    simy<-replicate(1, rZAGA(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma,nu=pred_nu))
    return(simy)
}

dZALN<-gamlss.inf::Zadj.d(family="LOGNO")
pZALN<-gamlss.inf::Zadj.p(family="LOGNO")
rZALN<-gamlss.inf::Zadj.r(family="LOGNO")
dens_ZALN<-function(y, model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
    pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
    pred_nu<-predict(model,parameter="xi0",newdata=newdata_nu,type="response")
    return(dZALN(x=y,mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
}

CDF_ZALN<-function(y, model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
    pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
    pred_nu<-predict(model,parameter="xi0",newdata=newdata_nu,type="response")
    return(pZALN(q=y,mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
}

mu_ZALN<-function(model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
    pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
    pred_nu<-predict(model,parameter="xi0",newdata=newdata_nu,type="response")
    mean<-(1-pred_nu)*exp(pred_mu+pred_sigma^2/2)
    return(mean)
}

sim_ZALN<-function(model, newdata){
    newdata_nu<-newdata
    newdata_nu$dev=as.numeric(as.character(newdata$dev))
    pred_mu<-predict(model,parameter="mu",newdata=newdata,type="response")
    pred_sigma<-predict(model,parameter="sigma",newdata=newdata,type="response")
    pred_nu<-predict(model,parameter="xi0",newdata=newdata_nu,type="response")
    pred_nu<-mean(pred_nu)
    simy<-replicate(1,rZALN(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma,xi0=pred_nu))
    return(simy)
}

##Calculate Gamma density for Gamma GAMLSS
dens_GA_Gamlss<-function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(dGA(x=y,mu=pred_mu,sigma=pred_sigma))
}

CDF_GA_Gamlss<-function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(pGA(q=y+tau,mu=pred_mu,sigma=pred_sigma))
}

mu_GA_Gamlss<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
    }))
    return(pred_mu - tau)
}

sim_GA_Gamlss<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    simy<-replicate(1, rGA(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
    simy[simy<0]<-0
    return(simy)
}

##Calculate Gamma density for Gamma GAMLSS
dens_LN_Gamlss<-function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(dLNO(x=y,mu=pred_mu,sigma=pred_sigma))
}

##Calculate Log-Normal CDF for Log-Normal GAMLSS
CDF_LN_Gamlss<-function(y, model, newdata){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    return(pLNO(q=y+tau,mu=pred_mu,sigma=pred_sigma))
}

mu_LN_Gamlss<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    mean<-exp(pred_mu+1/2*pred_sigma^2)
    return(mean - tau)
}

sim_LN_Gamlss<-function(model, newdata, tau){
    invisible(capture.output({
        pred_mu<-predict(model,what="mu",newdata=newdata,type="response")
        pred_sigma<-predict(model,what="sigma",newdata=newdata,type="response")
    }))
    simy<-replicate(1, rLNO(n=length(pred_mu),mu=pred_mu,sigma=pred_sigma)-tau)
    simy[simy<0]<-0
    return(simy)
}

################################################################################
### Special implementation for hierarchical models
################################################################################

vlookup_N_i <- function(N, data){

    N_rep <- c()
    occurence_years <- as.numeric(as.character(data$origin))

    for (i in occurence_years) {
        N_rep <- c(N_rep, N[i])
    }

    return(N_rep)
}

calc_PPCI <- function(N, data) {
    ##Fit a second model for paid loss
    ##Payment per Notified Claim
    ###Repeat each element in N_rep by the number of entries in each accident period
    N_rep<-vlookup_N_i(N,data)
    PPCI=data[order(as.numeric(data$origin)),]$aggregate_claims/N_rep

    return (PPCI)
}

dens_PPCI <- function(y, newdata, train_data, tri.size) {

    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    ###PPCI
    ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period
    PPCI <- calc_PPCI(N, train_data)
    fit_ODP_ppci<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=train_data)

    N_rep<-vlookup_N_i(N, newdata)
    pred_ODP<-predict.glm(fit_ODP_ppci,newdata=newdata,type="response",se.fit=TRUE)
    mu<-pred_ODP$fit*N_rep
    phi<-(pred_ODP$residual.scale)^2*(N_rep)^2
    dens<-c()
    for (i in 1:length(y)){
        dens[i]<-dODP(y=y[i],lambda=mu[i],phi=phi[i])
    }
    dens
}

CDF_PPCI<-function(y, newdata, train_data, tri.size){

    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    ###PPCI
    ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period
    PPCI <- calc_PPCI(N, train_data)
    fit_ODP_ppci<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=train_data)

    N_rep<-vlookup_N_i(N, newdata)
    pred_ODP<-predict.glm(fit_ODP_ppci,newdata=newdata,type="response",se.fit=TRUE)
    mu<-pred_ODP$fit*N_rep
    phi<-(pred_ODP$residual.scale)^2*(N_rep)^2

    return (pODP(y=y, lambda=mu,phi=phi))
}

mu_PPCI<-function(newdata, train_data, tri.size){
    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    ###PPCI
    ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period
    PPCI <- calc_PPCI(N, train_data)
    fit_ODP_ppci<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=train_data)

    N_rep<-vlookup_N_i(N, newdata)
    pred_ODP<-predict.glm(fit_ODP_ppci,newdata=newdata,type="response",se.fit=TRUE)
    mu<-pred_ODP$fit*N_rep
    return(mu)
}

sim_PPCI<-function(newdata, train_data, tri.size){

    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    ###PPCI
    ##Fit a ODP/Chain-Ladder model to predict total claim notification count at each accident period
    PPCI <- calc_PPCI(N, train_data)
    fit_ODP_ppci<-glm(PPCI~factor(dev),family=quasipoisson(link="log"),data=train_data)

    N_rep<-vlookup_N_i(N, newdata)
    pred_ODP<-predict.glm(fit_ODP_ppci,newdata=newdata,type="response",se.fit=TRUE)
    mu<-pred_ODP$fit*N_rep
    phi<-(pred_ODP$residual.scale)^2*(N_rep)^2
    simy<-replicate(1, rtweedie(length(mu),xi=1,mu=mu,phi=phi))
    return(simy)
}

calc_PPCF <- function(N, data) {
    PPCF=data$aggregate_claims/data$settle_count
    #Create a column for operation time
    ##Fit a second model for paid loss
    ##Payment per Notified Claim
    ###Repeat each element in N_rep by the number of entries in each accident period
    N_rep<-vlookup_N_i(N,data)
    OT=data$cum_settle_count/N_rep
    #Remove the cells with zero finalized claim count but with positive payments;
    PPCF_data <- data.frame(list(PPCF = PPCF, OT=OT))
    PPCF_data<-na.omit(PPCF_data[PPCF_data$PPCF!=Inf,])

    return (PPCF_data)
}

calc_pred_OT <- function(model_subCount, N, data, newdata, tri.size) {
    N_rep<-vlookup_N_i(N,newdata)
    # Predicted Finalised Claims for Each dataset
    # Assumes model only uses dev
    pred_F<-predict(model_subCount,newdata=data.frame(dev=1:tri.size),type="response")

    cum_F_count_train <-  matrix(0, nrow = tri.size, ncol = tri.size)
    for (i in 1:tri.size) {
        for (j in 1:tri.size) {
            # If cumulative settled count does not exist in train dataset
            # Add the most recent cumulative settled count from previous dev year and the predicted settle count for the year
            # Assumes that there is always a value for dev year 1
            cum_F_count_train[i, j] <- ifelse(nrow(data[(data$origin==i)&(data$dev==j), ]) == 0,
                                              cum_F_count_train[i, j-1] + pred_F[j],
                                              data[(data$origin==i)&(data$dev==j), 'cum_settle_count'])
        }
    }

    pred_F_cum <- unlist(mapply(function(i, j) cum_F_count_train[i, j],
                                as.numeric(as.character(newdata$origin)),
                                as.numeric(as.character(newdata$dev))
    ))

    pred_OT<-pred_F_cum/N_rep

    return(pred_OT)
}

dens_PPCF <- function(y, newdata, train_data, tri.size) {

    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    odp_FC<-glm(settle_count~factor(dev),data=train_data,family=quasipoisson(link="log"))
    PPCF_data <- calc_PPCF(N, train_data)
    ODP_PPCF<-glm(PPCF~OT,family=quasipoisson(link="log"),data=PPCF_data)

    pred_F<-predict(odp_FC,newdata=newdata,type="response")
    pred_OT <- calc_pred_OT(odp_FC, N, train_data, newdata, tri.size)

    pred_payment<-predict(ODP_PPCF,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
    mu<-pred_payment$fit*pred_F
    phi<-(pred_payment$residual.scale)^2*(pred_F)^2

    dens_PPCF<-c()
    for (i in 1:length(y)){
        dens_PPCF[i]<-dODP(y=y[i],lambda=mu[i],phi=phi[i])
    }
    dens_PPCF[dens_PPCF==0]=min(dens_PPCF[dens_PPCF!=0])
    dens_PPCF
}

CDF_PPCF<-function(y, newdata, train_data, tri.size){

    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    odp_FC<-glm(settle_count~factor(dev),data=train_data,family=quasipoisson(link="log"))
    PPCF_data <- calc_PPCF(N, train_data)
    ODP_PPCF<-glm(PPCF~OT,family=quasipoisson(link="log"),data=PPCF_data)

    pred_F<-predict(odp_FC,newdata=newdata,type="response")
    pred_OT <- calc_pred_OT(odp_FC, N, train_data, newdata, tri.size)

    pred_payment<-predict(ODP_PPCF,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
    mu<-pred_payment$fit*pred_F
    phi<-mean((pred_payment$residual.scale)^2*(pred_F)^2)

    CDF_PPCF <- pODP(y=y,lambda=mu,phi=phi)
    return (ifelse(CDF_PPCF == 0, min(CDF_PPCF[CDF_PPCF!=0]), CDF_PPCF))
}

mu_PPCF<- function(newdata, train_data, tri.size){

    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    odp_FC<-glm(settle_count~factor(dev),data=train_data,family=quasipoisson(link="log"))
    PPCF_data <- calc_PPCF(N, train_data)
    ODP_PPCF<-glm(PPCF~OT,family=quasipoisson(link="log"),data=PPCF_data)

    pred_F<-predict(odp_FC,newdata=newdata,type="response")
    pred_OT <- calc_pred_OT(odp_FC, N, train_data, newdata, tri.size)

    pred_payment<-predict(ODP_PPCF,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
    mu<-pred_payment$fit*pred_F
    return(mu)
}

sim_PPCF<-function(newdata, train_data, tri.size){

    ### Claims Count Information
    fit_nc<-glm(notif_count~factor(origin)+factor(dev),data=train_data,family=quasipoisson(link="log"))
    ##Obtain Estimated notification claims count for each accident period: Nk
    N<-c()
    N[1]<-sum(train_data[train_data$origin==1,]$notif_count)
    for (i in 2:tri.size){
        N[i]<- sum(train_data[train_data$origin==i,]$notif_count) +
            sum(round(predict(fit_nc,newdata=newdata[newdata$origin==i,],type="response"),0))
    }

    odp_FC<-glm(settle_count~factor(dev),data=train_data,family=quasipoisson(link="log"))
    PPCF_data <- calc_PPCF(N, train_data)
    ODP_PPCF<-glm(PPCF~OT,family=quasipoisson(link="log"),data=PPCF_data)

    pred_F<-predict(odp_FC,newdata=newdata,type="response")
    pred_OT <- calc_pred_OT(odp_FC, N, train_data, newdata, tri.size)

    pred_payment<-predict(ODP_PPCF,newdata=data.frame(OT=pred_OT),type="response",se.fit=TRUE)
    mu<-pred_payment$fit*pred_F
    phi<-(pred_payment$residual.scale)^2*(pred_F)^2
    simy<-replicate(1, rtweedie(length(mu), xi=1,mu=mu,phi=phi))
    return(simy)
}
