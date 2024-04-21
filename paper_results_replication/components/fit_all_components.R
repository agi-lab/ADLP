### Perform the simulations:

### Create a list to store the simulated component models:
components_list <- list()

### Create a matrix to store selected individual models' Log Score:
LS_ZAGA_mat <- matrix(NA, nrow = 780, ncol = n.sims)
LS_PPCF_mat <- matrix(NA, nrow = 780, ncol = n.sims)

suppressWarnings({
for (sim in 1:n.sims) {
    # Train on train, test on valid-test
    set.seed(20200130+sim)
    past_data <- read.csv(sprintf('simulation/triangle_%s-data/sim%s-past-data.csv', tri.size, sim))
    full_data <- read.csv(sprintf('simulation/triangle_%s-data/sim%s-full-data.csv', tri.size, sim))
    past_data$tau_Ga<-tau_Ga
    past_data$tau_LN<-tau_LN
    full_data$tau_Ga<-tau_Ga
    full_data$tau_LN<-tau_LN

    #### Partition the data into training and validation set #####
    #insample_data <- full_data[full_data$calendar <= 41,]
    insample_data <- past_data
    outsample_data <- full_data[full_data$calendar > 41,]
    train_val <- train_val_split(insample_data)

    ############################################################
    ## Create Fit Components ###################################
    ############################################################

    train_data <- train_val$train
    train_all_data <- rbind(train_val$train, train_val$valid)


    #### Fit all the component models used in the paper: ######
    invisible(capture.output({
        ODP_GLM<-glm(formula=aggregate_claims~factor(origin)+factor(dev),family=quasipoisson(link="log"),data=train_data)
        Ga_optimTau<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(origin)+factor(dev),data=train_data,family=GA(mu.link="log", sigma.link ="log"))
        LN_optimTau<-gamlss(formula=(aggregate_claims+tau_LN)~factor(origin)+factor(dev),data=train_data,family=LOGNO(mu.link="identity",sigma.link="log"))
        gamma_1<-gamlss(formula=aggregate_claims~factor(origin)+factor(dev),nu.formula=~as.numeric(as.character(dev)),data=train_data,family=ZAGA(mu.link="log",sigma.link = "log", nu.link = "logit"))
        LN_1<-gamlss.inf::gamlssZadj(y=aggregate_claims,mu.formula = ~factor(origin)+factor(dev),xi0.formula=~as.numeric(as.character(dev)),data=train_data,family=LOGNO(mu.link="identity",sigma.link="log"))
        glm_ODP_Ho_tr1<-glm(formula=aggregate_claims~factor(origin)+log(dev)+dev,family=quasipoisson(link="log"),data=train_data)
        glm_Ga_Ho_tr1<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(origin)+log(dev)+dev,data=train_data,family=GA(mu.link="log", sigma.link ="log"))
        glm_LN_Ho_tr1<-gamlss(formula=(aggregate_claims+tau_LN)~factor(origin)+log(dev)+dev,data=train_data,family=LOGNO(mu.link="identity",sigma.link="log"))
        glm_ODP_Cal_tr1<-glm(formula=aggregate_claims~factor(dev)+calendar,family=quasipoisson(link="log"),data=train_data)
        glm_Ga_Cal_tr1<-gamlss(formula=(aggregate_claims+tau_Ga)~factor(dev)+calendar,data=train_data,family=GA(mu.link="log", sigma.link ="log"))
        glm_LN_Cal_tr1<-gamlss(formula=(aggregate_claims+tau_LN)~factor(dev)+calendar,data=train_data,family=LOGNO(mu.link="identity",sigma.link="log"))
        sp_Normal<-gamlss(formula=aggregate_claims~scs(origin)+scs(dev),data=train_data,family=NO(),trace=FALSE)
        sp_Gamma<-gamlss(formula=(aggregate_claims+tau_Ga)~scs(origin)+scs(dev),data=train_data,family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
        sp_LN<-gamlss(formula=(aggregate_claims+tau_LN)~scs(origin)+scs(dev),data=train_data,family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)
        gamlss_GA<-gamlss(formula=(aggregate_claims+tau_Ga)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=train_data,sigma.formula=~cs(as.numeric(as.character(dev))),family=GA(mu.link="log", sigma.link ="log"),trace=FALSE)
        gamlss_LN<-gamlss(formula=(aggregate_claims+tau_LN)~scs(as.numeric(as.character(origin)))+scs(as.numeric(as.character(dev))),data=train_data,sigma.formula=~cs(as.numeric(as.character(dev))),family=LOGNO(mu.link="identity",sigma.link="log"),trace=FALSE)

        ODP_GLM_full <- update(ODP_GLM, data = train_all_data)
        Ga_optimTau_full <- update(Ga_optimTau, data = train_all_data)
        LN_optimTau_full <- update(LN_optimTau, data = train_all_data)
        gamma_1_full <- update(gamma_1, data = train_all_data)
        LN_1_full <- update(LN_1, data = train_all_data)
        glm_ODP_Ho_tr1_full <- update(glm_ODP_Ho_tr1, data = train_all_data)
        glm_Ga_Ho_tr1_full <- update(glm_Ga_Ho_tr1, data = train_all_data)
        glm_LN_Ho_tr1_full <- update(glm_LN_Ho_tr1, data = train_all_data)
        glm_ODP_Cal_tr1_full <- update(glm_ODP_Cal_tr1, data = train_all_data)
        glm_Ga_Cal_tr1_full <- update(glm_Ga_Cal_tr1, data = train_all_data)
        glm_LN_Cal_tr1_full <- update(glm_LN_Cal_tr1, data = train_all_data)
        sp_Normal_full <- update(sp_Normal, data = train_all_data)
        sp_Gamma_full <- update(sp_Gamma, data = train_all_data)
        sp_LN_full <- update(sp_LN, data = train_all_data)
        gamlss_GA_full <- update(gamlss_GA, data = train_all_data)
        gamlss_LN_full <- update(gamlss_LN, data = train_all_data)
    }))


    # Custom model class used for PPCI and PPCF
    ODP_PPCI<-custom_model(aggregate_claims~., data=train_data)
    ODP_PPCF<-custom_model(aggregate_claims~., data=train_data)

    ODP_PPCI_full <- update(ODP_PPCI, data = past_data)
    ODP_PPCF_full <- update(ODP_PPCF, data = past_data)


    ###########################################################
    # Create ADLP Components List #############################
    ###########################################################

    ODP_GLM_component = adlp_component(
        model_train = ODP_GLM,
        model_full = ODP_GLM_full,
        calc_dens = dens_ODP_GLM,
        calc_cdf = CDF_ODP_GLM,
        calc_mu = mu_ODP_GLM,
        sim_fun = sim_ODP_GLM
    )

    Ga_optimTau_component = adlp_component(
        model_train = Ga_optimTau,
        model_full = Ga_optimTau_full,
        calc_dens = dens_GA,
        calc_cdf = CDF_GA,
        calc_mu = mu_GA,
        sim_fun = sim_GA,
        tau = tau
    )

    LN_optimTau_component = adlp_component(
        model_train = LN_optimTau,
        model_full = LN_optimTau_full,
        calc_dens = dens_LN,
        calc_cdf = CDF_LN,
        calc_mu = mu_LN,
        sim_fun = sim_LN,
        tau = tau
    )

    gamma_1_component = adlp_component(
        model_train = gamma_1,
        model_full = gamma_1_full,
        calc_dens = dens_ZAGA,
        calc_cdf = CDF_ZAGA,
        calc_mu = mu_ZAGA,
        sim_fun =sim_ZAGA
    )

    LN_1_component = adlp_component(
        model_train = LN_1,
        model_full = LN_1_full,
        calc_dens = dens_ZALN,
        calc_cdf = CDF_ZALN,
        calc_mu = mu_ZALN,
        sim_fun = sim_ZALN
    )

    glm_ODP_Ho_tr1_component = adlp_component(
        model_train = glm_ODP_Ho_tr1,
        model_full = glm_ODP_Ho_tr1_full,
        calc_dens = dens_ODP_GLM,
        calc_cdf = CDF_ODP_GLM,
        calc_mu = mu_ODP_GLM,
        sim_fun = sim_ODP_GLM
    )

    glm_Ga_Ho_tr1_component = adlp_component(
        model_train = glm_Ga_Ho_tr1,
        model_full = glm_Ga_Ho_tr1_full,
        calc_dens = dens_GA,
        calc_cdf = CDF_GA,
        calc_mu = mu_GA,
        sim_fun = sim_GA,
        tau = tau
    )

    glm_LN_Ho_tr1_component = adlp_component(
        model_train = glm_LN_Ho_tr1,
        model_full = glm_LN_Ho_tr1_full,
        calc_dens = dens_LN,
        calc_cdf = CDF_LN,
        calc_mu = mu_LN,
        sim_fun = sim_LN,
        tau = tau
    )

    glm_ODP_Cal_tr1_component = adlp_component(
        model_train = glm_ODP_Cal_tr1,
        model_full = glm_ODP_Cal_tr1_full,
        calc_dens = dens_ODP_GLM,
        calc_cdf = CDF_ODP_GLM,
        calc_mu = mu_ODP_GLM,
        sim_fun = sim_ODP_GLM
    )

    glm_Ga_Cal_tr1_component = adlp_component(
        model_train = glm_Ga_Cal_tr1,
        model_full = glm_Ga_Cal_tr1_full,
        calc_dens = dens_GA,
        calc_cdf = CDF_GA,
        calc_mu = mu_GA,
        sim_fun = sim_GA,
        tau = tau
    )

    glm_LN_Cal_tr1_component = adlp_component(
        model_train = glm_LN_Cal_tr1,
        model_full = glm_LN_Cal_tr1_full,
        calc_dens = dens_LN,
        calc_cdf = CDF_LN,
        calc_mu = mu_LN,
        sim_fun = sim_LN,
        tau = tau
    )

    sp_Normal_component = adlp_component(
        model_train = sp_Normal,
        model_full = sp_Normal_full,
        calc_dens = dens_NO,
        calc_cdf = CDF_Normal,
        calc_mu = mu_Normal,
        sim_fun = sim_NO
    )

    sp_Gamma_component = adlp_component(
        model_train = sp_Gamma,
        model_full = sp_Gamma_full,
        calc_dens = dens_GA,
        calc_cdf = CDF_GA,
        calc_mu = mu_GA,
        sim_fun = sim_GA,
        tau = tau
    )

    sp_LN_component = adlp_component(
        model_train = sp_LN,
        model_full = sp_LN_full,
        calc_dens = dens_LN,
        calc_cdf = CDF_LN,
        calc_mu = mu_LN,
        sim_fun = sim_LN,
        tau = tau
    )

    gamlss_GA_component = adlp_component(
        model_train = gamlss_GA,
        model_full = gamlss_GA_full,
        calc_dens = dens_GA_Gamlss,
        calc_cdf = CDF_GA_Gamlss,
        calc_mu = mu_GA_Gamlss,
        sim_fun = sim_GA_Gamlss,
        tau = tau
    )

    gamlss_LN_component = adlp_component(
        model_train = gamlss_LN,
        model_full = gamlss_LN_full,
        calc_dens = dens_LN_Gamlss,
        calc_cdf = CDF_LN_Gamlss,
        calc_mu = mu_LN_Gamlss,
        sim_fun = sim_LN_Gamlss,
        tau = tau
    )


    ODP_PPCI_component = adlp_component(
        model_train = ODP_PPCI,
        model_full = ODP_PPCI_full,
        calc_dens = dens_PPCI,
        calc_cdf = CDF_PPCI,
        calc_mu = mu_PPCI,
        sim_fun = sim_PPCI,
        train_data = train_data,
        tri.size = 40
    )

    ODP_PPCF_component = adlp_component(
        model_train = ODP_PPCF,
        model_full = ODP_PPCF_full,
        calc_dens = dens_PPCF,
        calc_cdf = CDF_PPCF,
        calc_mu = mu_PPCF,
        sim_fun = sim_PPCF,
        train_data = train_data,
        tri.size = 40
    )

    LS_ZAGA_mat[, sim] <- log(dens_ZAGA(outsample_data$aggregate_claims,gamma_1_full,outsample_data)+1e-6)
    LS_PPCF_mat[, sim] <- log(dens_PPCF(y = outsample_data$aggregate_claims, newdata = outsample_data, train_data = train_all_data, tri.size = 40)+1e-6)

    components_list[[sim]] <- adlp_components(
        ODP_GLM = ODP_GLM_component,
        Ga_optimTau = Ga_optimTau_component,
        LN_optimTau = LN_optimTau_component,
        gamma_1 = gamma_1_component,
        LN_1 = LN_1_component,
        glm_ODP_Ho_tr1 = glm_ODP_Ho_tr1_component,
        glm_Ga_Ho_tr1 = glm_Ga_Ho_tr1_component,
        glm_LN_Ho_tr1 = glm_LN_Ho_tr1_component,
        glm_ODP_Cal_tr1 = glm_ODP_Cal_tr1_component,
        glm_Ga_Cal_tr1 = glm_Ga_Cal_tr1_component,
        glm_LN_Cal_tr1 = glm_LN_Cal_tr1_component,
        sp_Normal = sp_Normal_component,
        sp_Gamma = sp_Gamma_component,
        sp_LN = sp_LN_component,
        gamlss_GA = gamlss_GA_component,
        gamlss_LN = gamlss_LN_component,
        ODP_PPCI = ODP_PPCI_component,
        ODP_PPCF = ODP_PPCF_component
    )




}
  })

