decompaction <- function(phi_0,c,y,y1) {

    # This function decompacts sedimentary formations in order to correct for porosity loss during compaction.
    # Inputs = (Surface porosity, Porosity-depth coefficent, Present-day depth/thickness profile [top of fm., base of fm.], Starting depth after removal of overlying layer)
    # Outputs = (volume of porosity, volume of sediments, new depth-to-basement, thickness of decompacted fm., average porosity of fm.)


    yw <- (phi_0/c)*(exp(-c*y[1]) - exp(-c*y[2]))
    ytot <- y[2] - y[1]
    ys <- ytot - yw

    track <- 0
    tolerance <- 0.01
    numIterations = seq(ytot,ytot*1000,by=0.01)
    while(track == 0) {
        for(i in numIterations) {

            LHS <- i-y1
            RHS <- ys + (phi_0/c)*(exp(-c*y1) - exp(-c*i))

            if (abs(LHS-RHS) <= tolerance) {
                y2 = i
                track <- 1
                break
            }
        }
        tolerance <- tolerance + 0.5
    }

    thickness <- y2-y1
    avg_porosity <- ((phi_0/c)*(exp(-c*y1)-exp(-c*y2))) / (y2-y1)


    data = c(yw,ys,y2,thickness,avg_porosity)
    return(data)

}


backstripping <- function(data) {

    # This function removes the sediment loading effect assuing Airy Isostasy
    # Inputs = a data table with columnns:
        # 1. Stratigraphic Units, 2. Phi_0, 3. Porosity depth coefficient, 4. rho_sg, 5. age at top of fm.,
        # 6. age at base of fm., 7. Present day depth (top), 8. Present day depth (bottom), 9. Water depth correction, 10. Eustatic sea level change
    # Outputs = (1. Decompacted Thickness, 2. Average Porosity, 3. Bulk Density, 4. Bulk Density of column, 5. Water-loaded tectonic subsidence)

    ## Convert data into parameter arrays
    phi_0 <- data[,2]                                                                   # Surface porosity
    c <- data[,3]                                                                       # Porosity Depth Coefficient                                                [km^-1]
    c <- c/1000                                                                         # Porosity Depth coefficients                                               [m^-1]
    rho_sg <- data[,4]                                                                  # Sediment grain Density                                                    [kg/m^3]
    rho_w <- 1030                                                                       # Density of seawater                                                       [kg/m^3]
    rho_m <- 3300                                                                       # Density of the mantle                                                     [kg/m^3]
    upper_age <- data[,5]                                                               # Age at the top of formation                                               [Ma]
    base_age <- data[nrow(data),6]                                                      # Age at the base of stratigraphic section                                  [Ma]
    ages <- c(data[nrow(data),6],rev(data[,5]))                                         # Age array from base to top of column                                      [Ma]
    strat_units <- data[,1]                                                             # Names of the stratigraphic formations
    Wd <- data[,9]                                                                      # Water depth corrections                                                   [m]
    dSL <- data[,10]                                                                    # Eustatic sea level change                                                 [m]
    y <- matrix(0,nrow(data),2)
    for(i in seq(1,nrow(data))) {
        y[i,] <- c(data[i,7],data[i,8])                                                 # Present-day depths (top,base)                                             [m]
    }


    ## Decompact the sedimentary formations

    # Prepare matrices for results
    rows_dt <- c(strat_units,"Stratigraphic Thickness")
    columns_dt <- rev(upper_age)
    decompacted_thickness <- matrix(0,length(rows_dt),length(columns_dt))               # Thickness of each formation after each stage of decompaction              [m]
    horizon_thicknesses <- matrix(0,nrow(data),length(ages))                            # Stratigraphic thicknesses of each age horizon                             [m]
    avg_porosity <- matrix(0,length(strat_units),length(columns_dt))                    # Average porosities of each formation after each stage of decompaction     [m]
    rho_b <- matrix(0,length(strat_units),length(columns_dt))                           # Bulk density of each formation after each stage of decompaction           [kg/m^3]

    # Decompact the formations
    counter <- length(c)
    for(i in seq(1,length(columns_dt))) {
        y1 <- 0
        for(j in seq(1,length(c))) {
            if(j >= counter) {
                results <- decompaction(phi_0[j],c[j],y[j,],y1)
                decompacted_thickness[j,i] <- results[4]
                avg_porosity[j,i] <- results[5]
                rho_b[j,i] <- results[5]*rho_w + (1-results[5])*rho_sg[j]
                y1 <- results[3]
            }
            if(j == length(c)) {
                decompacted_thickness[length(rows_dt),i] <- sum(decompacted_thickness[1:length(c),i])
                counter <- counter-1
            }
        }
    }

    # Determine the decompacted depths at each age horizon
    for(i in seq(nrow(horizon_thicknesses),1)) {
        for(j in seq(1,ncol(decompacted_thickness))) {
            horizon_thicknesses[i,j+1] <- sum(decompacted_thickness[1:i,j])
        }
    }


    ## Calculate the water-loaded tectonic subsidence

    # Calculate the bulk density of the sedimantary column at each stage of deposition
    rho_b_column <- matrix(0,length(rows_dt),length(columns_dt))                        # Bulk Density of the sedimentary column                                    [kg/m^3]

    for(i in seq(1,ncol(decompacted_thickness))) {
        for(j in seq(1,nrow(rho_b))) {
            rho_b_column[j,i] <- rho_b[j,i]*decompacted_thickness[j,i]/decompacted_thickness[length(rows_dt),i]
        }
        rho_b_column[length(rows_dt),i] <- sum(rho_b_column[1:length(c),i])
    }


    # Calculate water-loaded subsidence
    yw <- matrix(0,length(strat_units),1)                                               # Water-loaded tectonic subsidence                                          [m]
    Wd <- rev(Wd)
    dSL <- rev(dSL)

    for(i in seq(1,length(strat_units))) {
        yw[i] <- decompacted_thickness[length(rows_dt),i]*((rho_m-rho_b_column[length(rows_dt),i])/(rho_m-rho_w)) - dSL[i]*(rho_m/(rho_m-rho_w)) + (Wd[i] - dSL[i])
    }


    dt_df <- data.frame(decompacted_thickness)
    rownames(dt_df) <- rows_dt
    colnames(dt_df) <- columns_dt

    ht_df <- data.frame(horizon_thicknesses)
    rownames(ht_df) <- strat_units
    colnames(ht_df) <- ages

    ap_df <- data.frame(avg_porosity)
    rownames(ap_df) <- strat_units
    colnames(ap_df) <- columns_dt

    rb_df <- data.frame(rho_b)
    rownames(rb_df) <- strat_units
    colnames(rb_df) <- columns_dt

    rbc_df <- data.frame(rho_b_column)
    rownames(rbc_df) <- c(strat_units,"Bulk Density of Column")
    colnames(rbc_df) <- columns_dt

    out <- list(dt_df,ht_df,ap_df,rb_df,rbc_df,yw)


    return(out)

}





thermal_subsidence <- function(B,duration) {


    ## This function calculates the post-rift thermal subsidence

    ## Inputs = (Stretching factor, # of years after rifting)



    ## Initial parameters:
    y_l <- 125000                                            # Initial lithospheric thickness in m          [m]
    rho_m0 <- 3330                                           # Density of the mantle at 0 degrees celcius   [kg/m^3]
    rho_w <- 1030                                            # Density of water                             [kg/m^3]
    alpha_v <- 3.28e-5                                       # volumetric coefficient of thermal expansion  [1/K]
    Tm <- 1333                                               # Temperature of the mantle                    [C]
    kappa <- 1e-6                                            # Thermal diffusivity                          [m^2/s]


    ## Calculate post-rift thermal subsidence
    t <- seq(0,duration,by=1)                                # Time since end of rifting                    [Myrs]
    time_s <- t*1e6*365*24*60*60                             # Time since end of rifting                    [s]

    E0 <- 4*y_l*rho_m0*alpha_v*Tm/(pi^2*(rho_m0-rho_w))
    tau <- y_l^2/(pi^2*kappa)                                # Thermal time constant                        [s]

    S_thermal <- E0*(B/pi)*sin(pi/B)*(1-exp(-time_s/tau))    # Thermal subsidence                           [m]


    ## Return results
    results <- rbind(t,S_thermal)
    return(results)



} # End Function



thermal_subsidence1 <- function(B,duration) {


    ## This function calculates the post-rift thermal subsidence

    ## Inputs = (Stretching factor, # of years after rifting)



    ## Initial parameters:
    y_l <- 125000                                            # Initial lithospheric thickness in m          [m]
    rho_m0 <- 3330                                           # Density of the mantle at 0 degrees celcius   [kg/m^3]
    rho_w <- 1030                                            # Density of water                             [kg/m^3]
    alpha_v <- 3.28e-5                                       # volumetric coefficient of thermal expansion  [1/K]
    Tm <- 1333                                               # Temperature of the mantle                    [C]
    kappa <- 1e-6                                            # Thermal diffusivity                          [m^2/s]


    ## Calculate post-rift thermal subsidence
    t <- seq(0,duration,by=0.001)                             # Time since end of rifting                    [Myrs]
    time_s <- t*1e6*365*24*60*60                             # Time since end of rifting                    [s]

    E0 <- 4*y_l*rho_m0*alpha_v*Tm/(pi^2*(rho_m0-rho_w))
    tau <- y_l^2/(pi^2*kappa)                                # Thermal time constant                        [s]

    S_thermal <- E0*(B/pi)*sin(pi/B)*(1-exp(-time_s/tau))    # Thermal subsidence                           [m]


    ## Return results
    results <- rbind(t,S_thermal)
    return(results)



} # End Function




beta_factor <- function(ages,yw) {

    ## This function performs linear regression on stratigraphic ages with their associated water-loaded tectonic subsidence
    ## to obtain a stretching factor for the stratigraphic section of interest.

    ## Inputs = (Ages [Ma], Water-loaded tectonic subsidence [m])




    ## Initial parameters:
    y_l <- 125000                               # Initial lithospheric thickness in m          [m]
    rho_m0 <- 3330                              # Density of the mantle at 0 degrees celcius   [kg/m^3]
    rho_w <- 1030                               # Density of water                             [kg/m^3]
    alpha_v <- 3.28e-5                          # volumetric coefficient of thermal expansion  [1/K]
    Tm <- 1333                                  # Temperature of the mantle                    [C]
    kappa <- 1e-6                               # Thermal diffusivity                          [m^2/s]



    ## Convert ages to time
    t <- matrix(0,length(ages),1)               # Time since end of rifting                    [Myrs]
    t[1] <- 0
    for(i in seq(2,length(ages))) {
        t_increase <- ages[i-1] - ages[i]
        t[i] <- t[i-1] + t_increase
    }



    ## Calculate lithospheric time constant (tau)
    tau <- y_l^2/(pi^2*kappa)                   # Lithospheric time constant                   [s]
    tau <- tau/(60*60*24*365*1e6)               # Lithospheric time constant                   [Myrs]



    ## Plot thermal subsidence as a function of 1-exp(-t/tau) and perform linear regression
    x <- 1-exp(-t/tau)
    b <- matrix(1,length(x),1)
    A <- cbind(b,x)
    y <- c(0,yw)
    l <- solve(t(A)%*%A) %*% t(A) %*% y

    plot(x,y,col="blue",type="p",pch=16,
         main="Linear Regression of Thermal Subsidence Data",
         xlab="1-exp(-t/tau)", ylab="Water-Loaded Thermal Subsidence [m]")
    xreg <- seq(0,1,by=0.01)
    yreg <- matrix(0,length(xreg),1)
    for(i in seq(1,length(xreg))) {
        yreg[i] <- l[2]*xreg[i] + l[1]
    }
    lines(xreg,yreg,type="l",col="black")
    legend("bottomright",
           legend = c("Data","Linear Fit"),
           col=c("blue","black"), lty=c(NA,1),pch=c(16,NA))


    ## Calculate stretch factor (B)
    E0 <- (4*y_l*rho_m0*alpha_v*Tm)/(pi^2*(rho_m0-rho_w))
    B <- 0
    tolerance <- 0.01
    while(B == 0) {
        for(i in seq(1,3,by=0.001)) {
            LHS <- (E0*i/pi)*sin(pi/i)
            if(abs(LHS-l[2]) <= tolerance) {
                B <- i
                break
            }
        }
        tolerance <- tolerance + 0.01
    }


    #out <- list(summary(linear_fit),B)
    out <- list(l,B)

    return(out)




} # End Function
