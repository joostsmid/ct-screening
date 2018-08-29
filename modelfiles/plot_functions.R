### Here we define functions designed to analyze the results of MainModel_STI.R
###
########################################################################################

### auxiliary function to choose nice colours from Color Brewer palette 
colours <- function(number){
  
  if(number == 1){
    return(c("blue"))
  }
  if(number == 2){
    return(c("blue","red"))
  }
  if(number > 2 && number < 10){
    return(brewer.pal(number,"Set1"))
  }
  if(number > 10){
    #return(sample(colours()[600:657],strain.number, replace = TRUE))
    #return(sample(rainbow(strain.number), replace = FALSE))
    #return(sample(terrain.colors(strain.number), replace = FALSE))
    return(sample(rainbow(number), replace = FALSE))
  }
}

### auxiliary function to choose nice colours from Color Brewer palette 
sim.cols <- function(number){
  
  if(number == 1){
    return(c("blue"))
  }
  if(number == 2){
    return(c("blue","red"))
  }
  if(number > 2 && number < 6){
    return(brewer.pal(number,"Set2"))
  }
  if(number > 5 && number < 16){
    half.numb <- floor(number/2)
    ### females 
    fems <- brewer.pal(half.numb,"Set2")
    ### males 
    mals <- brewer.pal(half.numb,"Pastel2")
    
    out <- as.vector(rbind(fems, mals))
    
    return(out)
  }
  if(number > 10){
    #return(sample(colours()[600:657],strain.number, replace = TRUE))
    #return(sample(rainbow(strain.number), replace = FALSE))
    #return(sample(terrain.colors(strain.number), replace = FALSE))
    return(sample(rainbow(number), replace = FALSE))
  }
}

# Multiple plot function
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

### plot the density distribution of individuals in a 3-dimensional tensor across age and activity classes, 
### separated by different genders
### normalization takes places over all indices sex, activity class and age!
plot.3tensor <- function(tensor.object, cex = 2, aggregatebySex = FALSE, aggregatebyAC = FALSE, normalize = TRUE, doublewindow = TRUE,plottype = "o"){
  
  ### extract the information about the tensor object
  genders <- gsub("(.*)", '\\1' ,dimnames(tensor.object)$sex)
  activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(tensor.object)$j))
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(tensor.object)$age))
  age.range <- range(ages)
  x.lim = age.range
  
  ##############################################################################
  # if aggregated by both sex and activity class (then 1 plot)
  if(aggregatebySex & aggregatebyAC){
    tensor.object <- margin.tensor(tensor.object,by=3)
    
    ### normalization (if required)
    if(normalize){
      tot <- sum(tensor.object)
      tensor.object <- tensor.object/tot
    }
    
    ### open plotting device
    par(mfrow = c(1,1))
    ymax <- max(tensor.object) + .2 * max(tensor.object)
    y.lim <- c(0,ymax)
    plot(min(age.range), 0, ylim = y.lim , xlim = x.lim, 
         type = "n", main = "all", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
    lines(ages, tensor.object, type = plottype, col = "red", lwd = cex)
    
    ### add axes
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
  }
  
  ##############################################################################
  # if aggregated by sex only (then 1 plot)
  else if(aggregatebySex &! aggregatebyAC){
    tensor.object <- margin.tensor(tensor.object,i=1)
    
    ### normalization (if required)
    if(normalize){
      tot <- sum(tensor.object)
      tensor.object <- tensor.object/tot
    }
    
    ### open plotting device
    par(mfrow = c(1,1))
    ymax <- max(tensor.object) + .2 * max(tensor.object)
    y.lim <- c(0,ymax)
    plot(min(age.range), 0, ylim = y.lim , xlim = x.lim, 
         type = "n", main = "all", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
    
    ### loop through activity classes to plot distribution across age per activity class
    cols = colours(length(activity.classes))
    for(j in activity.classes){
      lines(ages, tensor.object[j,], type = plottype, col = cols[j], lwd = cex)
    }
    
    ### add axes
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
    
  }
  
  ##############################################################################
  # if aggregated by activity class only (then 2 plots)
  else if(aggregatebyAC &! aggregatebySex){
    tensor.object <- margin.tensor(tensor.object,i=2)
    
    ### normalization (if required)
    if(normalize){
      tot <- sum(tensor.object)
      tensor.object <- tensor.object/tot
    }
    
    ### open plotting device
    if (doublewindow){
      par(mfrow = c(1,2), oma = c(2, 2, 1, 1))
    } else {par(mfrow = c(1,1))}
    ymax <- max(tensor.object) + .2 * max(tensor.object)
    y.lim <- c(0,ymax)
    plot(min(age.range), 0, ylim = y.lim , xlim = x.lim, 
         type = "n", main = "females", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
    lines(ages, tensor.object["F",], type = plottype, col = "red", lwd = cex)
    
    ### add axes
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
    
    ### proceed by plotting males 
    
    plot(min(age.range), 0, ylim = y.lim , xlim = x.lim, 
         type = "n", main = "males", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex, lwd = cex)
    lines(ages, tensor.object["M",], type = plottype, col = "red", lwd = cex)
    
    
    ### axes
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
    
    ### axis labels
    mtext("Frequency", side=2, outer=T, line= 0, cex = cex)
    mtext("Age [years]", side=1, outer=T, line= 0, cex = cex)
  }
  
  ##############################################################################
  # if no aggregations (then 2 plots)
  else {
    
    if(normalize){
      ### normalize the entries of the tensor object
      tot <- sum(tensor.object)
      tensor.object <- tensor.object/tot
    }
    
    ### open plotting device
    if (doublewindow){
      par(mfrow = c(1,2), oma = c(2, 2, 1, 1))
    } else {par(mfrow = c(1,1))}
    ymax <- max(tensor.object) + .2 * max(tensor.object)
    y.lim <- c(0,ymax)
    plot(min(age.range), 0, ylim = y.lim , xlim = x.lim, 
         type = "n", main = "females", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
    
    ### colors for each activity class 
    cols = colours(length(activity.classes))
    
    ### loop through activity classes to plot distribution across age per activity class
    for(j in activity.classes){
      
      lines(ages, tensor.object["F",j,], type = plottype, col = cols[j], lwd = cex)
      
    }
    
    ### add axes
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
    
    ### proceed by plotting males 
    
    plot(min(age.range), 0, ylim = y.lim , xlim = x.lim, 
         type = "n", main = "males", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex, lwd = cex)
    
    
    ### loop through activity classes to plot distribution across age per activity class
    for(j in activity.classes){
      
      lines(ages, tensor.object["M",j,], type = plottype, col = cols[j], lwd = cex)
      
    }
    
    ### axes
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
    
    ### axis labels
    mtext("Frequency", side=2, outer=T, line= 0, cex = cex)
    mtext("Age [years]", side=1, outer=T, line= 0, cex = cex)
  }
  
}

### test 
#plot.3tensor(I, aggregatebySex = T, aggregatebyAC = T, normalize = T)
#plot.3tensor(I, aggregatebySex = F, aggregatebyAC = F, normalize = T)
#plot.3tensor(ct, aggregatebySex = F, aggregatebyAC = F, normalize = F)
#plot.3tensor(ct, aggregatebySex = F, aggregatebyAC = F, normalize = F, doublewindow = F)

### This function aggrgates the sim.result matrix by age and/or activity class
aggregate_simresult <- function( sim.result, comb.ages = FALSE, comb.AC = FALSE, n_age = length(age.classes)-1, normalize = FALSE){
  
  if(comb.ages & comb.AC){
    sim.result.new <- as.data.frame(sim.result[,1])
    
    sim.result.new$U_M <- rowSums(sim.result[,grep("U_M",colnames(sim.result))] )
    sim.result.new$U_F <- rowSums(sim.result[,grep("U_F",colnames(sim.result))] )
    sim.result.new$S_M <- rowSums(sim.result[,grep("S_M",colnames(sim.result))] )
    sim.result.new$S_F <- rowSums(sim.result[,grep("S_F",colnames(sim.result))] )
    sim.result.new$I_A_M <- rowSums(sim.result[,grep("I_A_M",colnames(sim.result))] )
    sim.result.new$I_A_F <- rowSums(sim.result[,grep("I_A_F",colnames(sim.result))] )
    sim.result.new$I_S_M <- rowSums(sim.result[,grep("I_S_M",colnames(sim.result))] )
    sim.result.new$I_S_F <- rowSums(sim.result[,grep("I_S_F",colnames(sim.result))] )
    sim.result.new$R_M <- rowSums(sim.result[,grep("R_M",colnames(sim.result))] )
    sim.result.new$R_F <- rowSums(sim.result[,grep("R_F",colnames(sim.result))] )
    sim.result.new$D_M <- rowSums(sim.result[,grep("D_M",colnames(sim.result))] )
    sim.result.new$D_F <- rowSums(sim.result[,grep("D_F",colnames(sim.result))] )
    sim.result <- sim.result.new
    
    if(normalize){
      sim.result[,-1] <- sim.result[,-1]/rowSums(sim.result[,-1])
    }
    
  } else if(comb.ages &! comb.AC){
    
    sim.result.new <- as.data.frame(sim.result[,1])
    sim.result.new$U_M <-    rowSums(sim.result[,grep("U_M",colnames(sim.result))] )
    sim.result.new$U_F <-    rowSums(sim.result[,grep("U_F",colnames(sim.result))] )
    sim.result.new$S_M_j1 <- rowSums(sim.result[,grep("S_M_j1",colnames(sim.result))] )
    sim.result.new$S_M_j2 <- rowSums(sim.result[,grep("S_M_j2",colnames(sim.result))] )
    sim.result.new$S_F_j1 <- rowSums(sim.result[,grep("S_F_j1",colnames(sim.result))] )
    sim.result.new$S_F_j2 <- rowSums(sim.result[,grep("S_F_j2",colnames(sim.result))] )
    
    sim.result.new$I_A_M_j1 <- rowSums(sim.result[,grep("I_A_M_j1",colnames(sim.result))] )
    sim.result.new$I_A_M_j2 <- rowSums(sim.result[,grep("I_A_M_j2",colnames(sim.result))] )
    sim.result.new$I_A_F_j1 <- rowSums(sim.result[,grep("I_A_F_j1",colnames(sim.result))] )
    sim.result.new$I_A_F_j2 <- rowSums(sim.result[,grep("I_A_F_j2",colnames(sim.result))] )
    
    sim.result.new$I_S_M_j1 <- rowSums(sim.result[,grep("I_S_M_j1",colnames(sim.result))] )
    sim.result.new$I_S_M_j2 <- rowSums(sim.result[,grep("I_S_M_j2",colnames(sim.result))] )
    sim.result.new$I_S_F_j1 <- rowSums(sim.result[,grep("I_S_F_j1",colnames(sim.result))] )
    sim.result.new$I_S_F_j2 <- rowSums(sim.result[,grep("I_S_F_j2",colnames(sim.result))] )
    
    sim.result.new$R_M_j1 <- rowSums(sim.result[,grep("R_M_j1",colnames(sim.result))] )
    sim.result.new$R_M_j2 <- rowSums(sim.result[,grep("R_M_j2",colnames(sim.result))] )
    sim.result.new$R_F_j1 <- rowSums(sim.result[,grep("R_F_j1",colnames(sim.result))] )
    sim.result.new$R_F_j2 <- rowSums(sim.result[,grep("R_F_j2",colnames(sim.result))] )
    
    sim.result.new$D_M_j1 <- rowSums(sim.result[,grep("D_M_j1",colnames(sim.result))] )
    sim.result.new$D_M_j2 <- rowSums(sim.result[,grep("D_M_j2",colnames(sim.result))] )
    sim.result.new$D_F_j1 <- rowSums(sim.result[,grep("D_F_j1",colnames(sim.result))] )
    sim.result.new$D_F_j2 <- rowSums(sim.result[,grep("D_F_j2",colnames(sim.result))] )
    sim.result <- sim.result.new
    
    if(normalize){
      sim.result[,-1] <- sim.result[,-1]/rowSums(sim.result[,-1])
    }
    
  } else if(!comb.ages & comb.AC){
    sim.result.new <- as.data.frame(sim.result[,1])
    age<-paste0('age', age.class.groups)
    for (k in 1:length(age)){
      sim.result.new[[paste0("U_M_",age[k])]] <- sim.result[,grep(paste0("U_M_",age[k]),colnames(sim.result))]
      sim.result.new[[paste0("S_M_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*S_M)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      sim.result.new[[paste0("I_A_M_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*I_A_M)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      sim.result.new[[paste0("I_S_M_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*I_S_M)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      sim.result.new[[paste0("R_M_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*R_M)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      sim.result.new[[paste0("D_M_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*D_M)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      
      sim.result.new[[paste0("U_F_",age[k])]] <- sim.result[,grep(paste0("U_F_",age[k]),colnames(sim.result))]
      sim.result.new[[paste0("S_F_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*S_F)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)]) #http://stackoverflow.com/questions/869809/combine-regexp
      sim.result.new[[paste0("I_A_F_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*I_A_F)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      sim.result.new[[paste0("I_S_F_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*I_S_F)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      sim.result.new[[paste0("R_F_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*R_F)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
      sim.result.new[[paste0("D_F_",age[k])]] <- rowSums(sim.result[,grepl(paste0("(?=.*D_F)(?=.*",age[k],")"),colnames(sim.result), perl=TRUE)])
    }
    
    sim.result <- sim.result.new
    
    if(normalize){
      sim.result[,-1] <- sim.result[,-1]/rowSums(sim.result[,-1])
    }
  } else {
    if(normalize){
      sim.result[,-1] <- sim.result[,-1]/rowSums(sim.result[,-1])
    }
  }
  
  return(sim.result)
  
}

### plots the proportion of individuals in each infection state across ages for a particular simulation state (one row, minus time, of a simulation output)
# This function replaces the old plot.sim.state() function
plot.infstateprops.ages <- function(sim.state, par.cond = list(age_classes =age.classes[1:n.age.classes], nmb_act_classes = nJ, nmb_cantons = 1),
                                    aggregatebyAC = TRUE, y.lim = c(0, 1),
                                    path = NULL, 
                                    print_to_file = FALSE,
                                    file.name = "plots/plot_infstateprops_ages.pdf", 
                                    width = 6, 
                                    height = 6, 
                                    cex = 1.5,
                                    plottype = "o"){
  
  ### sim.end.state should not contain time information in the first entry of the vector
  
  #### plot an entry of the results vector
  tensors <- initvec.to.tensor.cpp(sim.state, par.cond = par.cond)
  
  U.mod <- tensors$U
  S.mod <- tensors$S
  I.mod <- tensors$I
  R.mod <- tensors$R
  V.mod <- tensors$V
  
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(U)$age))
  
  age.range <- range(ages)
  age.range <- c(10,65)
  scaling.factor <- 0.8
  x.point <- floor((age.range[1]+age.range[2])/2)
  
  ### separate the two genders into different panels
  par(mfrow = c(5,2), oma = c (5, 6, 1, 1)) #oma determines the outer margins
  par(mar=c(1,1,1,1)) # number of lines of margin to be specified on the four sides of the plot
  
  S.marg <- margin.tensor(S.mod,i=2) #activity classes are aggregated here
  I.marg <- margin.tensor(I.mod,i=2)
  R.marg <- margin.tensor(R.mod,i=2)
  V.marg <- margin.tensor(V.mod,i=2)
  
  #total number of individuals of sex g and age a
  tot <- U.mod+S.marg+I.marg+R.marg+V.marg
  
  if (aggregatebyAC){ #if aggregated by activity class
    
    ### prepare the
    if(print_to_file){
      pdf(file = paste(path, file.name, sep = ""), width = width, height = height)
    }
    
    #proportion of individuals of sex g and age a in different infection classes
    
    U.mod <- U.mod / tot #proportion in infection class over all infection classes, given specific sex and age
    S.marg <- S.marg / tot
    I.marg <- I.marg / tot
    R.marg <- R.marg / tot
    V.marg <- V.marg / tot
    
    plot(ages, U.mod["F",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "virgin females")
    plot(ages, U.mod["M",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "virgin males")
    
    plot(ages, S.marg["F",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "susceptible females")
    plot(ages, S.marg["M",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "susceptible males")
    
    plot(ages, I.marg["F",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "infected females")
    plot(ages, I.marg["M",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "infected males")
    
    plot(ages, R.marg["F",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "recovered females")
    plot(ages, R.marg["M",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "recovered males")
    
    plot(ages, V.marg["F",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "vaccinated females")
    plot(ages, V.marg["M",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "vaccinated males")
    
  } else { #if not aggregated by activity class
    
    ### prepare the
    if(print_to_file){
      pdf(file = paste(path, file.name, sep = ""), width = width, height = height)
    }
    
    U.mod <- U.mod / tot
    plot(ages, U.mod["F",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "virgin females")
    plot(ages, U.mod["M",], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "virgin males")
    
    activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(S.mod)$j))
    cols = colours(length(activity.classes))
    
    #susceptible females
    
    S.acloop <- S.mod[,1,] / tot
    plot(ages, S.acloop[2,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "susceptible females")
    for(j in activity.classes[-1]){
      S.acloop <- S.mod[,j,] / tot
      lines(ages, S.acloop[2,], type = plottype, col = cols[j], lwd = cex)
    }
    
    #susceptible males
    S.acloop <- S.mod[,1,] / tot
    plot(ages, S.acloop[1,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "susceptible males")
    for(j in activity.classes[-1]){
      S.acloop <- S.mod[,j,] / tot
      lines(ages, S.acloop[1,], type = plottype, col = cols[j], lwd = cex)
    }
    
    #infected females
    I.acloop <- I.mod[,1,] / tot
    plot(ages, I.acloop[2,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "infected females")
    for(j in activity.classes[-1]){
      I.acloop <- I.mod[,j,] / tot
      lines(ages, I.acloop[2,], type = plottype, col = cols[j], lwd = cex)
    }
    
    #infected  males
    I.acloop <- I.mod[,1,] / tot
    plot(ages, I.acloop[1,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "infected males")
    for(j in activity.classes[-1]){
      I.acloop <- I.mod[,j,] / tot
      lines(ages, I.acloop[1,], type = plottype, col = cols[j], lwd = cex)
    }
    
    #recovered females
    R.acloop <- R.mod[,1,] / tot
    plot(ages, R.acloop[2,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "recovered females")
    for(j in activity.classes[-1]){
      R.acloop <- R.mod[,j,] / tot
      lines(ages, R.acloop[2,], type = plottype, col = cols[j], lwd = cex)
    }
    
    #recovered  males
    R.acloop <- R.mod[,1,] / tot
    plot(ages, R.acloop[1,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "recovered males")
    for(j in activity.classes[-1]){
      R.acloop <- R.mod[,j,] / tot
      lines(ages, R.acloop[1,], type = plottype, col = cols[j], lwd = cex)
    }
    
    #vaccinated females
    V.acloop <- V.mod[,1,] / tot
    plot(ages, V.acloop[2,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "vaccinated females")
    for(j in activity.classes[-1]){
      V.acloop <- V.mod[,j,] / tot
      lines(ages, V.acloop[2,], type = plottype, col = cols[j], lwd = cex)
    }
    
    #vaccinated  males
    V.acloop <- V.mod[,1,] / tot
    plot(ages, V.acloop[1,], ylim = y.lim, cex = cex, type=plottype)
    text(x = x.point, y = scaling.factor*y.lim[2], labels = "vaccinated males")
    for(j in activity.classes[-1]){
      V.acloop <- V.mod[,j,] / tot
      lines(ages, V.acloop[1,], type = plottype, col = cols[j], lwd = cex)
    }
    
  }
  if (print_to_file){
    dev.off()
  }
}

# last.sim.state <- sim.res[dim(sim.res)[1],2:dim(sim.res)[2]]
# plot.infstateprops.ages(last.sim.state, aggregatebyAC = T)
# plot.infstateprops.ages(last.sim.state, aggregatebyAC = F)



### This function plots the time course of simulation outcomes of the 
### generic sti model. It loops through all compartments. 
plot.sim.out <- function(	sim.result, combine.ages = TRUE, combine.AC = FALSE,
                          n_age = length(age.classes)-1, n_ac = nJ, normalize = TRUE,
                          y.lim = NULL,
                          path = NULL, 
                          print_to_file = FALSE,
                          file.name = "plots/genericsti_sim_run.pdf", 
                          main = "",
                          width = 6, 
                          height = 6, 
                          do_legend = TRUE, 
                          cex = 1.5,
                          returnstats=FALSE){
  
  ### prepare the print to file
  if(print_to_file){
    pdf(file = paste(path, file.name, sep = ""), width = width, height = height)
  }
  
  sim.result <- aggregate_simresult( sim.result, comb.ages = combine.ages, comb.AC = combine.AC, normalize = normalize)
  
  #calculate number of rows in data (how many simulation steps)
  L <- dim(sim.result)[1]
  n <- numeric(L)
  
  ## number of types
  nmb.types <- dim(sim.result)[2]-1
  ## number of columns
  nmb.cols <- dim(sim.result)[2]
  
  ### determine plotting range for x
  #x.lim <- c(0, sim.result[L,1])
  #x.lim <- c(0, sim.result[L,1]) * ar #JS270116
  sim.time <- sim.result[,1]
  x.lim <- c(0,max(sim.time))
  
  ### build the sum of all possible "compartments", to check whether they remain constant
  for (i in 1:L){
    n[i] <- sum(sim.result[i, 2:nmb.cols])
  }
  
  ### determine plotting range for y
  #y_lim <- c(0,max(n))
  if(is.null(y.lim)){
    y.lim <- c(0,ceiling(max(sim.result[,-1])))
  }
  
  ## open the plotting device: by plotting two extreme points which span up the plotting area 
  plot(range(sim.time),  #range(sim.result[,1]) , JS270116
       y.lim, 
       main = main, 
       type = "n", 
       ylim = y.lim,
       xlim = x.lim,
       axes = F, 
       xlab = "Time [years]", 
       ylab = "Population Density", 
       cex.axis = cex, 
       cex.lab = cex
  )
  
  leg.txt <- character(nmb.types)
  colrs <- sim.cols(nmb.types)
  
  for (k in 1:nmb.types){
    lines(sim.time, #sim.result[,1], #JS270116
          sim.result[,k+1], 
          col = colrs[k], #sim.cols(k), 
          lty = 1, 
          ylim = y.lim, 
          xlim = x.lim,
          lwd = 2)
    #leg.txt[k+1] = paste(legend.names(k), sep = "")
  }
  
  
  ## add the total population
  lines(sim.time, 
        n, 
        col = 1, 
        lty = 1, 
        ylim = y.lim, 
        xlim = x.lim)
  
  ## add this information to the 
  leg.txt[1] = "dunno"
  #colrs[1]  = 1
  
  scal.fac <- 0.8
  legend_position <- c(scal.fac*max(x.lim), scal.fac*max(y.lim))
  
  if(do_legend){
    legend(legend_position[1],legend_position[2], 
           colnames(sim.result)[-1], 
           col = colrs, 
           lty = rep(1,nmb.types), ## line type
           cex = .7 )
  }
  
  axis(1, cex.axis = cex); axis(2, cex.axis = cex)
  
  if (print_to_file){
    dev.off()
  }
  if(returnstats){
    return(tail(sim.result,n=1))
  }
}

### plots the proportion of one state of a particular simulation state (one row, minus time, of a simulation output) among all other states
### for different age groups. Default: plot proportion of infecteds (prevalence) among all states
#OLD NAME (before 200616): plot.prevalence.ages(): 
plot.props.ages <- function(sim.state_OR_epidata,
                            infect.comp = "I",
                            use_epidata=F, #if this is false, the prevalence is computed from sim.state. If TRUE, the prevalence is computed from epidata
                            par.cond = list(age_classes =age.classes[1:n.age.classes], nmb_act_classes = nJ),
                            aggregatebyAC = TRUE,
                            path = NULL, 
                            print_to_file = FALSE,
                            file.name = "plots/genericsti_sim_endstate.pdf", 
                            width = 6, 
                            height = 6, 
                            cex = 1.5,
                            plottype = "o",
                            additionalmain = "",
                            returnstats=T,
                            doplot=T){
  
  
  #### plot an entry of the results vector
  if (use_epidata){ 
    #use epidata to compute the prevalence
    #epidata is a tensor with one dimension being the numerator (number positives) and another dimension being the denominator (number tested)
    
    if (aggregatebyAC){
      sim.state_OR_epidata <- margin.tensor(sim.state_OR_epidata,i=2)
      prevalence <- sim.state_OR_epidata[,,1]/sim.state_OR_epidata[,,2]
      dimnames(prevalence)<-list(sex=c("M","F"), age=paste0('age', age.class.groups))
    } else {
    
      prevalence <- sim.state_OR_epidata[,,,1]/sim.state_OR_epidata[,,,2]
      dimnames(prevalence)<-list(sex=c("M","F"), j=paste0('j',1:nJ_data),age=paste0('age', age.class.groups))
    }
    
  } else { 
    
    ### sim.end.state should not contain time information in the first entry of the vector
    
    #### plot an entry of the results vector
    tensors <- initvec.to.tensor.cpp(sim.state_OR_epidata, par.cond = par.cond)
    
    # B.mod <- tensors$B
    U.mod <- tensors$U
    S.mod <- tensors$S
    I_A.mod <- tensors$I_A
    I_S.mod <- tensors$I_S
    R.mod <- tensors$R
    
    if (aggregatebyAC){
      
      S.mod <- margin.tensor(S.mod,i=2)
      I_A.mod <- margin.tensor(I_A.mod,i=2)
      I_S.mod <- margin.tensor(I_S.mod,i=2)
      R.mod <- margin.tensor(R.mod,i=2)

    }
    
    if (infect.comp=="I"){
      prevalence <- (I_A.mod + I_S.mod) / (S.mod + I_A.mod + I_S.mod + R.mod)
    } else {
      prevalence <- eval(parse( text=paste0(infect.comp,".mod"))) / (S.mod + I_A.mod + I_S.mod + R.mod)
    }
    
    
  }
  
  ### prepare the
  if(print_to_file & doplot){
    pdf(file = paste(path, file.name, sep = ""), width = width, height = height)
  }
  
  activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(prevalence)$j))
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(prevalence)$age))
  age.range <- range(ages)
  #age.range <- c(10,65)
  
  # DO THE PLOTS
  if (doplot){
    #par(mfrow = c(1,2), oma = c(2, 2, 1, 1))
    ymax <- max(prevalence) + .2 * max(prevalence)
    plot(min(age.range), 0, ylim = c(0,ymax) , xlim = age.range, type = "n", main = paste0("females",additionalmain), xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
  }
  
  if (aggregatebyAC){
    if (doplot){
      lines(ages, prevalence["F",], type = plottype, col = "red", lwd = cex)
    }
    
    if(returnstats){
      maxprev.F <- max(prevalence["F",])
      age.maxprev.F <- ages[prevalence["F",]==maxprev.F]
    }
    
  } else{
    if (doplot){
      cols = colours(length(activity.classes))
      for(j in activity.classes){
        lines(ages, prevalence["F",j,], type = plottype, col = cols[j], lwd = cex)
      }
      legend("topright",c("low","high") ,lty=rep(1,length(activity.classes)), col = cols,bty = "n")
      
    }
    if(returnstats){
      maxprev.F1 <- max(prevalence["F",1,])
      age.maxprev.F1 <- ages[prevalence["F",1,]==maxprev.F1]
      maxprev.F2 <- max(prevalence["F",2,])
      age.maxprev.F2 <- ages[prevalence["F",2,]==maxprev.F2]
    }
  }
  
  if (doplot){
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
    
    ### plotting males 
    plot(min(age.range), 0, ylim = c(0,ymax) , xlim = age.range, type = "n", main = paste0("males",additionalmain), xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
  }
  
  if (aggregatebyAC){
    if (doplot){
      lines(ages, prevalence["M",], type = plottype, col = "red", lwd = cex)
      
    }
    
    if(returnstats){
      maxprev.M <- max(prevalence["M",])
      age.maxprev.M <- ages[prevalence["M",]==maxprev.M]
    }
    
  } else{
    if (doplot){
      cols = colours(length(activity.classes))
      for(j in activity.classes){
        lines(ages, prevalence["M",j,], type = plottype, col = cols[j], lwd = cex)
        legend("topright",c("low","high") ,lty=rep(1,length(activity.classes)), col = cols,bty = "n")
      }
    }
    
    if(returnstats){
      maxprev.M1 <- max(prevalence["M",1,])
      age.maxprev.M1 <- ages[prevalence["M",1,]==maxprev.M1]
      maxprev.M2 <- max(prevalence["M",2,])
      age.maxprev.M2 <- ages[prevalence["M",2,]==maxprev.M2]
    }
    
  }
  
  if (doplot){
    axis(1, cex.axis = cex, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = cex, las = 1)
    
    ### axis labels
    mtext("Frequency", side=2, outer=T, line= 0, cex = cex)
    mtext("Age [years]", side=1, outer=T, line= 0, cex = cex)
  }
  
  if(returnstats){
    if (aggregatebyAC){
      return(list(prevalence=prevalence, maxprev.F=maxprev.F, age.maxprev.F=age.maxprev.F, maxprev.M=maxprev.M, age.maxprev.M=age.maxprev.M ))
    } else {
      return(list(prevalence=prevalence, maxprev.F1=maxprev.F1, age.maxprev.F1=age.maxprev.F1, maxprev.F2=maxprev.F2, age.maxprev.F2=age.maxprev.F2, 
                  maxprev.M1=maxprev.M1, age.maxprev.M1=age.maxprev.M1 , maxprev.M2=maxprev.M2, age.maxprev.M2=age.maxprev.M2 ))
    }
  }
  
  if (print_to_file & doplot){
    dev.off()
  }
}

### plots the proportion of one state of a particular simulation state (one row, minus time, of a simulation output) among all other states
### for different age groups (1 gender). Default: plot proportion of infecteds (prevalence) among all states
#OLD NAME (before 200616): plot.prevalence.ages(): 
plot.props.ages.onesex <- function(sim.state_OR_epidata,
                                        infect.comp = "I",
                                        use_epidata=F, #if this is false, the prevalence is computed from sim.state. If TRUE, the prevalence is computed from epidata
                                        par.cond = list(age_classes =age.class.groups, nmb_act_classes = nJ),# nmb_cantons = 1),
                                        aggregatebyAC = TRUE,
                                        path = NULL, 
                                        print_to_file = FALSE,
                                        file.name = "plots/genericsti_sim_endstate.pdf", 
                                        width = 6, 
                                        height = 6, 
                                        mycex = 1.5,
                                        my.lwd = 1.5,
                                        mycex.axis=1,
                                        plottype = "o",
                                        plotsymbol = 1,
                                        returnstats=FALSE,
                                        mysex="M",
                                        mycol="red",
                                        mylty=1,
                                        ymax=NULL,
                                        plotlineonly=FALSE){
  

  #### plot an entry of the results vector
  if (use_epidata){ 
    #use epidata to compute the prevalence
    #epidata is a tensor with one dimension being the numerator (number positives) and another dimension being the denominator (number tested)
    
    if (aggregatebyAC){
      #sim.state_OR_epidata <- margin.tensor(sim.state_OR_epidata,i=2)
      prevalence <- sim.state_OR_epidata[,,1]/sim.state_OR_epidata[,,2]
      dimnames(prevalence)<-list(sex=c("M","F"), age=paste0('age', age.class.groups))
    } else {
      
      prevalence <- sim.state_OR_epidata[,,,1]/sim.state_OR_epidata[,,,2]
      dimnames(prevalence)<-list(sex=c("M","F"), j=paste0('j',1:nJ_data),age=paste0('age', age.class.groups))
    }
    
  } else { 
    # use sim.state to compute the prevalence
    # sim.state is a vector, obtained from the ode solution
    # sim.end.state should not contain time information in the first entry of the vector
    
    tensors <- initvec.to.tensor.cpp(sim.state_OR_epidata, par.cond = par.cond)
    
    U.mod <- tensors$U
    S.mod <- tensors$S
    I.mod <- tensors$I
    R.mod <- tensors$R
    V.mod <- tensors$V
    
    if (aggregatebyAC){
      
      S.mod <- margin.tensor(S.mod,i=2)
      I.mod <- margin.tensor(I.mod,i=2)
      R.mod <- margin.tensor(R.mod,i=2)
      V.mod <- margin.tensor(V.mod,i=2)
    }
    
    ### calculate the tensor of interest
    prevalence <- eval(parse( text=paste0(infect.comp,".mod"))) / (S.mod + I.mod + R.mod + V.mod)
  }
  
  ### prepare the
  if(print_to_file){
    pdf(file = paste(path, file.name, sep = ""), width = width, height = height)
  }
  
  activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(prevalence)$j))
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(prevalence)$age))
  age.range <- range(ages)
  #   age.range <- c(10,65)
  #   age.range <- c(16,60)
  
  ### plotting
  if (is.null(ymax)){
    ymax <- max(prevalence) + .2 * max(prevalence)
  }

  if (!plotlineonly){
    plot(min(age.range), 0, ylim = c(0,ymax) , xlim = age.range, type = "n", xlab = "", ylab = "", lwd=my.lwd, cex = mycex, axes = FALSE, cex.main = mycex, pch = plotsymbol,lty=mylty)
    }
  
  if (aggregatebyAC){
    lines(ages, prevalence[mysex,], type = plottype, col = mycol, lwd=my.lwd, pch = plotsymbol,lty=mylty)
    
    if(returnstats){
      maxprev <- max(prevalence[mysex,])
      age.maxprev <- ages[prevalence[mysex,]==maxprev]
    }
    
  } else{
    cols = colours(length(activity.classes))
    for(j in activity.classes){
      lines(ages, prevalence[mysex,j,], type = plottype, col = cols[j], lwd=my.lwd, pch = plotsymbol,lty=mylty)
    }
    if(returnstats){
      maxprev1 <- max(prevalence[mysex,1,])
      age.maxprev1 <- ages[prevalence[mysex,1,]==maxprev1]
      maxprev2 <- max(prevalence[mysex,2,])
      age.maxprev2 <- ages[prevalence[mysex,2,]==maxprev2]
    }
  }
  
  if (!plotlineonly){
    axis(1, cex.axis = mycex.axis, mgp  = c(3,1.5,0)) 
    axis(2, cex.axis = mycex.axis, las = 1)
    
    ### axis labels
    #mtext("Prevalence", side=2, outer=T, line= 0, cex = mycex)
    #mtext("Age [years]", side=1, outer=T, line= 0, cex = mycex)
  }
  
  if(returnstats){
    if (aggregatebyAC){
      return(list(prev=prevalence[mysex,], maxprev=maxprev, age.maxprev=age.maxprev))
    } else {
      return(list(prev=prevalence[mysex,], maxprev1=maxprev1, age.maxprev1=age.maxprev1, maxprev2=maxprev2, age.maxprev2=age.maxprev2 ))
    }
  }
  
  if (print_to_file){
    dev.off()
  }
}

plotprevs <- function(simulation, epidata=NULL, ymax=NULL){
  
  # function to plot prevalences in 2000 and 2011
  # epidata should be given to the function as a list with elements epidata$time_0 and epidata$time_11
  
  ### MODEL COMPUTED PREVALENCES
  
  ###  2000 #########################################################################
  prescreen.sim.state <- simulation[burnintime,2:dim(simulation)[2]]
  tensors <- initvec.to.tensor.cpp(prescreen.sim.state, par.cond = list(age_classes =age.classes[1:n.age.classes], nmb_act_classes = nJ))
  
  U.mod <- tensors$U
  S.mod <- tensors$S
  I_A.mod <- tensors$I_A
  I_S.mod <- tensors$I_S
  R.mod <- tensors$R
  I_A.mod2 <- tensors$I_A2
  I_S.mod2 <- tensors$I_S2
  
  S.mod <- margin.tensor(S.mod,i=2)
  I_A.mod <- margin.tensor(I_A.mod,i=2)
  I_S.mod <- margin.tensor(I_S.mod,i=2)
  R.mod <- margin.tensor(R.mod,i=2)
  I_A.mod2 <- margin.tensor(I_A.mod2,i=2)
  I_S.mod2 <- margin.tensor(I_S.mod2,i=2)
  
  
  
  prevalence <- (I_A.mod + I_S.mod + I_A.mod2 + I_S.mod2) / (S.mod + I_A.mod + I_S.mod + R.mod + I_A.mod2 + I_S.mod2) # model prevalence
  
  activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(prevalence)$j))
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(prevalence)$age))
  age.range <- range(ages)
  
  par(mfrow=c(2,2))
  if (is.null(ymax)){
    ymax <- max(prevalence) + .2 * max(prevalence)
  }
  cex = 1.5
  plottype = "o"
  plot(min(age.range), 0,ylim=c(0,ymax) , xlim = age.range, type = "n", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
  lines(ages, prevalence["M",], type = plottype, col = "red", lwd = cex)
  
  if (!is.null(epidata)){
    library(binom)
    CI2M <- binom.confint(as.numeric(epidata$time_0[1,,1]),as.numeric(epidata$time_0[1,,2]), conf.level=.95,method="wilson")
    errbar(age.class.groups, epidata$time_0[1,,1] / epidata$time_0[1,,2], CI2M$upper, CI2M$lower, add=T, pch=1, cap=.015)
  }
  title(main = "male prevalence - 2000")
  
  axis(1, at=age.class.groups, labels = FALSE)
  text(age.class.groups, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = paste0(age.classes[1:n.age.classes], "-",age.classes[2:(n.age.classes+1)]-1), 
       srt = 45, pos = 1, xpd = TRUE)
  axis(2, cex.axis = 1, las = 1)
  
  ### plotting females 
  plot(min(age.range), 0,ylim=c(0,ymax) , xlim = age.range, type = "n", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
  lines(ages, prevalence["F",], type = plottype, col = "red", lwd = cex)
  if (!is.null(epidata)){
    CI2F <- binom.confint(as.numeric(epidata$time_0[2,,1]),as.numeric(epidata$time_0[2,,2]), conf.level=.95,method="wilson")
    errbar(age.class.groups, epidata$time_0[2,,1] / epidata$time_0[2,,2], CI2F$upper, CI2F$lower, add=T, pch=1, cap=.015)
  }
  title(main = "female prevalence - 2000")
  
  axis(1, at=age.class.groups, labels = FALSE)
  text(age.class.groups, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = paste0(age.classes[1:n.age.classes], "-",age.classes[2:(n.age.classes+1)]-1), 
       srt = 45, pos = 1, xpd = TRUE) 
  axis(2, cex.axis = 1, las = 1)
  
  ### axis labels
  mtext("Frequency", side=2, outer=T, line= 0, cex = cex)
  mtext("Age [years]", side=1, outer=T, line= 0, cex = cex)
  
  ###  2012 #########################################################################
  last.sim.state <- simulation[max.simtime,2:dim(simulation)[2]]
  tensors <- initvec.to.tensor.cpp(last.sim.state, par.cond = list(age_classes =age.classes[1:n.age.classes], nmb_act_classes = nJ))
  
  # B.mod <- tensors$B
  U.mod <- tensors$U
  S.mod <- tensors$S
  I_A.mod <- tensors$I_A
  I_S.mod <- tensors$I_S
  R.mod <- tensors$R
  
  S.mod <- margin.tensor(S.mod,i=2)
  I_A.mod <- margin.tensor(I_A.mod,i=2)
  I_S.mod <- margin.tensor(I_S.mod,i=2)
  R.mod <- margin.tensor(R.mod,i=2)
  
  prevalence <- (I_A.mod + I_S.mod) / (S.mod + I_A.mod + I_S.mod + R.mod)
  
  activity.classes <- as.numeric(gsub("j(.*)", '\\1' ,dimnames(prevalence)$j))
  ages <- as.numeric(gsub("age(.*)", '\\1' ,dimnames(prevalence)$age))
  age.range <- range(ages)
  
  if (is.null(ymax)){
    ymax <- max(prevalence) + .2 * max(prevalence)
  }
  cex = 1.5
  plottype = "o"
  plot(min(age.range), 0,ylim=c(0,ymax) , xlim = age.range, type = "n", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
  lines(ages, prevalence["M",], type = plottype, col = "red", lwd = cex)
  if (!is.null(epidata)){
    CI3M <- binom.confint(as.numeric(epidata$time_11[1,,1]),as.numeric(epidata$time_11[1,,2]), conf.level=.95,method="wilson")
    errbar(age.class.groups, epidata$time_11[1,,1] / epidata$time_11[1,,2], CI3M$upper, CI3M$lower, add=T, pch=1, cap=.015)
  }
  title(main = "male prevalence - 2012")
  
  axis(1, at=age.class.groups, labels = FALSE)
  text(age.class.groups, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = paste0(age.classes[1:n.age.classes], "-",age.classes[2:(n.age.classes+1)]-1), 
       srt = 45, pos = 1, xpd = TRUE) 
  axis(2, cex.axis = 1, las = 1)
  
  ### plotting males 
  plot(min(age.range), 0,ylim=c(0,ymax) , xlim = age.range, type = "n", xlab = "", ylab = "", cex = cex, axes = FALSE, cex.main = cex)
  lines(ages, prevalence["F",], type = plottype, col = "red", lwd = cex)
  if (!is.null(epidata)){
    CI3F <- binom.confint(as.numeric(epidata$time_11[2,,1]),as.numeric(epidata$time_11[2,,2]), conf.level=.95,method="wilson")
    errbar(age.class.groups, epidata$time_11[2,,1] / epidata$time_11[2,,2], CI3F$upper, CI3F$lower, add=T, pch=1, cap=.015)
  }
  title(main = "female prevalence - 2012")
  
  axis(1, at=age.class.groups, labels = FALSE)
  text(age.class.groups, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = paste0(age.classes[1:n.age.classes], "-",age.classes[2:(n.age.classes+1)]-1), 
       srt = 45, pos = 1, xpd = TRUE) 
  axis(2, cex.axis = 1, las = 1)
  
  ### axis labels
  mtext("Frequency", side=2, outer=T, line= 0, cex = cex)
  mtext("Age [years]", side=1, outer=T, line= 0, cex = cex)
}

plotdiagdata <- function(simulation){
  
  require("RColorBrewer")
  mycol <- brewer.pal(9,"Set1")
  mycol <- c(mycol[9],mycol[1:8])
  
  age.catdata <- c(15,20,25,35,45)
  years.in.categories <- to.tensor(NA, dims = list(age=paste0('age', age.classes[1:(n.age.classes)]) ,screencatage=paste0('diagcatage', age.catdata[1:(length(age.catdata)-1)])))
  for (k in 1:n.age.classes){
    for (y in 1:(length(age.catdata)-1)){
      years.in.categories[k,y] <- length(intersect( seq(age.classes[k],age.classes[k+1]-1), seq(age.catdata[y],age.catdata[y+1]-1) ))
    }
  }
  norm.years.in.categories <- years.in.categories / margin.tensor(years.in.categories, i=1) # fraction of time that individuals in age.classes[k] are in the screencatage[j] class
  
  # tests
  source(file.path(path_in,"CT_screendata.R"))
  uptake.minmax <- get.CT.screendata(nat,age.catdata)
  names_uptake <- dimnames(uptake.minmax)
  names_uptake$screening_y <- names_uptake$screening_y[1:12]
  uptake.minmax <- uptake.minmax[,,,1:12]
  dimnames(uptake.minmax) <- names_uptake
  
  size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416
  
  par(mfrow=c(2,2))
  plot(1, type="n", xaxt = "n", xlab="years", ylab="# tests pp", main="men total", xlim=range(burnintime:max.simtime), ylim=c(0, 0.25))
  axis(1, at=burnintime:max.simtime, labels = FALSE)
  text(burnintime:max.simtime, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = 2000:2011, 
       srt = 45, pos = 1, xpd = TRUE) 
  for (k in 1:(length(age.catdata)-1)){
    points(burnintime:max.simtime,(uptake.minmax[1,k,1,]+uptake.minmax[1,k,2,])/2, col=mycol[k], pch=16)
    segments(x0=burnintime:max.simtime,y0=uptake.minmax[1,k,1,], x1=burnintime:max.simtime, y1=uptake.minmax[1,k,2,], col=mycol[k])
  }
  legend("topleft",paste0(age.classes[1:n.age.classes][-2], "-",age.classes[2:(n.age.classes+1)][-1]), lty=1, col=mycol, bty = "n")
  
  plot(1, type="n", xaxt = "n", xlab="years", ylab="# tests pp", main="women total", xlim=range(burnintime:max.simtime), ylim=c(0, 0.45))
  axis(1, at=burnintime:max.simtime, labels = FALSE)
  text(burnintime:max.simtime, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = 2000:2011, 
       srt = 45, pos = 1, xpd = TRUE) 
  for (k in 1:(length(age.catdata)-1)){
    points(burnintime:max.simtime,(uptake.minmax[2,k,1,]+uptake.minmax[2,k,2,])/2, col=mycol[k], pch=16)
    segments(x0=burnintime:max.simtime,y0=uptake.minmax[2,k,1,], x1=burnintime:max.simtime, y1=uptake.minmax[2,k,2,], col=mycol[k])
  }
  legend("topleft",paste0(age.classes[1:n.age.classes][-2], "-",age.classes[2:(n.age.classes+1)][-1]), lty=1, col=mycol, bty = "n")
  
  # diagnoses
  source(file.path(path_in,"CT_diagdata.R"))
  epidata_diag <- get.CT.diagdata(nat,age.catdata)
  names_diag <- dimnames(epidata_diag)
  names_diag$screening_y <- names_diag$screening_y[1:12]
  epidata_diag <- epidata_diag[,,,1:12]
  dimnames(epidata_diag) <- names_diag
  
  diagtot_M <- simulation[,grep("D_M_j1_age",names(simulation))] + simulation[,grep("D_M_j2_age",names(simulation))]
  
  size.age.classes <- to.tensor(diff(age.classes), dims=list(age=paste0('age', age.classes[1:n.age.classes]))) #JS040416
  
  #rename the row names data frame
  diagtot_M <- to.tensor(unlist(diagtot_M), dims= list(sim=paste0('sim',1:dim(diagtot_M)[1]), age=paste0('age', age.classes[1:(n.age.classes)]) ))
  # convert to per capita number of diagnoses in the same age groups as epidata_diag
  diagtot_M <- diagtot_M / size.age.classes * norm.years.in.categories
  diagtot_M <- margin.tensor(diagtot_M, i=2) # per capita number of diagnoses in age classes in data
  dimnames(diagtot_M) <- list(sim=paste0('sim',1:dim(diagtot_M)[1]), age=paste0('age', age.classes[1:(n.age.classes-1)]) )
  
  maxD <- max(diagtot_M)
  plot(1, type="n", xaxt = "n", xlab="years", ylab="# diagnosed cases", main="men total", xlim=range(burnintime:max.simtime), ylim=c(0, 0.025))
  axis(1, at=burnintime:max.simtime, labels = FALSE)
  text(burnintime:max.simtime, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = 2000:2011, 
       srt = 45, pos = 1, xpd = TRUE) 
  for (k in 1:(length(age.catdata)-1)){
    lines(burnintime:max.simtime,diagtot_M[burnintime:max.simtime,k], col=mycol[k])
    points(burnintime:max.simtime,((epidata_diag[1,k,1,] + epidata_diag[1,k,2,])/2)/100000, col=mycol[k], pch=16)
    segments(x0=burnintime:max.simtime,y0=epidata_diag[1,k,1,]/100000, x1=burnintime:max.simtime, y1=epidata_diag[1,k,2,]/100000, col=mycol[k])
  }
  legend("topleft",paste0(age.classes[1:n.age.classes][-2], "-",age.classes[2:(n.age.classes+1)][-1]), lty=1, col=mycol, bty = "n")
  
  diagtot_F <- simulation[,grep("D_F_j1_age",names(simulation))] + simulation[,grep("D_F_j2_age",names(simulation))]
  
  diagtot_F <- to.tensor(unlist(diagtot_F), dims = list(sim=paste0('sim',1:dim(diagtot_F)[1]), age=paste0('age', age.classes[1:(n.age.classes)]) ))
  diagtot_F <- diagtot_F / size.age.classes * norm.years.in.categories
  diagtot_F <- margin.tensor(diagtot_F, i=2) # per capita number of diagnoses in age classes in data
  dimnames(diagtot_F) <- list(sim=paste0('sim',1:dim(diagtot_F)[1]), age=paste0('age', age.classes[1:(n.age.classes-1)]) )
  
  maxD <- max(diagtot_F)
  plot(1, type="n", xaxt = "n", xlab="years", ylab="# diagnosed cases", main="women total", xlim=range(burnintime:max.simtime), ylim=c(0, 0.04))
  axis(1, at=burnintime:max.simtime, labels = FALSE)
  text(burnintime:max.simtime, par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]), labels = 2000:2011, 
       srt = 45, pos = 1, xpd = TRUE) 
  for (k in 1:(length(age.catdata)-1)){
    lines(burnintime:max.simtime,diagtot_F[burnintime:max.simtime,k], col=mycol[k])
    points(burnintime:max.simtime,((epidata_diag[2,k,1,] + epidata_diag[2,k,2,])/2)/100000, col=mycol[k], pch=16)
    segments(x0=burnintime:max.simtime,y0=epidata_diag[2,k,1,]/100000, x1=burnintime:max.simtime, y1=epidata_diag[2,k,2,]/100000, col=mycol[k])
  }
  legend("topleft",paste0(age.classes[1:n.age.classes][-2], "-",age.classes[2:(n.age.classes+1)][-1]), lty=1, col=mycol, bty = "n")
  
}

plot.pos.rate <- function(simulation){
  
  source(file.path(path_in,"CT_screendata.R"))
  uptake.minmax <- get.CT.screendata(nat,age.classes)
  
  agecategories <- age.classes[1:n.age.classes]
  years <- 2000+1:modelled_NCSP_years
  
  posrate.mean <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.classes[1:n.age.classes]),screening_y = years))
  posrate.low <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.classes[1:n.age.classes]),screening_y = years))
  posrate.upp <- to.tensor(NA, dims = list(sex=c("M","F"), age=paste0('age', age.classes[1:n.age.classes]),screening_y = years))
  
  for (y in 1:modelled_NCSP_years){
    
    tensors <- initvec.to.tensor.cpp(simulation[burnintime+y,-1], par.cond = list(age_classes =age.classes[1:n.age.classes], nmb_act_classes = nJ))
    
    # per capita number of diagnoses
    Npostest.y <- margin.tensor(tensors$D,i=2) / size.age.classes
    
    # per capita expected number of tests
    Ntest.y.mean <- to.tensor(as.vector((uptake.minmax[,,1,y] + uptake.minmax[,,2,y])/2), dims = list(sex=c("M","F"), age=paste0('age', age.classes[1:n.age.classes])))
    Ntest.y.low <- to.tensor(as.vector(uptake.minmax[,,1,y]), dims = list(sex=c("M","F"), age=paste0('age', age.classes[1:n.age.classes])))
    Ntest.y.upp <- to.tensor(as.vector(uptake.minmax[,,2,y]), dims = list(sex=c("M","F"), age=paste0('age', age.classes[1:n.age.classes])))
    
    posrate.mean[,,y] <- Npostest.y/Ntest.y.mean
    posrate.low[,,y] <- Npostest.y/Ntest.y.upp
    posrate.upp[,,y] <- Npostest.y/Ntest.y.low
    
  }
  
  par(mfrow=c(n.age.classes,2))
  for (a in 1:n.age.classes){
    plot(years, posrate.mean[1,a,], type="p",xlab="years",ylab=paste0('Pos. rate M, age', age.classes[a]), pch=16)
    segments(x0=years,y0=posrate.low[1,a,], x1=years, y1=posrate.upp[1,a,], col="black")
    
    plot(years, posrate.mean[2,a,], type="p", col="red",xlab="years",ylab=paste0('Pos. rate F, age', age.classes[a]), pch=16)
    segments(x0=years,y0=posrate.low[2,a,], x1=years, y1=posrate.upp[2,a,], col="red")
  }
  
}