

#' @export


#######################################################################################

#### ESCV functions

# from Yu et al
escv.glmnet <- function(X,Y,k=10, nlambda=100) {
  require(glmnet)
  n = dim(X)[1]
  p = dim(X)[2]
  if (length(Y)!=n) {warning("length(Y) != nrow(X)")}

  glmnet.list = list()
  beta.list = list()
  Xb.list = list()

  #Get overall solution path

  orig_glmnet = glmnet::glmnet(X,Y, family="gaussian", alpha=1,
                                nlambda=nlambda, standardize=TRUE, intercept=TRUE)

  lambdavec = orig_glmnet$lambda
  orig_betahat = as.matrix(orig_glmnet$beta)

  #Get pseudo solutions

  grps = cut(1:n,k,labels=FALSE)[sample(n)]

  for (i in 1:k) {
    omit = which(grps==i)
    glmnet.list[[i]] = glmnet(X[-omit,,drop=FALSE],Y[-omit], family="gaussian", alpha=1, lambda=lambdavec, standardize=TRUE, intercept=TRUE)
    beta.list[[i]] = as.matrix(glmnet.list[[i]]$beta)
    Xb.list[[i]] = (X %*% beta.list[[i]])
  }

  ## CV

  l2pe = matrix(0,length(lambdavec),k)

  for (i in 1:k) {
    omit = which(grps==i)
    pe = Y[omit] - X[omit,,drop=FALSE] %*% beta.list[[i]]
    l2pe[,i] = apply(pe^2,2,mean)
  }

  CV = apply(l2pe,1,mean)
  CV.index = which.min(CV)

  #Compute ES metric

  beta.sum = beta.list[[1]]
  for (i in 2:k) {
    beta.sum = beta.sum + beta.list[[i]]
  }
  beta.mean <- beta.sum /k
  Xb.mean <- (X %*% beta.mean)

  l2diffX.matrix = matrix(0,nrow=k,ncol=length(lambdavec))
  for (i in 1:k) {
    l2diffX.matrix[i,]=apply((Xb.list[[i]] - Xb.mean)^2,2,sum)
  }
  l2diffX <- apply(l2diffX.matrix,2,mean)

  ES = l2diffX/apply(Xb.mean^2,2,sum)

  # check for convergence
  ESgrad = ES[2:length(ES)]-ES[1:(length(ES)-1)]
  con.index = which(ESgrad < 0)
  if(length(con.index)==0) {
    con.index = 1
  } else {con.index = con.index[1]}

  # Find min after convergence before CV choice
  if(con.index > CV.index) {
    con.index = CV.index
  }

  ESCV.index = which.min(ES[con.index:CV.index]) + con.index - 1

  #browser()


  # BIC, extended BIC
  l0norm <- function(x){return(sum(x!=0))}
  df = apply(orig_betahat,2,l0norm)
  RSSMat = X %*% orig_betahat
  ncol2 = dim(RSSMat)[2]
  for (i in 1:ncol2) {
    RSSMat[,i] = Y - RSSMat[,i]
  }
  RSS = apply( (RSSMat)^2 , 2 ,sum )

  BIC = n*log(RSS) + log(n)*df
  EBIC1 = n*log(RSS) + log(n)*df + 2 * 0.5 * log(choose(p,df)) #Gamma = 0.5 for high dimensional case
  EBIC2 = n*log(RSS) + log(n)*df + 2 * 1 * log(choose(p,df)) #Gamma = 1 for ultra high dimensional case

  BIC.index = which.min(BIC)
  EBIC.index = c(which.min(EBIC1),which.min(EBIC2))

  results = list()
  results$glmnet <- orig_glmnet
  results$selindex <- c(ESCV.index,CV.index,BIC.index, EBIC.index)
  names(results$selindex) <- c("ESCV","CV","BIC","EBIC1","EBIC2")
  results$ES = ES
  results$CV = CV
  return(results)

}

# from Yu et al
#### repeated CV for v=2 repeated 5 times

rcv.glmnet <- function(X,Y, k=2, r=5, nlambda=100) {

  require(glmnet)
  n = dim(X)[1]
  p = dim(X)[2]
  if (length(Y)!=n) {warning("length(Y) != nrow(X)")}

  glmnet.list = list()
  beta.list = list()
  Xb.list = list()

  #Get overall solution path

  orig_glmnet = glmnet::glmnet(X,Y, family="gaussian", alpha=1, nlambda=nlambda, standardize=TRUE, intercept=FALSE)

  lambdavec = orig_glmnet$lambda
  orig_betahat = as.matrix(orig_glmnet$beta)

  #Get pseudo solutions
  l2per = matrix(0,length(lambdavec),k*r)

  for (ii in 1:r) {

    grps = cut(1:n,k,labels=FALSE)[sample(n)]

    for (i in 1:k) {
      omit = which(grps==i)
      glmnet.list[[i]] = glmnet(X[-omit,,drop=FALSE],Y[-omit], family="gaussian", alpha=1, lambda=lambdavec, standardize=FALSE, intercept=FALSE)
      beta.list[[i]] = as.matrix(glmnet.list[[i]]$beta)
      Xb.list[[i]] = (X %*% beta.list[[i]])
    }

    ## CV

    l2pe = matrix(0,length(lambdavec),k)

    for (i in 1:k) {
      omit = which(grps==i)
      pe = Y[omit] - X[omit,,drop=FALSE] %*% beta.list[[i]]
      l2pe[,i] = apply(pe^2,2,mean)
    }

    l2per[,((ii-1)*k+1):(ii*k)]=l2pe

  }


  rCV = apply(l2per,1,mean)
  rCV.index = which.min(rCV)

  results = list()
  results$glmnet <- orig_glmnet
  results$selindex <- rCV.index
  return(results)
}

# from Wyneken and Yang

getESCVCoef_Normal = function(xMat,yVec) {
  #browser()
  escvRes = escv.glmnet(X=scale(xMat,center=TRUE,scale=TRUE),Y = yVec)
  lassoModel = glmnet::glmnet(x=xMat,y=yVec,
                      intercept=TRUE,standardize=TRUE)
  tempCoef = as.numeric(coef(lassoModel,s = escvRes$glmnet$lambda[escvRes$selindex["ESCV"]]))
  escvCoef = tempCoef[-1]

  muVec = cbind(1,xMat) %*% tempCoef
  sigmaError = sqrt(mean((yVec - muVec)^2))

  selectedGenes = which(escvCoef != 0)
  if (length(selectedGenes) > 0) {
    tempX = xMat[,selectedGenes]
    tempGLM = try(glm(yVec ~ tempX,family="gaussian"))
    if (class(tempGLM)[1] != "try-error") {
      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(escvCoef != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
        rownames(tempCI) <- as.character(which(escvCoef != 0))
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])
      if (class(pVals)[1] != "try-error") {
        names(pVals) = as.character(which(escvCoef != 0))
      }
      else {
        pVals = "summary Failed"
      }
    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
    }
  }
  else {
    tempCI = NULL
    pVals = NULL
  }

  resList = list(mu = muVec,
                 sigmaError = sigmaError,
                 coef = escvCoef,
                 activeSet = which(escvCoef != 0),
                 tempCI = tempCI,
                 pVals = pVals,
                 escvRes = escvRes)
  return(resList)
}

#################################################################################

getFSAICCoef_Normal <- function(xMat,yVec,sigmaError) {
  n = dim(xMat)[1]
  p = dim(xMat)[2]

  #### filter on marginal correlations first to reduce the data down
  #### to a max p of 100
  if (p <= 500) {
    numSteps = 120
    fsfit = seletiveInference::fs(xMat,yVec) #,maxsteps=numSteps)
    out.seq = selectiveInference::fsInf(fsfit,alpha=0.05,sigma=sigmaError,
                    type = "aic") #,k=numSteps-1)
    #browser()
    ### get an estimate of error sigma
    betaMat = fsfit$beta
    numSteps = dim(betaMat)[2]
    beta = betaMat[,numSteps]
    numActive = length(out.seq$vars)

    finalCoef = beta
    tempActiveSet = out.seq$vars
    for (i in 1:numActive) {
      tempActiveVar = tempActiveSet[i]
      finalCoef[tempActiveVar] = ifelse(out.seq$pv[i] < 0.05,finalCoef[tempActiveVar],0)
    }
    activeSetBeta = beta
  }
  else {
    corFunc = function(x) {
      res = abs(cor(x,yVec))
      return(res)
    }
    corVec = apply(xMat,2,corFunc)

    cor100th = sort(corVec,decreasing=TRUE)[100]

    largeCorInds = which(corVec >= cor100th)
    xMat_Star = xMat[,largeCorInds]
    fsfit = selectiveInference::fs(xMat_Star,yVec)
    out.seq = selectiveInference::fsInf(fsfit,alpha=0.05,sigma=NULL,
                    type = "active",k=NULL)

    betaMat = fsfit$beta
    numSteps = dim(betaMat)[2]
    beta = betaMat[,numSteps]
    numActive = length(which(beta != 0))

    finalCoef = beta
    tempActiveSet = which(beta != 0)
    for (i in 1:numActive) {
      tempActiveVar = tempActiveSet[i]
      finalCoef[tempActiveVar] = ifelse(out.seq$pv[i] < 0.05,finalCoef[tempActiveVar],0)
    }

    finalCoef2 = rep(0,p)
    for (i in 1:numActive) {
      tempActiveVar = tempActiveSet[i]
      tempOrigInd = largeCorInds[i]
      finalCoef2[tempOrigInd] = finalCoef[tempActiveVar]
    }
    finalCoef = finalCoef2

    activeSetBeta = rep(0,p)
    for (i in 1:numActive) {
      tempActiveVar = tempActiveSet[i]
      tempOrigInd = largeCorInds[i]
      activeSetBeta[tempOrigInd] = beta[tempActiveVar]
    }
  }

  resList = list(coef = finalCoef,
                 activeSet = out.seq$vars,
                 pvalues = out.seq$pv,
                 tempCI = out.seq$ci,
                 fsfit = fsfit)
  return(resList)

}

##################################################################################


getISISCoef_Normal = function(xMat,yVec) {
  tempSeed = as.numeric(sample(1:1000,size=1))
  isisOfficialRes = SIS::SIS(x = xMat,y = yVec,family="gaussian",
                        penalty = "lasso",tune = "cv",seed=tempSeed)
  tempCoef = isisOfficialRes$coef.est
  isisOfficialVars = as.numeric(gsub(x=names(tempCoef),pattern="X",replacement=""))
  isisOfficialVars = isisOfficialVars[which(isisOfficialVars != 0)]

  p = dim(xMat)[2]

  #isisOfficialCoef = isisOfficialCoef[-1]



  selectedGenes = isisOfficialVars
  if (length(selectedGenes) > 0) {
    tempX = xMat[,selectedGenes]
    tempGLM = try(lm(yVec ~ tempX))
    if (class(tempGLM)[1] != "try-error") {
      tempCoef2 = coef(tempGLM)
      muVec = cbind(1,tempX) %*% tempCoef2
      sigmaError = sqrt(mean((yVec - muVec)^2))

      isisOfficialCoef = rep(0,p)
      isisOfficialCoef[isisOfficialVars] = tempCoef2[-1]

      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(isisOfficialCoef != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
        rownames(tempCI) <- as.character(which(isisOfficialCoef != 0))
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])
      if (class(pVals)[1] != "try-error") {
        names(pVals) = as.character(which(isisOfficialCoef != 0))
      }
      else {
        pVals = "summary Failed"
      }
    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
      muVec = "glm failed"
      sigmaError = "glm failed"
      isisOfficialCoef = "glm failed"
    }
  }
  else {
    tempCI = NULL
    pVals = NULL
    muVec = NULL
    sigmaError = NULL
    isisOfficialCoef = NULL
  }

  resList = list(coef = isisOfficialCoef,
                 activeSet = which(isisOfficialCoef != 0),
                 mu = muVec,
                 sigmaError = sigmaError,
                 tempCI = tempCI,
                 pVals = pVals,
                 isisOfficialRes = isisOfficialRes)
  return(resList)
}

#############################################################################

getLassoCoef_Normal = function(xMat,yVec) {
  lassoCV = glmnet::cv.glmnet(x=xMat,y=yVec,
                      family="gaussian",standardize=TRUE,intercept=TRUE)
  lassoCoef = as.numeric(coef(lassoCV))
  muVec = cbind(1,xMat) %*% lassoCoef
  sigmaError = sqrt(mean((muVec - yVec)^2))

  tempX = xMat[,which(lassoCoef[-1] != 0)]
  #tempM = glm(yVec ~ tempX,family="gaussian")

  selectedGenes = which(lassoCoef[-1] != 0)
  if (length(selectedGenes) > 0) {
    tempX = xMat[,selectedGenes]
    tempGLM = try(lm(yVec ~ tempX))
    if (class(tempGLM)[1] != "try-error") {
      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(lassoCoef[-1] != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
        rownames(tempCI) <- as.character(which(lassoCoef[-1] != 0))
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])
      if (class(pVals)[1] != "try-error") {
        names(pVals) = as.character(which(lassoCoef[-1] != 0))
      }
      else {
        pVals = "summary Failed"
      }
    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
    }
  }
  else {
    tempCI = NULL
    pVals = NULL
  }


  resList = list(mu = muVec,
                 sigmaError = sigmaError,
                 coef = lassoCoef[-1],
                 activeSet = which(lassoCoef[-1] != 0),
                 tempCI = tempCI,
                 pVals = pVals,
                 lassoCV = lassoCV)
  return(resList)
}

#################################################################################

getLassoInfCoef_Normal <- function(xMat,yVec,sigma=NULL,lambdaStrategy="Negahban",lambdaVal=NULL) {
  N = dim(xMat)[1]
  p = dim(xMat)[2]


  cvResult = glmnet::cv.glmnet(x=xMat,y=yVec,family="gaussian")
  predValues = predict(cvResult,newx=xMat)
  sigmaHat_1 = sqrt(mean((yVec - predValues)^2))

  if (lambdaStrategy == "InfNorm") {

    # choose lam with simulation
    NSim = 100
    maxValVec = rep(NA,NSim)
    for (i in 1:NSim) {
      tempEps = rnorm(n=N,mean=0,sd=sigmaHat_1)
      tempVec = as.numeric(t(xMat) %*% tempEps)

      maxValVec[i] = max(abs(tempVec))
    }
    lam = 2 * mean(maxValVec)
  }
  if (lambdaStrategy == "Negahban") {
    lam = N*sigmaHat_1*sqrt(2* log(p)/N)
  }
  if (lambdaStrategy == "Custom") {
    lam = lambdaVal
  }
  ##


  ## copied from the documentation
  gfit = glmnet::glmnet(xMat,yVec,standardize=FALSE,thresh=1e-10,lambda.min.ratio = 0.0001)
  #browser()
  # extract coef for a given lambda; note the 1/n factor!
  # (and we don't save the intercept term)
  #lambda = 2
  beta = coef(gfit, x=xMat, y=yVec, s=lam/N, exact=TRUE)[-1]

  #browser()

  if (sum(as.numeric(beta) > 0)) {
    if (is.null(sigma)) {
      predY = as.numeric(xMat %*% beta)
      sigma = sqrt(mean((predY - yVec)^2))
    }

    out = selectiveInference::fixedLassoInf(xMat,yVec,beta,lam,sigma=sigma,alpha=0.05)

    tempCI = out$ci
    pVals = out$pv
    if (is.null(tempCI)) {
      finalCoef = rep(0,p)
      activeSet = NULL

    }
    else {
      activeSet = out$vars

    }
    resList = list(coef = beta,
                   activeSet = activeSet,
                   tempCI = tempCI,
                   pVals = pVals)
  }
  else {
    resList = list(coef = rep(0,p),
                   activeSet = NULL,
                   tempCI = NULL,
                   pVals = NULL)
  }

  return(resList)
}

###################################################################################

getMCPCoef_Normal = function(xMat,yVec) {
  mcpCV = ncvreg::cv.ncvreg(X=xMat,y=yVec,
                    family="gaussian",standardize=TRUE,intercept=TRUE)
  mcpCoef = as.numeric(coef(mcpCV))
  muVec = cbind(1,xMat) %*% mcpCoef
  sigmaError = sqrt(mean((muVec - yVec)^2))

  selectedGenes = which(mcpCoef[-1] != 0)
  if (length(selectedGenes) > 0) {
    tempX = xMat[,selectedGenes]
    tempGLM = try(glm(yVec ~ tempX,family="gaussian"))
    if (class(tempGLM)[1] != "try-error") {
      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(mcpCoef[-1] != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
        rownames(tempCI) <- as.character(which(mcpCoef[-1] != 0))
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])
      if (class(pVals)[1] != "try-error") {
        names(pVals) = as.character(which(mcpCoef[-1] != 0))
      }
      else {
        pVals = "summary Failed"
      }
    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
    }
  }
  else {
    tempCI = NULL
    pVals = NULL
  }


  resList = list(mu = muVec,
                 sigmaError = sigmaError,
                 activeSet = which(mcpCoef[-1] != 0),
                 coef = mcpCoef[-1],
                 tempCI = tempCI,
                 pVals = pVals)
  return(resList)
}

#############################################################################

getSCADCoef_Normal = function(xMat,yVec) {
  scadCV = ncvreg::cv.ncvreg(X=xMat,y=yVec,
                     family="gaussian",standardize=TRUE,intercept=TRUE,
                     penalty = "SCAD")
  scadCoef = as.numeric(coef(scadCV))
  muVec = cbind(1,xMat) %*% scadCoef
  sigmaError = sqrt(mean((muVec - yVec)^2))

  selectedGenes = which(scadCoef[-1] != 0)
  if (length(selectedGenes) > 0) {
    tempX = xMat[,selectedGenes]
    tempGLM = try(glm(yVec ~ tempX,family="gaussian"))
    if (class(tempGLM)[1] != "try-error") {
      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(scadCoef[-1] != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
        rownames(tempCI) <- as.character(which(scadCoef[-1] != 0))
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])
      if (class(pVals)[1] != "try-error") {
        names(pVals) = as.character(which(scadCoef[-1] != 0))
      }
      else {
        pVals = "summary Failed"
      }
    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
    }
  }
  else {
    tempCI = NULL
    pVals = NULL
  }


  resList = list(mu = muVec,
                 sigmaError = sigmaError,
                 activeSet = which(scadCoef[-1] != 0),
                 coef = scadCoef[-1],
                 tempCI = tempCI,
                 pVals = pVals)
  return(resList)
}

####################################################################################

getSOILHighCoef_Normal = function(xMat,yVec) {
  #browser()
  lassoCV = cv.glmnet(x=xMat,y=yVec,
                      family="gaussian",standardize=TRUE,intercept=TRUE)
  lassoCoef = as.numeric(coef(lassoCV))

  activeVars = which(lassoCoef[-1] != 0)
  numActive = length(activeVars)
  p = dim(xMat)[2]
  if (numActive > 0) {
    soilRes = SOIL(x=xMat,y=yVec,psi=1,family="gaussian")

    #soilCutoff = sort(as.numeric(soilRes$importance),decreasing=TRUE)[numActive]
    soilCutoff = 0.5

    whichVars = which(as.numeric(soilRes$importance) >= soilCutoff)
    tempX = xMat[,whichVars]

    tempGLM = try(lm(yVec ~ tempX))
    if (class(tempGLM)[1] != "try-error") {
      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(lassoCoef[-1] != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])

    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
    }
    finalCoef = rep(0,p)
    finalCoef[whichVars] = lassoCoef[whichVars+1]
    activeSet = whichVars
  }
  else {
    finalCoef = lassoCoef
    activeSet = NULL
    tempCI = "none selected"
    pVals = "none selected"
  }

  resList = list(coef = finalCoef,
                 activeSet = activeSet,
                 pvalues = pVals,
                 tempCI = tempCI)
  return(resList)
}

######################################################################################

getSOILMedCoef_Normal = function(xMat,yVec) {
  #browser()
  lassoCV = cv.glmnet(x=xMat,y=yVec,
                      family="gaussian",standardize=TRUE,intercept=TRUE)
  lassoCoef = as.numeric(coef(lassoCV))

  activeVars = which(lassoCoef[-1] != 0)
  numActive = length(activeVars)
  p = dim(xMat)[2]
  if (numActive > 0) {
    soilRes = SOIL(x=xMat,y=yVec,psi=1,family="gaussian")

    #soilCutoff = sort(as.numeric(soilRes$importance),decreasing=TRUE)[numActive]
    soilCutoff = 0.25

    whichVars = which(as.numeric(soilRes$importance) >= soilCutoff)
    tempX = xMat[,whichVars]

    tempGLM = try(lm(yVec ~ tempX))
    if (class(tempGLM)[1] != "try-error") {
      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(lassoCoef[-1] != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
        #rownames(tempCI) <- as.character(whichVars)
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])
      #if (class(pVals)[1] != "try-error") {
      #  names(pVals) = as.character(which(lassoCoef[-1] != 0))
      #}
      #else {
      #  pVals = "summary Failed"
      #}
    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
    }
    finalCoef = rep(0,p)
    finalCoef[whichVars] = lassoCoef[whichVars+1]
    activeSet = whichVars
  }
  else {
    finalCoef = lassoCoef
    activeSet = NULL
    tempCI = "none selected"
    pVals = "none selected"
  }

  resList = list(coef = finalCoef,
                 activeSet = activeSet,
                 pvalues = pVals,
                 tempCI = tempCI)
  return(resList)
}

################################################################################3

getSOILLowCoef_Normal = function(xMat,yVec) {
  #browser()
  lassoCV = cv.glmnet(x=xMat,y=yVec,
                      family="gaussian",standardize=TRUE,intercept=TRUE)
  lassoCoef = as.numeric(coef(lassoCV))

  activeVars = which(lassoCoef[-1] != 0)
  numActive = length(activeVars)
  p = dim(xMat)[2]
  if (numActive > 0) {
    soilRes = SOIL(x=xMat,y=yVec,psi=1,family="gaussian")

    #soilCutoff = sort(as.numeric(soilRes$importance),decreasing=TRUE)[numActive]
    soilCutoff = 0.1

    whichVars = which(as.numeric(soilRes$importance) >= soilCutoff)
    tempX = xMat[,whichVars]

    tempGLM = try(lm(yVec ~ tempX))
    if (class(tempGLM)[1] != "try-error") {
      tempCI = try(confint(tempGLM))
      if (class(tempCI)[1] != "try-error") {
        tempCI = tempCI[-1,]
        if (sum(lassoCoef[-1] != 0) == 1) {
          tempCI = matrix(tempCI,nrow=1,ncol=2)
        }
        #rownames(tempCI) <- as.character(whichVars)
      }
      else {
        tempCI = "confint Failed"
      }

      pVals = try(coef(summary(tempGLM))[-1,4])

    }
    else {
      tempCI = "glm Failed"
      pVals = "glm Failed"
    }
    finalCoef = rep(0,p)
    finalCoef[whichVars] = lassoCoef[whichVars+1]
    activeSet = whichVars
  }
  else {
    finalCoef = lassoCoef
    activeSet = NULL
    tempCI = "none selected"
    pVals = "none selected"
  }

  resList = list(coef = finalCoef,
                 activeSet = activeSet,
                 pvalues = pVals,
                 tempCI = tempCI)
  return(resList)
}
