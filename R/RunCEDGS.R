
#' @export
RunCEDGS_Normal = function(vsMethods = c("getESCVCoef_Normal",
                                         "getFSAICCoef_Normal",
                                         "getISISCoef_Normal",
                                         "getLassoCoef_Normal",
                                         "getLassoInfCoef_Normal"),
                           NReps = 10,
                           x,y) {



  ### this proceeds in a "round-robin" style.
  ### each method is fit to the data, and the estimated fit is the "home scenario"
  ### new data is then simulated from this fit
  ### then, all of the vs-methods are re-fit on the simulated data
  ### VS-diagnostics are calculated against the simulation ground truth

  # x is the matrix of predictor variables (n by p)
  # y is the continuous reponse vector (n by 1)

  # the methods in vsMethods must perform variable selection
  # and follow the rules for "X_Normal" methods
  # this means returning the set of selected variables, CI's and p-values

  # maybe add a "varsToEmph" argument later for specific results on certain variables

  ####
  RepVec = c()
  HomeVec = c()
  AwayVec = c()

  FStatVec = c()
  GStatVec = c()
  VSDVec = c()
  VSDPlusVec = c()
  VSDMinusVec = c()


  FullCoverageStarVec = c()
  FullCoverageVec = c()
  PartialCoverageVec = c()
  FullCoverageActiveVec = c()
  FullCoverageInactiveVec = c()
  PartialCoverageActiveVec = c()
  PartialCoverageInactiveVec = c()

  MinCILengthVec = c()
  MedCILengthVec = c()


  ### extract basic info about the data
  N = dim(x)[1]
  p = dim(x)[2]


  for (i in 1:NReps) {
    for (homeMethod in vsMethods) {
      homeMethodFunc = get(homeMethod)

      homeMethodFit = homeMethodFunc(x,y)
      tempBeta = homeMethodFit$coef
      tempMeanFunc = as.numeric(x %*% tempBeta)
      tempSigma = sqrt(mean((y - tempMeanFunc)^2,na.rm=TRUE))

      homeActiveSet = homeMethodFit$activeSet
      numHomeActive = length(homeActiveSet)

      ### generate simulated data
      simY = rep(0,N)
      for (j in 1:N) {
        simY[j] = rnorm(n=1,mean=tempMeanFunc[j],sd=tempSigma)
      }

      # use the original X matrix

      for (awayMethod in vsMethods) {
        print(sprintf("%s - %s",homeMethod,awayMethod))
        awayMethodFunc = get(awayMethod)

        tempAwayFit = awayMethodFunc(x,y)

        awayActiveSet = tempAwayFit$activeSet


        RepVec = c(RepVec,i)
        HomeVec = c(HomeVec,homeMethod)
        AwayVec = c(AwayVec,awayMethod)

        numAwayActive = length(awayActiveSet)

        if (numHomeActive > 0) {
          FStatVec = c(FStatVec,FFunc(homeActiveSet,awayActiveSet))
          GStatVec = c(GStatVec,GFunc(homeActiveSet,awayActiveSet))
          VSDMinusVec = c(VSDMinusVec,VSDMinusFunc(homeActiveSet,awayActiveSet))
          VSDPlusVec = c(VSDPlusVec,VSDPlusFunc(homeActiveSet,awayActiveSet))
          VSDVec = c(VSDVec,VSDFunc(homeActiveSet,awayActiveSet))

          if (numAwayActive > 0) {
            tempCI = tempAwayFit$tempCI

            ### Full Target Coverage
            numCovered = 0
            numCoveredActive = 0
            numCoveredInactive = 0
            numTotalActive = 0
            numTotalInactive = 0
            for (j in 1:numAwayActive) {
              tempVar = awayActiveSet[j]
              numCovered = ifelse((tempCI[j,1] <= tempBeta[tempVar]) & (tempCI[j,2] >= tempBeta[tempVar]),
                                  numCovered + 1,numCovered)

              numTotalActive = ifelse(tempBeta[tempVar] != 0,numTotalActive +1,
                                      numTotalActive)
              if (tempBeta[tempVar] != 0) {
                numCoveredActive = ifelse((tempCI[j,1] <= tempBeta[tempVar]) & (tempCI[j,2] >= tempBeta[tempVar]),
                                          numCoveredActive +1,numCoveredActive)
              }
              numTotalInactive = ifelse(tempBeta[tempVar] == 0,numTotalInactive+1,
                                        numTotalInactive)
              if (tempBeta[tempVar] == 0) {
                numCoveredInactive = ifelse((tempCI[j,1] <= tempBeta[tempVar]) & (tempCI[j,2] >= tempBeta[tempVar]),
                                            numCoveredInactive + 1,numCoveredInactive)
              }

            }
            FullCoverageStarVec = c(FullCoverageStarVec,
                                    numCovered/(numAwayActive + length(setdiff(homeActiveSet,awayActiveSet))))
            FullCoverageVec = c(FullCoverageVec,numCovered/numAwayActive)
            FullCoverageActiveVec = c(FullCoverageActiveVec,numCoveredActive/numTotalActive)
            FullCoverageInactiveVec = c(FullCoverageInactiveVec,numCoveredInactive/numTotalInactive)


            ### Partial Target Coverage
            tempX = x[,awayActiveSet]
            tempPartialTarget = as.numeric(solve(t(tempX) %*% tempX) %*% tempX %*% tempMeanFunc)
            numCovered = 0
            numCoveredActive = 0
            numCoveredInactive = 0
            numTotalActive = 0
            numTotalInactive = 0
            for (j in 1:numAwayActive) {
              tempVar = awayActiveSet[j]
              numCovered = ifelse((tempCI[j,1] <= tempPartialTarget[tempVar]) & (tempCI[j,2] >= tempPartialTarget[tempVar]),
                                  numCovered + 1,numCovered)

              numTotalActive = ifelse(tempBeta[tempVar] != 0,numTotalActive +1,
                                      numTotalActive)
              if (tempBeta[tempVar] != 0) {
                numCoveredActive = ifelse((tempCI[j,1] <= tempPartialTarget[tempVar]) & (tempCI[j,2] >= tempPartialTarget[tempVar]),
                                          numCoveredActive +1,numCoveredActive)
              }
              numTotalInactive = ifelse(tempBeta[tempVar] == 0,numTotalInactive+1,
                                        numTotalInactive)
              if (tempBeta[tempVar] == 0) {
                numCoveredInactive = ifelse((tempCI[j,1] <= tempPartialTarget[tempVar]) & (tempCI[j,2] >= tempPartialTarget[tempVar]),
                                            numCoveredInactive + 1,numCoveredInactive)
              }
            }
            PartialCoverageVec = c(PartialCoverageVec,numCovered/numAwayActive)
            PartialCoverageActiveVec = c(PartialCoverageActiveVec,numCoveredActive/numTotalActive)
            PartialCoverageInactiveVec = c(PartialCoverageInactiveVec,numCoveredInactive/numTotalInactive)


            ciLengthVec = as.numeric(tempCI[,2] - tempCI[,1])
            MinCILengthVec = c(MinCILengthVec,min(ciLengthVec,na.rm=TRUE))
            MedCILengthVec = c(MedCILengthVec,median(ciLengthVec,na.rm=TRUE))

          }
          else {
            FullCoverageStarVec = c(FullCoverageStarVec,0)
            FullCoverageVec = c(FullCoverageVec,NA)
            FullCoverageActiveVec = c(FullCoverageActiveVec,NA)
            FullCoverageInactiveVec = c(FullCoverageInactiveVec,NA)

            PartialCoverageVec = c(PartialCoverageVec,NA)
            PartialCoverageActiveVec = c(PartialCoverageActiveVec,NA)
            PartialCoverageInactiveVec = c(PartialCoverageInactiveVec,NA)

            MinCILengthVec = c(MinCILengthVec,NA)
            MedCILengthVec = c(MedCILengthVec,NA)

          }

        }
        else {
          FStatVec = c(FStatVec,NA)
          GStatVec = c(GStatVec,NA)
          VSDMinusVec = c(VSDMinusVec,VSDMinusFunc(homeActiveSet,awayActiveSet))
          VSDPlusVec = c(VSDPlusVec,NA)
          VSDVec = c(VSDVec,NA)

          if (numAwayActive > 0) {
            tempCI = tempAwayFit$tempCI

            ### Full Target Coverage
            numCovered = 0
            numCoveredActive = 0
            numCoveredInactive = 0
            numTotalActive = 0
            numTotalInactive = 0
            for (j in 1:numAwayActive) {
              tempVar = awayActiveSet[j]
              numCovered = ifelse((tempCI[j,1] <= tempBeta[tempVar]) & (tempCI[j,2] >= tempBeta[tempVar]),
                                  numCovered + 1,numCovered)

              numTotalActive = ifelse(tempBeta[tempVar] != 0,numTotalActive +1,
                                      numTotalActive)
              if (tempBeta[tempVar] != 0) {
                numCoveredActive = ifelse((tempCI[j,1] <= tempBeta[tempVar]) & (tempCI[j,2] >= tempBeta[tempVar]),
                                          numCoveredActive +1,numCoveredActive)
              }
              numTotalInactive = ifelse(tempBeta[tempVar] == 0,numTotalInactive+1,
                                        numTotalInactive)
              if (tempBeta[tempVar] == 0) {
                numCoveredInactive = ifelse((tempCI[j,1] <= tempBeta[tempVar]) & (tempCI[j,2] >= tempBeta[tempVar]),
                                            numCoveredInactive + 1,numCoveredInactive)
              }

            }
            FullCoverageStarVec = c(FullCoverageStarVec,
                                    numCovered/(numAwayActive + length(setdiff(homeActiveSet,awayActiveSet))))
            FullCoverageVec = c(FullCoverageVec,numCovered/numAwayActive)
            FullCoverageActiveVec = c(FullCoverageActiveVec,numCoveredActive/numTotalActive)
            FullCoverageInactiveVec = c(FullCoverageInactiveVec,numCoveredInactive/numTotalInactive)


            ### Partial Target Coverage
            tempX = x[,awayActiveSet]
            tempPartialTarget = as.numeric(solve(t(tempX) %*% tempX) %*% tempX %*% tempMeanFunc)
            numCovered = 0
            numCoveredActive = 0
            numCoveredInactive = 0
            numTotalActive = 0
            numTotalInactive = 0
            for (j in 1:numAwayActive) {
              tempVar = awayActiveSet[j]
              numCovered = ifelse((tempCI[j,1] <= tempPartialTarget[tempVar]) & (tempCI[j,2] >= tempPartialTarget[tempVar]),
                                  numCovered + 1,numCovered)

              numTotalActive = ifelse(tempBeta[tempVar] != 0,numTotalActive +1,
                                      numTotalActive)
              if (tempBeta[tempVar] != 0) {
                numCoveredActive = ifelse((tempCI[j,1] <= tempPartialTarget[tempVar]) & (tempCI[j,2] >= tempPartialTarget[tempVar]),
                                          numCoveredActive +1,numCoveredActive)
              }
              numTotalInactive = ifelse(tempBeta[tempVar] == 0,numTotalInactive+1,
                                        numTotalInactive)
              if (tempBeta[tempVar] == 0) {
                numCoveredInactive = ifelse((tempCI[j,1] <= tempPartialTarget[tempVar]) & (tempCI[j,2] >= tempPartialTarget[tempVar]),
                                            numCoveredInactive + 1,numCoveredInactive)
              }
            }
            PartialCoverageVec = c(PartialCoverageVec,numCovered/numAwayActive)
            PartialCoverageActiveVec = c(PartialCoverageActiveVec,numCoveredActive/numTotalActive)
            PartialCoverageInactiveVec = c(PartialCoverageInactiveVec,numCoveredInactive/numTotalInactive)

            ciLengthVec = as.numeric(tempCI[,2] - tempCI[,1])
            MinCILengthVec = c(MinCILengthVec,min(ciLengthVec,na.rm=TRUE))
            MedCILengthVec = c(MedCILengthVec,median(ciLengthVec,na.rm=TRUE))


          }
          else {

            FullCoverageStarVec = c(FullCoverageStarVec,NA)
            FullCoverageVec = c(FullCoverageVec,NA)
            FullCoverageActiveVec = c(FullCoverageActiveVec,NA)
            FullCoverageInactiveVec = c(FullCoverageInactiveVec,NA)

            PartialCoverageVec = c(PartialCoverageVec,NA)
            PartialCoverageActiveVec = c(PartialCoverageActiveVec,NA)
            PartialCoverageInactiveVec = c(PartialCoverageInactiveVec,NA)

            MinCILengthVec = c(MinCILengthVec,NA)
            MedCILengthVec = c(MedCILengthVec,NA)

          }

        }




      }

    }
  }


  resDF = data.frame(Rep = RepVec,
                     HomeMethod = HomeVec,
                     AwayMethod = AwayVec,
                     FStat = FStatVec,
                     GStat = GStatVec,
                     VSDPlus = VSDPlusVec,
                     VSDMinus = VSDMinusVec,
                     VSD = VSDVec,
                     FullCoverageStar = FullCoverageStarVec,
                     FullCoverage = FullCoverageVec,
                     FullCoverageActive = FullCoverageActiveVec,
                     FullCoverageInactive = FullCoverageInactiveVec,
                     PartialCoverage = PartialCoverageVec,
                     PartialCoverageActive = PartialCoverageActiveVec,
                     PartialCoverageInactive = PartialCoverageInactiveVec,
                     MinCILength = MinCILengthVec,
                     MedCILength = MedCILengthVec)
  return(resDF)
}
