#' @title Average mean absolute error (aMAE)
#' @description Calculates average mean absolute error (aMAE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average mean absolute error (aMAE) under multiple, different weighting schemes
#' @details aMAE for weighting scheme # => mean value of the MAEs for specified variables
#' under weighting scheme # => mean value of MAEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVEMAEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVEMAEw
AVEMAEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    avemaevec=vector(length=(length(Weights)+1))
  #calculate unweighted statistic and store value
    avemaeuw=mean(apply(abs(mapply('-', Actual, Survey)), 2, mean))
    avemaevec[1]=avemaeuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avemae=mean(apply(abs(mapply('-', Actual, wtsvy)), 2, mean))
      avemaevec[i+1]=avemae
    }
  #convert to matrix and label
    avemaematrix=noquote(cbind(matrix(format(avemaevec, signif=7))))
    rownames(avemaematrix)=paste("   wgt scheme",
      seq(from=0, along=avemaematrix), " => ")
    rownames(avemaematrix)[1]=paste("   unweighted    => ")
    colnames(avemaematrix)=c("aMAE")
  #return results
    avemaematrix
}

#' @title Average mean squared error (aMSE) with bias-variance decomposition
#' @description Calculates average mean squared error (aMSE) with bias-variance decomposition
#' under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average mean squared error (aMSE) with bias-variance decomposition under multiple, different weighting schemes
#' @details aMSE for weighting scheme # => mean value of the MSEs for specified variables
#' under weighting scheme # => mean value of MSEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVEMSEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVEMSEw
AVEMSEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vectors
    avemsevec=vector(length=(length(Weights)+1))
    avebias2vec=vector(length=(length(Weights)+1))
    avevarvec=vector(length=(length(Weights)+1))
    addvec=vector(length=(length(Weights)+1))
    addvec[1:(length(Weights)+1)]=" + "
    eqvec=vector(length=(length(Weights)+1))
    eqvec[1:(length(Weights)+1)]=" => "
  #calculate unweighted statistics and store values
    avemseuw=mean(apply(((mapply('-', Actual, Survey))^2), 2, mean))
    avemsevec[1]=avemseuw
    avebias2uw=mean((apply((mapply('-', Actual, Survey)), 2, mean))^2)
    avebias2vec[1]=avebias2uw
    avevaruw=mean(apply(((mapply('-', (Actual-Survey),
      (apply((Actual-Survey), 2, mean))))^2), 2, mean))
    avevarvec[1]=avevaruw
  #weight data, calculate statistics, store values
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avemse=mean(apply(((mapply('-', Actual, wtsvy))^2), 2, mean))
      avemsevec[[i+1]]=avemse
      avebias2=mean((apply((mapply('-', Actual, wtsvy)), 2, mean))^2)
      avebias2vec[[i+1]]=avebias2
      avevar=mean(apply(((mapply('-', (Actual-wtsvy),
        (apply((Actual-wtsvy), 2, mean))))^2), 2, mean))
      avevarvec[[i+1]]=avevar
    }
  #convert to matrix and label
    avemsematrix=noquote(matrix(cbind(format(avemsevec, signif=7),
      eqvec, format(avebias2vec, signif=7), addvec, format(avevarvec,
      signif=7)), ncol=5))
    rownames(avemsematrix)=paste("   wgt scheme",
      seq(from=0, along=avemsematrix[, 1]), " => ")
    rownames(avemsematrix)[1]=paste("   unweighted    => ")
    colnames(avemsematrix)=c("aMSE", " ", "aBias^2", " ", "aVar")
  #return results
    avemsematrix
}

#' @title Average root mean squared error (aRMSE)
#' @description Calculates average root mean squared error (aRMSE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average root mean squared error (aRMSE) under multiple, different weighting schemes
#' @details aRMSE for weighting scheme # => mean value of the RMSEs for specified variables
#' under weighting scheme # => mean value of RMSEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVERMSEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVERMSEw
AVERMSEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    avermsevec=vector(length=length(Weights))
  #calculate unweighted statistic and store value
    avermseuw=mean(sqrt(apply(((mapply('-', Actual, Survey))^2), 2, mean)))
    avermsevec[1]=avermseuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avermse=mean(sqrt(apply(((mapply('-', Actual, wtsvy))^2), 2, mean)))
      avermsevec[[i+1]]=avermse
    }
  #convert to matrix and label
    avermsematrix=noquote(cbind(matrix(format(avermsevec, signif=7))))
    rownames(avermsematrix)=paste("   wgt scheme",
      seq(from=0, along=avermsematrix), " => ")
    rownames(avermsematrix)[1]=paste("   unweighted    => ")
    colnames(avermsematrix)=c("aRMSE")
  #return results
    avermsematrix
}

#' @title Average mean squared logarithmic error (aMSLE)
#' @description Calculates average mean squared logarithmic error (aMSLE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average mean squared logarithmic error (aMSLE) under multiple, different weighting schemes
#' @details aMSLE for weighting scheme # => mean value of the aMSLEs for specified variables
#' under weighting scheme # => mean value of aMSLEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVEMSLEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVEMSLEw
AVEMSLEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    avemslevec=vector(length=length(Weights))
  #calculate unweighted statistic and store value
    avemsleuw=mean(apply(((mapply('-', log(Actual+1),
      log(Survey+1)))^2), 2, mean))
    avemslevec[1]=avemsleuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avemsle=mean(apply(((mapply('-', log(Actual+1),
        log(wtsvy+1)))^2), 2, mean))
      avemslevec[[i+1]]=avemsle
    }
  #convert to matrix and label
    avemslematrix=noquote(cbind(matrix(format(avemslevec, signif=7))))
    rownames(avemslematrix)=paste("   wgt scheme",
      seq(from=0, along=avemslematrix), " => ")
    rownames(avemslematrix)[1]=paste("   unweighted    => ")
    colnames(avemslematrix)=c("aMSLE")
  #return results
    avemslematrix
}

#' @title Average root mean squared logarithmic error (aRMSLE)
#' @description Calculates average root mean squared logarithmic error (aRMSLE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average root mean squared logarithmic error (aRMSLE) under multiple, different weighting schemes
#' @details aRMSLE for weighting scheme # => mean value of the aRMSLEs for specified variables
#' under weighting scheme # => mean value of aRMSLEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVERMSLEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVERMSLEw
AVERMSLEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    avermslevec=vector(length=length(Weights))
  #calculate unweighted statistic and store value
    avermsleuw=mean(sqrt(apply(((mapply('-', log(Actual+1),
      log(Survey+1)))^2), 2, mean)))
    avermslevec[1]=avermsleuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avermsle=mean(sqrt(apply(((mapply('-', log(Actual+1),
        log(wtsvy+1)))^2), 2, mean)))
      avermslevec[[i+1]]=avermsle
    }
  #convert to matrix and label
    avermslematrix=noquote(cbind(matrix(format(avermslevec, signif=7))))
    rownames(avermslematrix)=paste("   wgt scheme",
      seq(from=0, along=avermslematrix), " => ")
    rownames(avermslematrix)[1]=paste("   unweighted    => ")
    colnames(avermslematrix)=c("aRMSLE")
  #return results
    avermslematrix
}

#' @title Full scale-dependent statistics
#' @description Calculates full scale-dependent statistics
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Full scale-dependent statistics
#' @examples FULLSDw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export FULLSDw
FULLSDw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vectors
    avemaevec=vector(length=(length(Weights)+1))
    avemsevec=vector(length=(length(Weights)+1))
    avermsevec=vector(length=(length(Weights)+1))
    avemslevec=vector(length=(length(Weights)+1))
    avermslevec=vector(length=(length(Weights)+1))
    spacevec=vector(length=(length(Weights)+1))
    spacevec[1:(length(Weights)+1)]=" "
  #calculate unweighted statistics and store values
    avemaeuw=mean(apply(abs(mapply('-', Actual, Survey)), 2, mean))
    avemaevec[1]=avemaeuw
    avemseuw=mean(apply(((mapply('-', Actual, Survey))^2), 2, mean))
    avemsevec[1]=avemseuw
    avermseuw=mean(sqrt(apply(((mapply('-', Actual, Survey))^2), 2, mean)))
    avermsevec[1]=avermseuw
    avemsleuw=mean(apply(((mapply('-', log(Actual+1),
      log(Survey+1)))^2), 2, mean))
    avemslevec[1]=avemsleuw
    avermsleuw=mean(sqrt(apply(((mapply('-', log(Actual+1),
      log(Survey+1)))^2), 2, mean)))
    avermslevec[1]=avermsleuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avemae=mean(apply(abs(mapply('-', Actual, wtsvy)), 2, mean))
      avemaevec[[i+1]]=avemae
      avemse=mean(apply(((mapply('-', Actual, wtsvy))^2), 2, mean))
      avemsevec[[i+1]]=avemse
      avermse=mean(sqrt(apply(((mapply('-', Actual, wtsvy))^2), 2, mean)))
      avermsevec[[i+1]]=avermse
      avemsle=mean(apply(((mapply('-', log(Actual+1),
        log(wtsvy+1)))^2), 2, mean))
      avemslevec[[i+1]]=avemsle
      avermsle=mean(sqrt(apply(((mapply('-', log(Actual+1),
        log(wtsvy+1)))^2), 2, mean)))
      avermslevec[[i+1]]=avermsle
    }
  #convert to matrix and label
    fullsdmatrix=noquote(matrix(cbind(format(avemaevec, signif=7), spacevec,
      format(avemsevec, signif=7), spacevec, format(avermsevec, signif=7),
      spacevec, format(avemslevec, signif=7), spacevec,
      format(avermslevec, signif=7)), ncol=9))
    rownames(fullsdmatrix)=paste("   wgt scheme",
      seq(from=0, along=fullsdmatrix[, 1]), " => ")
    rownames(fullsdmatrix)[1]=paste("   unweighted    => ")
    colnames(fullsdmatrix)=c("aMAE", " ", "aMSE", " ", "aRMSE", " ",
      "aMSLE", " ", "aRMSLE")
  #return results
    fullsdmatrix
}

#' @title Average mean absolute percentage error (aMAPE)
#' @description Calculates average mean absolute percentage error (aMAPE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average mean absolute percentage error (aMAPE) under multiple, different weighting schemes
#' @details aMAPE for weighting scheme # => mean value of the aMAPEs for specified variables
#' under weighting scheme # => mean value of aMAPEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVEMAPEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVEMAPEw
AVEMAPEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    avemapevec=vector(length=(length(Weights)+1))
  #calculate unweighted statistic and store value
    avemapeuw=mean(apply(abs(mapply('-', Actual, Survey))/Actual, 2, mean))
    avemapevec[1]=avemapeuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avemape=mean(apply(abs(mapply('-', Actual, wtsvy))/Actual, 2, mean))
      avemapevec[[i+1]]=avemape
    }
  #convert to matrix and label
    avemapematrix=noquote(cbind(matrix(format(avemapevec, signif=7))))
    rownames(avemapematrix)=paste("   wgt scheme",
      seq(from=0, along=avemapematrix), " => ")
    rownames(avemapematrix)[1]=paste("   unweighted    => ")
    colnames(avemapematrix)=c("aMAPE")
  #return results
    avemapematrix
}

#' @title Average symmetric mean absolute percentage error (aSMAPE)
#' @description Calculates average symmetric mean absolute percentage error (aSMAPE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average symmetric mean absolute percentage error (aSMAPE) under multiple, different weighting schemes
#' @details aSMAPE for weighting scheme # => mean value of the aSMAPEs for specified variables
#' under weighting scheme # => mean value of aSMAPEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVESMAPEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVESMAPEw
AVESMAPEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    avesmapevec=vector(length=(length(Weights)+1))
  #calculate unweighted statistic and store value
    avesmapeuw=mean((apply((abs(mapply('-', Actual, Survey))/
      (abs(Actual)+abs(Survey))), 2, mean))*2)
    avesmapevec[1]=avesmapeuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avesmape=mean((apply((abs(mapply('-', Actual, wtsvy))/
        (abs(Actual)+abs(wtsvy))), 2, mean))*2)
      avesmapevec[[i+1]]=avesmape
    }
  #convert to matrix and label
    avesmapematrix=noquote(cbind(matrix(format(avesmapevec, signif=7))))
    rownames(avesmapematrix)=paste("   wgt scheme",
      seq(from=0, along=avesmapematrix), " => ")
    rownames(avesmapematrix)[1]=paste("   unweighted    => ")
    colnames(avesmapematrix)=c("aSMAPE")
  #return results
    avesmapematrix
}

#' @title Average relative absolute error (aRAE)
#' @description Calculates average relative absolute error (aRAE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average relative absolute error (aRAE) under multiple, different weighting schemes
#' @details aRAE for weighting scheme # => mean value of the aRAEs for specified variables
#' under weighting scheme # => mean value of aRAEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVERAEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVERAEw
AVERAEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    averaevec=vector(length=(length(Weights)+1))
  #calculate unweighted statistic and store value
    averaeuw=mean((apply(abs(mapply('-', Actual, Survey)), 2, sum))/
      apply((abs(mapply('-', Actual, apply(Actual, 2, mean)))), 2, sum))
    averaevec[1]=averaeuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      averae=mean((apply(abs(mapply('-', Actual, wtsvy)), 2, sum))/
        apply((abs(mapply('-', Actual, apply(Actual, 2, mean)))), 2, sum))
      averaevec[[i+1]]=averae
    }
  #convert to matrix and label
    averaematrix=noquote(cbind(matrix(format(averaevec, signif=7))))
    rownames(averaematrix)=paste("   wgt scheme",
      seq(from=0, along=averaematrix), " => ")
    rownames(averaematrix)[1]=paste("   unweighted    => ")
    colnames(averaematrix)=c("aRAE")
  #return results
    averaematrix
}

#' @title Average relative squared error (aRSE)
#' @description Calculates average relative squared error (aRSE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average relative squared error (aRSE) under multiple, different weighting schemes
#' @details aRSE for weighting scheme # => mean value of the aRSEs for specified variables
#' under weighting scheme # => mean value of aRSEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVERSEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVERSEw
AVERSEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    aversevec=vector(length=(length(Weights)+1))
  #calculate unweighted statistic and store value
    averseuw=mean((apply((mapply('-', Actual, Survey))^2, 2, sum))/
      apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum))
    aversevec[1]=averseuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      averse=mean((apply((mapply('-', Actual, wtsvy))^2, 2, sum))/
        apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum))
      aversevec[[i+1]]=averse
    }
  #convert to matrix and label
    aversematrix=noquote(cbind(matrix(format(aversevec, signif=7))))
    rownames(aversematrix)=paste("   wgt scheme",
      seq(from=0, along=aversematrix), " => ")
    rownames(aversematrix)[1]=paste("   unweighted    => ")
    colnames(aversematrix)=c("aRSE")
  #return results
    aversematrix
}

#' @title Average root relative squared error (aRRSE)
#' @description Calculates average root relative squared error (aRRSE) under multiple, different weighting schemes
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Average root relative squared error (aRRSE) under multiple, different weighting schemes
#' @details aRRSE for weighting scheme # => mean value of the aRRSEs for specified variables
#' under weighting scheme # => mean value of aRRSEs for objects in Survey=data.frame() *
#' objects in Weights=data.frame()
#' @examples AVERRSEw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export AVERRSEw
AVERRSEw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vector
    averrsevec=vector(length=(length(Weights)+1))
  #calculate unweighted statistic and store value
    averrseuw=mean(sqrt((apply((mapply('-', Actual, Survey))^2, 2, sum))/
      apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum)))
    averrsevec[1]=averrseuw
  #weight data, calculate statistic, store value
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      averrse=mean(sqrt((apply((mapply('-', Actual, wtsvy))^2, 2, sum))/
        apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum)))
      averrsevec[[i+1]]=averrse
    }
  #convert to matrix and label
    averrsematrix=noquote(cbind(matrix(format(averrsevec, signif=7))))
    rownames(averrsematrix)=paste("   wgt scheme",
      seq(from=0, along=averrsematrix), " => ")
    rownames(averrsematrix)[1]=paste("   unweighted    => ")
    colnames(averrsematrix)=c("aRRSE")
  #return results
    averrsematrix
}

#' @title Full scale-independent statistics
#' @description Calculates full scale-independent statistics
#' @param Actual data from a "gold standard" survey; objects are variable columns
#' from "gold standard" survey that corruspond to variable columns Survey
#' @param Survey data from a survey; objects are variable columns from a
#' survey that corruspond to variable columns from Actual
#' @param Weights weights to be applied to Survey data; objects are weights columns
#' @return Full scale-independent statistics
#' @examples FULLSIw(Actual=data.frame(TESTWGT$A1, TESTWGT$A2),
#' Survey=data.frame(TESTWGT$Q1, TESTWGT$Q2),
#' Weights=data.frame(TESTWGT$W1, TESTWGT$W2))
#' @note Make sure to properly order inputs, per the example: Actual=data.frame() objects
#' and corrusponding Survey=data.frame() objects must be given in the same order as each other;
#' and Weights=data.frame() objects must be given in sequence of weighting scheme #.
#' @export FULLSIw
FULLSIw=function(Actual=data.frame(), Survey=data.frame(), Weights=data.frame()){
  #create storage vectors
    avemapevec=vector(length=(length(Weights)+1))
    avesmapevec=vector(length=(length(Weights)+1))
    averaevec=vector(length=(length(Weights)+1))
    aversevec=vector(length=(length(Weights)+1))
    averrsevec=vector(length=(length(Weights)+1))
    spacevec=vector(length=(length(Weights)+1))
    spacevec[1:(length(Weights)+1)]=" "
  #calculate unweighted statistics and store values
    avemapeuw=mean(apply(abs(mapply('-', Actual, Survey))/Actual, 2, mean))
    avemapevec[1]=avemapeuw
    avesmapeuw=mean((apply((abs(mapply('-', Actual, Survey))/
      (abs(Actual)+abs(Survey))), 2, mean))*2)
    avesmapevec[1]=avesmapeuw
    averaeuw=mean((apply(abs(mapply('-', Actual, Survey)), 2, sum))/
      apply((abs(mapply('-', Actual, apply(Actual, 2, mean)))), 2, sum))
    averaevec[1]=averaeuw
    averseuw=mean((apply((mapply('-', Actual, Survey))^2, 2, sum))/
      apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum))
    aversevec[1]=averseuw
    averrseuw=mean(sqrt((apply((mapply('-', Actual, Survey))^2, 2, sum))/
      apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum)))
    averrsevec[1]=averrseuw
  #weight data, calculate statistics, store values
    for(i in 1:length(Weights)){
      wtsvy=data.frame(lapply(Survey, '*', Weights[[i]]))
      avemape=mean(apply(abs(mapply('-', Actual, wtsvy))/Actual, 2, mean))
      avemapevec[[i+1]]=avemape
      avesmape=mean((apply((abs(mapply('-', Actual, wtsvy))/
        (abs(Actual)+abs(wtsvy))), 2, mean))*2)
      avesmapevec[[i+1]]=avesmape
      averae=mean((apply(abs(mapply('-', Actual, wtsvy)), 2, sum))/
        apply((abs(mapply('-', Actual, apply(Actual, 2, mean)))), 2, sum))
      averaevec[[i+1]]=averae
      averse=mean((apply((mapply('-', Actual, wtsvy))^2, 2, sum))/
        apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum))
      aversevec[[i+1]]=averse
      averrse=mean(sqrt((apply((mapply('-', Actual, wtsvy))^2, 2, sum))/
        apply(((mapply('-', Actual, apply(Actual, 2, mean)))^2), 2, sum)))
      averrsevec[[i+1]]=averrse
    }
  #convert to matrix and label
    fullsimatrix=noquote(matrix(cbind(format(avemapevec, signif=7), spacevec,
      format(avesmapevec, signif=7), spacevec, format(averaevec, signif=7),
      spacevec, format(aversevec, signif=7), spacevec, format(averrsevec,
      signif=7)), ncol=9))
    rownames(fullsimatrix)=paste("   wgt scheme",
      seq(from=0, along=fullsimatrix[, 1]), " => ")
    rownames(fullsimatrix)[1]=paste("   unweighted    => ")
    colnames(fullsimatrix)=c("aMAPE", " ", "aSMAPE", " ", "aRAE", " ",
      "aRSE", " ", "aRRSE")
    #return results
      fullsimatrix
}
