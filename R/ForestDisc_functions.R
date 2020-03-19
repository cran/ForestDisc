
ForestDisc=function(data,id_target,ntree=50,max_splits=10,opt_meth="NelderMead")
{
  X=data[,-id_target]
  Y=data[,id_target]
  ntree=ntree
  max_splits=max_splits
  opt_meth=opt_meth
  SelectedTREES=RF2Selectedtrees(X=X,Y=Y,ntree=ntree)

  cont_splits=Extract_cont_splits(SelectedTREES)
  cont_variables=cont_splits$continuous_var
  SS=Select_cont_splits(cont_splits=cont_splits,max_splits=max_splits,opt_meth)
  allsplits=SS$All_splits
  splitsselection=SS$Selected_splits

  data_disc=data
  for (i in unique(splitsselection[,1]))
  {
    splitpoints=as.numeric(splitsselection[which(splitsselection[,1]==i),2])
    data_disc[,i]=cut(data_disc[,i],c(-Inf,as.numeric(splitpoints),Inf))

  }

  C=splitsselection
  cutp=vector("list", length(unique(splitsselection$var)))

  V=unique(C$var)

  for (i in 1:length(V))
  {
    if (length(V)==1)
    {
      varm=V
    }else
    {
      varm=V[i]
    }

    match=which(C[,"var"] == varm)
    if (length(match)<1)
    {
      cutp[[i]]=unique(C[match,"Split"])
    }else
    {
      for (j in 1:(length(match)))
      {
        cutp[[i]]=c(cutp[[i]],C[match[j],"Split"])
      }
      cutp[[i]]=unique(cutp[[i]])
    }

  }

  return(list(data_disc=data_disc,cut_points=splitsselection,opt_results=allsplits,cont_variables=cont_variables,listcutp=cutp))
}


RF2Selectedtrees=function(X,Y,ntree,max_TreeRules = 'default',min_RuleSupport = 'default')
{

  if (max_TreeRules != 'default' & min_RuleSupport !='default')
  {
    rf = randomForest::randomForest(x=X,y=Y, ntree=ntree,nodesize=round(min_RuleSupport*nrow(X),0), keep.forest=TRUE,norm.votes=FALSE,maxnodes=max_TreeRules)
  } else if (max_TreeRules != 'default' & min_RuleSupport =='default')
  {
    rf = randomForest::randomForest(x=X,y=Y, ntree=ntree, keep.forest=TRUE,norm.votes=FALSE,maxnodes=max_TreeRules)
  } else if (max_TreeRules == 'default' & min_RuleSupport !='default')
  {
    rf = randomForest::randomForest(x=X,y=Y, ntree=ntree,nodesize=round(min_RuleSupport*nrow(X),0), keep.forest=TRUE,norm.votes=FALSE)
  } else {
    rf = randomForest::randomForest(x=X,y=Y, ntree=ntree, keep.forest=TRUE,norm.votes=FALSE)
  }
  # rf$forest$xlevels

  ID.trees.total=1:ntree
  Selectedtrees = NULL
  Selectedtrees$ntree = length(ID.trees.total)
  Selectedtrees$list = vector("list",Selectedtrees$ntree)
  for(i in 1:Selectedtrees$ntree)
  {
    Selectedtrees$list[[i]] = randomForest::getTree(rf,k=as.integer(ID.trees.total[i]),labelVar=FALSE)
  }
  Selectedtrees$RF=rf
  Selectedtrees$xlevels=rf$forest$xlevels

  continuous_var=vector()
  categorical_var=vector()
  for (i in 1:length(Selectedtrees$xlevels))
  {
    if (is.numeric(Selectedtrees$xlevels[[i]]))
    {
      continuous_var=c(continuous_var,i)
    }else
    {
      categorical_var=c(categorical_var,i)
    }
  }

  allt=as.data.frame(Selectedtrees$list[1])
  if (Selectedtrees$ntree > 1)
  {
    for (i in 2:Selectedtrees$ntree)
    {
      allt=rbind(allt,as.data.frame(Selectedtrees$list[i]))
    }
  }
  for (var in continuous_var)
  {
    Selectedtrees$xlevels[[var]]=sort(c(-Inf,unique(round(allt[which(allt[,3]== var),4],2)),Inf))
  }
  for (var in categorical_var)
  {
    Selectedtrees$xlevels[[var]]= gsub("\\s","",Selectedtrees$xlevels[[var]])
  }
  Selectedtrees$continuous_var=continuous_var
  Selectedtrees$categorical_var=categorical_var
  # -----------------
  return(Selectedtrees)
}


Extract_cont_splits= function(SelectedTREES)
{
  allt=as.data.frame(SelectedTREES$list[1])
  if (SelectedTREES$ntree > 1)
  {
    for (i in 2:SelectedTREES$ntree)
    {
      allt=rbind(allt,as.data.frame(SelectedTREES$list[i]))
    }
  }

  continuous_var=SelectedTREES$continuous_var
  categorical_var=SelectedTREES$categorical_var

  if (length(continuous_var)>0)
  {
   continuous_splits=allt[which(allt$split.var %in% continuous_var),c("split.var", "split.point")]
  }
  return(list(continuous_var=continuous_var,continuous_splits=continuous_splits))
}

Select_cont_splits= function(cont_splits,max_splits,opt_meth)
{
  continuous_var = cont_splits$continuous_var
  continuous_splits= cont_splits$continuous_splits

  if (length(continuous_var)>0)
  {
    All_splits= vector("list",length(continuous_var))
    names(All_splits)= continuous_var
    i=1

    for (var in c(continuous_var))
    {
      All_splits[i]= var
      All_splits[[i]]=continuous_splits[which(continuous_splits$split.var == var),"split.point"]
      i=i+1
    }
    # remove continious variable not used in RF splitting
    removeattrib=which(lengths(All_splits)<1)
    if (length(removeattrib)>0)
    {
    All_splits=All_splits[-removeattrib]
    continuous_var=continuous_var[-removeattrib]
    }

    Selected_splits= vector("list",length(All_splits)) # splits to be found by optimisation pb
    names(Selected_splits)= names(All_splits)
    summary_find_splits=data.frame(matrix(ncol=20 + 5))
    names(summary_find_splits)=c("var","nsplits","optim_value","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","score","uniq_split")

    i=1
    all_Split_dist=data.frame()

    for (j in 1: length(All_splits))
    {

      RF_Splits=All_splits[[j]]
      w=max(max(RF_Splits),-min(RF_Splits),1)
      decimal_count=vector()
      for (k in 1: length(RF_Splits))
      {
        decimal_count[k]=nchar(strsplit(format(RF_Splits[k], scientific = FALSE),".",fixed=TRUE)[[1]][2])
      }
      decimal_count[which(is.na(decimal_count))]=0
      max.decimal=max(decimal_count,na.rm=TRUE)
      S_moments=moments::all.moments(RF_Splits, order.max = 4, central = FALSE, absolute = FALSE, na.rm = FALSE)
      Em=S_moments[2:5]
      S_min=min(RF_Splits)
      S_max=max(RF_Splits)
      S_mean= mean(RF_Splits)
      S_med= stats::median(RF_Splits)
      W=c(1/w^2,1/w^4,1/w^6,1/w^8)
      for (nsplit in 2:max_splits)
      # for (nsplit in 1:max_splits)
      {
        n=nsplit
        x0.splits=c(rep(S_min,n),rep(0,n))
        # x0.splits[1]=S_med
        x0.splits[2]=S_mean
        x0.splits[(n+1):(n+2)]=1/2
        x0.splits[n+1]=1
        fn.splits <- function(x)
        {
          W[1]*(sum(x[1:n]*x[(n+1):(2*n)])-Em[1])^2 + W[2]*(sum(x[1:n]^2*x[(n+1):(2*n)])-Em[2])^2 +W[3]*(sum(x[1:n]^3*x[(n+1):(2*n)])-Em[3])^2 +W[4]*(sum(x[1:n]^4*x[(n+1):(2*n)])-Em[4])^2 + (1-sum(x[(n+1):(2*n)]))^2
        }
        Lower=c(rep(S_min,n),rep(0,n))
        Upper=c(rep(S_max,n),rep(1,n))
        heq.splits <- function(x)
        {
          h = (1-sum(x[(n+1):(2*n)]))
          return(h)
        }

        if (opt_meth=="SLSQP")
        {
          find_best_Splits = nloptr::slsqp(x0.splits, fn = fn.splits,heq=heq.splits,lower=Lower, upper = Upper,
                                   control = list(xtol_rel = 1e-25, check_derivatives = TRUE, maxeval=4000))
        }else if (opt_meth=="NelderMead")
        {
          find_best_Splits = nloptr::neldermead(x0.splits, fn = fn.splits,lower=Lower, upper = Upper,
                                        control = list(xtol_rel = 1e-25, maxeval=4000))
        }else if(opt_meth=="directL") #
        {
          find_best_Splits=nloptr::directL(x0.splits, fn = fn.splits,
                                   lower=Lower, upper = Upper,nl.info = TRUE,
                                   control = list(xtol_rel = 1e-25, maxeval=4000))
        }

        summary_find_splits[i,1]=continuous_var[j]
        summary_find_splits[i,2]=n
        SP=data.frame(find_best_Splits$par[1:n],find_best_Splits$par[(n+1):(2*n)])
        SP=SP[order(SP[,2],decreasing = TRUE),]
        SP[,1]=round(SP[,1],max.decimal)
        SP[,2]=round(SP[,2],4)
        summary_find_splits[i,3]=find_best_Splits$value
        summary_find_splits[i,4:(3+n)] = SP[,1]
        summary_find_splits[i,14:(13+n)] = SP[,2]

        U_P=which(summary_find_splits[i,14:(13+n)]>0) # probability>0
        U_S=as.data.frame(summary_find_splits[i,U_P+3]) # splits correcponding to proba>0
        U_S=U_S[1,which(U_S[1,]>0 | U_S[1,]<=0)] #removing NA
        summary_find_splits[i,"uniq_split"]=length(unique(t(U_S)))
        summary_find_splits[i,"score"]= summary_find_splits[i,"optim_value"]*10^(summary_find_splits[i,"uniq_split"]-2)
        i=i+1

      }

      # Selected_splits (minimize opt value + select unique splits which probability>0)
      A=summary_find_splits[which(summary_find_splits[,1]==continuous_var[j]),]
      A=A[which(A[,"optim_value"]==min(A[,"optim_value"])),] # min optimum value for variable j
      A=A[1,] # to select one line if there is more than one row whith the same opt value
      B=which(A[1,14:(13+n)]>0) # probability>0
      Selected_splits[[j]]=unique(A[1,B+3]) # unique splits
      S=A[1,B+3] # splits
      P=A[1,B+13] # probabilities
      SP=data.frame(t(S),t(P))
      Split_dist=stats::aggregate(list(freq=SP[,2]), by=list(Split=SP[,1]), FUN=sum) # unique splits
      Split_dist[,"var"]=continuous_var[j]
      Split_dist[,"opt_value"]=A$optim_value
      Split_dist=Split_dist[,c(3,1,2,4)]
      Split_dist=Split_dist[order(Split_dist[,2]),]
      all_Split_dist=rbind(all_Split_dist,Split_dist)

    }
  }

  return(list(All_splits=summary_find_splits,Selected_splits=all_Split_dist))
}





