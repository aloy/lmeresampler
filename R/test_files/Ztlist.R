Ztlist.lme <- function(model) {
  
  Z <- Matrix(extract.lmeDesign(model)$Z)
  Zt <- t(Z)
  
  getInds <- function(i) {
    n <- model$dims$ngrps[i] *  model$dims$qvec[i]           ## number of elements in this block
    #     n <- diff(object@Gp)[i]      ## number of elements in this block
    nt <- model$dims$qvec[i]           ## number of REs
    inds <- lapply(seq(nt),seq, to = n, by = nt)  ## pull out individual RE indices
    inds <- lapply(inds, function(x) x )  ## add group offset
  }
  
  inds <- do.call(c,lapply(seq_len(model$dims$Q),getInds))
  
  Ztlist <- lapply(inds,function(i) Zt[i,])
  
  return(Ztlist)
}



## for lme4

tnames <- function(object,diag.only = FALSE,old = TRUE,prefix = NULL) {
  pfun <- mkPfun(diag.only = diag.only, old = old, prefix = prefix)
  c(unlist(mapply(pfun, names(object@cnms), object@cnms)))
}

mkPfun <- function(diag.only = FALSE, old = TRUE, prefix = NULL){
  local({
    function(g,e) {
      mm <- outer(e,e,paste,sep = ".")
      if(old) {
        diag(mm) <- e
      } else {
        mm[] <- paste(mm,g,sep = "|")
        if (!is.null(prefix)) mm[] <- paste(prefix[2],mm,sep = "_")
        diag(mm) <- paste(e,g,sep = "|")
        if (!is.null(prefix))  diag(mm) <- paste(prefix[1],diag(mm),sep = "_")
      }
      mm <- if (diag.only) diag(mm) else mm[lower.tri(mm,diag = TRUE)]
      if(old) paste(g,mm,sep = ".") else mm
    }
  })
}

"Ztlist" =
{
  getInds <- function(i) {
    n2 <- diff(object@Gp)[i]      ## number of elements in this block
    nt2 <- length(cnms[[i]]) ## number of REs
    inds2 <- lapply(seq(nt2),seq,to = n2,by = nt2)  ## pull out individual RE indices
    inds2 <- lapply(inds2,function(x) x + object@Gp[i])  ## add group offset
  }
  inds2 <- do.call(c,lapply(seq_along(cnms),getInds))
  setNames(lapply(inds,function(i) PR$Zt[i,]),
           tnames(object,diag.only = TRUE))
},