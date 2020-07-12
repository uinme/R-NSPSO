crowding_distance_sort <- function (objValues, opt) {
  # この関数は，各個体の目的関数値が格納された行列objValuesに混雑距離ソートを適用
  # し，その結果をベクトルで返す．
  # 
  # 引数：
  #   objValuesは，次のような行列を想定している．
  # 
  #               |f11, f21, ..., fk1|
  #   objValues = |f12, f22, ..., fk2|
  #               | : ,  : , ...,  : |
  #               |f1p, f2p, ..., fkp|
  #   
  #   fkpは，p番目の個体におけるk番目の目的関数の値．
  # 
  #   optは，最小化問題の場合は"min"に，最大化問題の場合は"max"に設定する．
  # 
  # 戻り値： 
  #   混雑距離が大きい順に並べ替えられたobjValuesの行番号をベクトルで返す．
  #

    numObj <- ncol(objValues)
    numPop <- nrow(objValues)

    popMin <- apply(objValues, 2, min)
    popMax <- apply(objValues, 2, max)

    if (opt == "min" || opt == 0) {
        ixdes <- apply(objValues, 2, order, decreasing=TRUE)
    } else if (opt == "max" || opt == 1) {
        ixdes <- apply(objValues, 2, order, decreasing=FALSE)
    }

    cDist <- matrix(Inf, nrow=numPop, ncol=numObj)
    for (m in 1:numObj) {
        ix <- 2:(numPop-1)
        iplus  <- ixdes[ix+1, m]
        iminus <- ixdes[ix-1, m]
        cDist[ixdes[ix, m], m] <- (objValues[iminus, m] - objValues[iplus, m]) / (popMax[m] - popMin[m])
    }
    cDist <- rowSums(cDist)
    ixcd <- order(cDist, decreasing=TRUE)

    return(ixcd)
}
