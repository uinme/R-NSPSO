non_dominated <- function (objValues, opt) {
  # この関数は，目的関数の値が格納されている行列objValuesに
  # おいて，非劣解となる行番号が格納されているベクトルixndomを返す．
  # 
  # 引数:
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
  #   objValuesのなかで，非劣解となる個体の行番号ベクトルを返す．

  numObj <- ncol(objValues)  # 目的関数の数
  numPop <- nrow(objValues)  # 集団の数

  ixa <- 1:numPop            # 比較対象Aのindex
  ixndom <- NULL             # 非劣解のindexを保存するベクトル

  # 2目的最小化問題の場合における，非劣解の抽出方法 
  # [個体1]
  #    行列a          行列b                      行列                   ベクトル
  # |f11, f21|     |f12, f22|   各成分の比較    |F, F|   各行の論理積     |F|
  # |f11, f21|  >  |f13, f23| --結果が行列で--> |T, F| --をとる．    -->  |F|
  # |f11, f21|     |f14, f24|   出力される．    |T, T|                    |T|
  # (TおよびFは，それぞれ，TRUEおよびFALSE)
  # 
  # ベクトルの全ての要素がFであれば，個体1は，他の個体に支配されていないことになる．
  # 
  # （わざわざ，行列を用いて比較を行う理由は，処理速度を向上させるため．）
  for (i in ixa) {
    ixb <- setdiff(ixa, i)  # 比較対象Bのindex
    lenix <- length(ixb)

    a <- matrix(rep(objValues[i, ], lenix), nrow=lenix, ncol=numObj, byrow=TRUE)
    b <- objValues[ixb, ]

    if (opt == "min" || opt == 0) {
      cmpRes <- rowSums(a > b)
    }
    else if (opt == "max" || opt == 1) {
      cmpRes <- rowSums(a < b)
    }

    numNondominated <- sum(cmpRes >= numObj)  # 他の解に支配されている数

    if (numNondominated == 0) {             # どの解にも支配されていなければ，それは非劣解
      ixndom <- c(ixndom, i)
    }
    else {
      ixa <- setdiff(ixa, i)                # 支配されている解は，これ以降，比較する必要が無いため，比較対象から取り除く
    }
  }

  return(ixndom)
}

