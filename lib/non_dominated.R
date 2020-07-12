non_dominated <- function (objValues, opt) {
  # ���̊֐��́C�ړI�֐��̒l���i�[����Ă���s��objValues��
  # �����āC�����ƂȂ�s�ԍ����i�[����Ă���x�N�g��ixndom��Ԃ��D
  # 
  # ����:
  #   objValues�́C���̂悤�ȍs���z�肵�Ă���D
  # 
  #               |f11, f21, ..., fk1|
  #   objValues = |f12, f22, ..., fk2|
  #               | : ,  : , ...,  : |
  #               |f1p, f2p, ..., fkp|
  # 
  #   fkp�́Cp�Ԗڂ̌̂ɂ�����k�Ԗڂ̖ړI�֐��̒l�D
  # 
  #   opt�́C�ŏ������̏ꍇ��"min"�ɁC�ő剻���̏ꍇ��"max"�ɐݒ肷��D
  # 
  # �߂�l�F
  #   objValues�̂Ȃ��ŁC�����ƂȂ�̂̍s�ԍ��x�N�g����Ԃ��D

  numObj <- ncol(objValues)  # �ړI�֐��̐�
  numPop <- nrow(objValues)  # �W�c�̐�

  ixa <- 1:numPop            # ��r�Ώ�A��index
  ixndom <- NULL             # ������index��ۑ�����x�N�g��

  # 2�ړI�ŏ������̏ꍇ�ɂ�����C�����̒��o���@ 
  # [��1]
  #    �s��a          �s��b                      �s��                   �x�N�g��
  # |f11, f21|     |f12, f22|   �e�����̔�r    |F, F|   �e�s�̘_����     |F|
  # |f11, f21|  >  |f13, f23| --���ʂ��s���--> |T, F| --���Ƃ�D    -->  |F|
  # |f11, f21|     |f14, f24|   �o�͂����D    |T, T|                    |T|
  # (T�����F�́C���ꂼ��CTRUE�����FALSE)
  # 
  # �x�N�g���̑S�Ă̗v�f��F�ł���΁C��1�́C���̌̂Ɏx�z����Ă��Ȃ����ƂɂȂ�D
  # 
  # �i�킴�킴�C�s���p���Ĕ�r���s�����R�́C�������x�����コ���邽�߁D�j
  for (i in ixa) {
    ixb <- setdiff(ixa, i)  # ��r�Ώ�B��index
    lenix <- length(ixb)

    a <- matrix(rep(objValues[i, ], lenix), nrow=lenix, ncol=numObj, byrow=TRUE)
    b <- objValues[ixb, ]

    if (opt == "min" || opt == 0) {
      cmpRes <- rowSums(a > b)
    }
    else if (opt == "max" || opt == 1) {
      cmpRes <- rowSums(a < b)
    }

    numNondominated <- sum(cmpRes >= numObj)  # ���̉��Ɏx�z����Ă��鐔

    if (numNondominated == 0) {             # �ǂ̉��ɂ��x�z����Ă��Ȃ���΁C����͔���
      ixndom <- c(ixndom, i)
    }
    else {
      ixa <- setdiff(ixa, i)                # �x�z����Ă�����́C����ȍ~�C��r����K�v���������߁C��r�Ώۂ����菜��
    }
  }

  return(ixndom)
}
