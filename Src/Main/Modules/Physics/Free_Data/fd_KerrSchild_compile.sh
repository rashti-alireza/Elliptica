## compile Analytic Kerr-Schild free data:

#!/bin/bash

cpi fd_KerrSchild_c_k0_k1_k2_kt.math > fd_KerrSchild_c_k0_k1_k2_kt.c && \
  sed -i '/Welcome to Cpi/,$d' fd_KerrSchild_c_k0_k1_k2_kt.c && \
  sed -i -E -f fd_KerrSchild_replace_derives.sed \
               fd_KerrSchild_c_k0_k1_k2_kt.c \
& \
cpi fd_KerrSchild_g_dg_ddg_dddg.math > fd_KerrSchild_g_dg_ddg_dddg.c && \
  sed -i '/Welcome to Cpi/,$d' fd_KerrSchild_g_dg_ddg_dddg.c && \
  sed -i -E 's/\bk0\b\(/fd_ks_k0\(/g' fd_KerrSchild_g_dg_ddg_dddg.c && \
  sed -i -E 's/\bk1\b\(/fd_ks_k1\(/g' fd_KerrSchild_g_dg_ddg_dddg.c && \
  sed -i -E 's/\bk2\b\(/fd_ks_k2\(/g' fd_KerrSchild_g_dg_ddg_dddg.c && \
  sed -i -E 's/\bc\b\(/fd_ks_c\(/g'   fd_KerrSchild_g_dg_ddg_dddg.c && \
  sed -i -E -f fd_KerrSchild_replace_derives.sed \
               fd_KerrSchild_g_dg_ddg_dddg.c \
& \
cpi fd_KerrSchild_H_K0_K1_K2.math > fd_KerrSchild_H_K0_K1_K2.c && \
  sed -i '/Welcome to Cpi/,$d' fd_KerrSchild_H_K0_K1_K2.c && \
  sed -i -E 's/\bX\b\(/fd_ks_X\(/g' fd_KerrSchild_H_K0_K1_K2.c &&\
  sed -i -E 's/\bY\b\(/fd_ks_Y\(/g' fd_KerrSchild_H_K0_K1_K2.c &&\
  sed -i -E 's/\bZ\b\(/fd_ks_Z\(/g' fd_KerrSchild_H_K0_K1_K2.c &&\
  sed -i -E 's/\bR\b\(/fd_ks_R\(/g' fd_KerrSchild_H_K0_K1_K2.c &&\
  sed -i -E -f fd_KerrSchild_replace_derives.sed \
               fd_KerrSchild_H_K0_K1_K2.c \
& \
cpi fd_KerrSchild_R_X_Y_Z.math > fd_KerrSchild_R_X_Y_Z.c && \
  sed -i '/Welcome to Cpi/,$d' fd_KerrSchild_R_X_Y_Z.c && \
  sed -i -E 's/Blank\((\w+)\)/\1/g' fd_KerrSchild_R_X_Y_Z.c && \
  sed -i -E 's/\bXX\b\(/fd_ks_X\(/g' fd_KerrSchild_R_X_Y_Z.c && \
  sed -i -E 's/\bYY\b\(/fd_ks_Y\(/g' fd_KerrSchild_R_X_Y_Z.c && \
  sed -i -E 's/\bZZ\b\(/fd_ks_Z\(/g' fd_KerrSchild_R_X_Y_Z.c && \
  sed -i -E 's/\bList\b//g' fd_KerrSchild_R_X_Y_Z.c && \
  sed -i -E -f fd_KerrSchild_replace_derives.sed \
               fd_KerrSchild_R_X_Y_Z.c


               
