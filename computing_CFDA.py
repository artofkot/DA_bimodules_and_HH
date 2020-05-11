# -*- coding: utf-8 -*- 
import sys
# sys.path.append('../')

from algebraic_structures.specific_algebras_modules_bimodules import torus_A
from algebraic_structures.basics import  Bunch_of_arrows,in_red
from algebraic_structures.da_bimodule import  (
    da_randomly_cancel_until_possible, are_equal_smart_da,
    da_check_df_is_0, composition)
from algebraic_structures.tensoring import (
    da_da_box_tensor_many_efficient_cancelations, 
    da_da_box_tensor_many_no_cancelations)
from algebraic_structures.hochschild_homology import is_bounded, CH, dimHH
# from algebraic_structures.visual import draw_DA_bimodule, draw_chain_complex


from input_DA_bimodules import (
    ID1,ID2,ID3,M_RHD,M_LHD,L_RHD,L_LHD,
    g2_ID, g2_ID_bounded, g2_M_RHD, g2_M_LHD, g2_L_RHD,
    g2_L_LHD, g2_K_LHD, g2_K_RHD, g2_N_LHD, g2_N_RHD,
    g2_T_RHD, g2_T_LHD)
from itertools import permutations
import os
import timeit

####################### genus one algebra computations
# MORPHISMS
F2=Bunch_of_arrows([
        # from x
        (                         ID2.gen_by_name.x,(),
                1,ID1.gen_by_name.x),
        (                           ID2.gen_by_name.y,(),
                1,ID1.gen_by_name.y),
        (              ID2.gen_by_name.z2,(),
                torus_A.gen_by_name.r3,ID1.gen_by_name.y)
    ])
#morphism from Tova's paper
THETA=Bunch_of_arrows([
        # from x
        (                         ID1.gen_by_name.x,(),
                1,M_RHD.gen_by_name.p),
        (                           ID1.gen_by_name.x,(),
                torus_A.gen_by_name.r3,M_RHD.gen_by_name.r),
        (              ID1.gen_by_name.x,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r1,M_RHD.gen_by_name.r),
        (              ID1.gen_by_name.x,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r1,M_RHD.gen_by_name.q),
        # from y
        (              ID1.gen_by_name.y,(),
                1,M_RHD.gen_by_name.q),
        (              ID1.gen_by_name.y,(torus_A.gen_by_name.r2,),
                1,M_RHD.gen_by_name.r),
        (              ID1.gen_by_name.y,(torus_A.gen_by_name.r23,),
                1,M_RHD.gen_by_name.q)
    ])


########## COMPUTATIONS ##########

# F2_THETA=composition(F2,THETA,A)
# F2_THETA.show()
# da_check_df_is_0(ID2,M_RHD,F2_THETA)

# ########## Pseudo_Anosov mapping class
# start_time = timeit.default_timer()
# X=da_da_box_tensor_many_efficient_cancelations(L_RHD,M_LHD)
# X=da_da_box_tensor_many_efficient_cancelations(X,X,X)
# X=da_da_box_tensor_many_no_cancelations(ID3,X,ID3)
# print "\ndim(HH)=" + str(dimHH(X))
# elapsed = timeit.default_timer() - start_time
# print elapsed


########### Boundary Dehn twist
# AB=da_da_box_tensor_many_no_cancelations(M_RHD,L_RHD)
# BD=da_randomly_cancel_until_possible(da_da_box_tensor_many_no_cancelations(AB,AB,AB))
# BD=da_randomly_cancel_until_possible(da_da_box_tensor_many_no_cancelations(BD,BD))
# BD_bounded=da_da_box_tensor_many_no_cancelations(ID3,BD,ID3)
# # draw_DA_bimodule(BD)
# print "\ndim(HH)=" + str(dimHH(BD_bounded))

########### Any
# X=da_da_box_tensor_many_efficient_cancelations(L_RHD,M_LHD)
# X=da_da_box_tensor_many_efficient_cancelations(X,X,X)
# X=da_da_box_tensor_many_no_cancelations(ID3,X,ID3)
# HC=CH(X)
# HC.show()
# print is_bounded(X)
# print "\n2 BASEPOINT CASE(+1 added here)! dim(HH)=" + str(dimHH(X)+1)
# print "\ndim(HH)=" + str(dimHH(X))


####################### genus two algebra computations 
# Bimodules: 
# g2_ID, g2_ID_bounded
# g2_K_RHD,g2_K_LHD
# g2_N_RHD,g2_N_LHD
# g2_T_RHD,g2_T_LHD
# g2_M_RHD,g2_M_LHD
# g2_L_RHD,g2_L_LHD

A_=g2_K_RHD
B_=g2_N_RHD
C_=g2_T_RHD
D_=g2_L_RHD
E_=g2_M_RHD
A_inv=g2_K_LHD
B_inv=g2_N_LHD
C_inv=g2_T_LHD
D_inv=g2_L_LHD
E_inv=g2_M_LHD

##### checking relations
# X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_,E_,E_,D_,C_,B_,A_)
# Z_=da_da_box_tensor_many_efficient_cancelations(E_,D_,C_,B_)
# Z=da_da_box_tensor_many_efficient_cancelations(Z_,Z_,Z_,Z_,Z_)
# X.show()
# print(are_equal_smart_da(X,Z))

##### computing HH 
start_time = timeit.default_timer()
X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_,C_,A_,B_,E_)
X=da_da_box_tensor_many_efficient_cancelations(X,X,X)
X=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X,g2_ID_bounded)
X.show_short()
print("\ndim(HH)=" + str(dimHH(X)))
elapsed = timeit.default_timer() - start_time
print(elapsed)

##### experiment, that shows that order of elements in the product matters for HH
##### (cyclic order doesnt matter due to conjugation invariance)
# a=[A_,B_,C_,A_,B_,C_]
# perms=permutations(a,len(a))
# for perm in perms:
#     X=g2_ID
#     for bim in perm:
#         X=da_da_box_tensor_many_efficient_cancelations(X,bim)
#     X=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X,g2_ID_bounded)
#     print "\ndim(HH)=" + str(dimHH(X))

##### experiment, that doesn't show dependence on substituting RHD by LHD
# a=[(A_,A_inv),(B_,B_inv),(C_,C_inv),(A_,A_inv),(B_,B_inv)]
# perms=permutations(a,len(a))
# for perm in perms:
#     X1=g2_ID
#     X2=g2_ID
#     for bim in perm:
#         X1=da_da_box_tensor_many_efficient_cancelations(X1,bim[0])
#         X2=da_da_box_tensor_many_efficient_cancelations(X2,bim[1])
#     X1=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X1,g2_ID_bounded)
#     X2=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X2,g2_ID_bounded)
#     print "\ndim(HH)=" + str(dimHH(X1)),
#     print " dim(HH)=" + str(dimHH(X2))

#### experiment, where one computes possible ranks
# k=0
# a=[A_,A_inv,B_,B_inv,C_,C_inv,A_,A_inv,B_,B_inv,g2_ID]
# for i in range(11):
#     for j in range(11):
#         for k in range(11):
#                 X=da_da_box_tensor_many_efficient_cancelations(a[i],a[j],a[k])
#                 Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X,g2_ID_bounded)
#                 hf=dimHH(Y)
#                 if hf==1: 
#                     in_red('YAY, 1 dimensional homology!')
#                     i.show_short()
#                     j.show_short()
#                     k.show_short()
#                     l.show_short()
#                     k=k+1
#                 else:
#                     print 'Ok: homology is {}'.format(str(hf))
# print k





