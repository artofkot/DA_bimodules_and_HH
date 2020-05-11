# -*- coding: utf-8 -*- 
import sys

from checking_da_T_LHD.arcslides_tensor_n_check import exec
from algebraic_structures.specific_algebras_modules_bimodules import torus_A
from algebraic_structures.basics import  Bunch_of_arrows,in_red
from algebraic_structures.da_bimodule import  (
    da_randomly_cancel_until_possible, are_equal_smart_da,
    da_check_df_is_0, composition)
from algebraic_structures.tensoring import (
    da_da_box_tensor_many_efficient_cancelations, 
    da_da_box_tensor_many_no_cancelations)
from algebraic_structures.hochschild_homology import is_bounded, CH, dimHH
from algebraic_structures.chain_complex import homology_dim
# from algebraic_structures.visual import draw_DA_bimodule, draw_chain_complex


from input_DA_bimodules import (
    ID1,ID2,ID3,M_RHD,M_LHD,L_RHD,L_LHD,
    g2_ID, g2_ID_bounded, g2_M_RHD, g2_M_LHD, g2_L_RHD,
    g2_L_LHD, g2_K_LHD, g2_K_RHD, g2_N_LHD, g2_N_RHD,
    g2_T_RHD, g2_T_LHD)
from itertools import permutations
import os
import timeit

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

########## SHOWCASES ##########

def showcase0():
    print("Here we initialize the main ten bimodules of the paper:")
    for Y in [A_,B_,C_,D_,E_,A_inv,B_inv,C_inv,D_inv,E_inv]:
        Y.show_short()
    print('\nType "yes" and press enter if you want more details about the bimodules; type "no" otherwise.')
    input1= input()
    if input1=='yes':
        for Y in [A_,B_,C_,D_,E_,A_inv,B_inv,C_inv,D_inv,E_inv]:
            Y.show() 
    print('To perform various sanity checks (by verifying the MCG relations) run the command "python3 run.py showcase3".')
    

def showcase1():
    print("\n5_1=Torus(5,2) knot has genus 2, and is fibered with monodromy τ_A᛫τ_B᛫τ_C᛫τ_D (see knotinfo.math.indiana.edu). Its knot Floer homology in the second to lowest Alexander grading has dimension 1 (see http://katlas.math.toronto.edu/wiki/Heegaard_Floer_Knot_Homology). Lets see this using Hochschild homology computation:")
    print("Computing HH(N(τ_A᛫τ_B᛫τ_C᛫τ_D))...")
    start_time = timeit.default_timer()
    X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_)
    Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X,g2_ID_bounded)
    print("dim(HH(N(τ_A᛫τ_B᛫τ_C᛫τ_D))=" + str(dimHH(Y)))
    elapsed = timeit.default_timer() - start_time
    print("Time it took to compute: "+str(elapsed)+"seconds")

    print("\n6_2 knot has genus 2, and is fibered with monodromy τ_A᛫τ_B᛫τ_C᛫τ_D^{-1}. Its knot Floer homology in the second to lowest Alexander grading has dimension 3. Lets see this using Hochschild homology computation:")
    print("Computing HH(N(τ_A᛫τ_B᛫τ_C᛫τ_D^{-1}))...")
    start_time = timeit.default_timer()
    X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_inv)
    Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X,g2_ID_bounded)
    print("dim(HH(N(τ_A᛫τ_B᛫τ_C᛫τ_D^{-1}))=" + str(dimHH(Y)))
    elapsed = timeit.default_timer() - start_time
    print("Time it took to compute: "+str(elapsed)+"seconds")

    print("\n7_6 knot has genus 2, and is fibered with monodromy τ_A᛫τ_B᛫τ_C^{-1}᛫τ_D. Its knot Floer homology in the second to lowest Alexander grading has dimension 5. Lets see this using Hochschild homology computation:")
    print("Computing HH(N(τ_A᛫τ_B᛫τ_C^{-1}᛫τ_D))...")
    start_time = timeit.default_timer()
    X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_inv,D_)
    Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X,g2_ID_bounded)
    print("dim(HH(N(τ_A᛫τ_B᛫τ_C^{-1}᛫τ_D^))=" + str(dimHH(Y)))
    elapsed = timeit.default_timer() - start_time
    print("Time it took to compute: "+str(elapsed)+"seconds")


def showcase2():
    print("ψ=τ_A᛫τ_B^{-1}")
    for n in range(5):
        print("\nComputing HH(N(ψ^"+str(n+1)+"))...")
        if n==4:
            print("It takes 57 seconds to finish the last computation. If you do not want to wait so long, you may press ctrl+c to stop the program.")
        start_time = timeit.default_timer()
        X=da_da_box_tensor_many_efficient_cancelations(A_,B_inv)
        Y=X
        for i in range(n):
            Y=da_da_box_tensor_many_efficient_cancelations(Y,X)
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        print("dim(HH(N(ψ^"+str(n+1)+")))=" + str(dimHH(Y)))
        elapsed = timeit.default_timer() - start_time
        print("Time it took to compute: "+str(elapsed)+"seconds")

def showcase3():
    X=da_da_box_tensor_many_efficient_cancelations(A_,A_inv)
    Z=da_da_box_tensor_many_efficient_cancelations(g2_ID)
    print("We first show that the tensor product of N(τ_A) with N(τ_A^{-1}) is isomorphic to the identity bimodule.\n")
    print("The following bimodule is N(τ_A᛫τ_A^{-1}), maximally canceled:")  
    X.show_short()
    print("The following bimodule is N(id):")  
    Z.show_short()
    are_equal_smart_da(X,g2_ID)
    print('\nType "yes" and press enter if you want more details about the bimodules and isomorphisms between them; type "no" to proceed to commuting relation.')
    input1= input()
    if input1=='yes':
        print("The following bimodule is N(τ_A᛫τ_C):")  
        X.show()
        print("The following bimodule is N(τ_C᛫τ_A):")  
        Z.show()
        are_equal_smart_da(X,Z,verbose=True)
        print("End of the first sanity check. Now off to the commuting relation:\n\n")  


    X=da_da_box_tensor_many_efficient_cancelations(A_,C_)
    Z=da_da_box_tensor_many_efficient_cancelations(C_,A_)
    print("Below is an example of bimodules satisfying commuting relation from MCG presentation (2.1) in the paper.\n")
    print("The following bimodule is N(τ_A᛫τ_C):")  
    X.show_short()
    print("The following bimodule is N(τ_C᛫τ_A):")  
    Z.show_short()
    are_equal_smart_da(X,Z)
    print('\nType "yes" and press enter if you want more details about the bimodules and isomorphisms between them; type "no"  to proceed to braid relation.')
    input1= input()
    if input1=='yes':
        print("The following bimodule is N(τ_A᛫τ_C):")  
        X.show()
        print("The following bimodule is N(τ_C᛫τ_A):")  
        Z.show()
        are_equal_smart_da(X,Z,verbose=True)
        print("End of commuting relation showcase. Now off to the braid relation:\n\n")  

    X=da_da_box_tensor_many_efficient_cancelations(A_,B_,A_)
    Z=da_da_box_tensor_many_efficient_cancelations(B_,A_,B_)
    print("Below is an example of bimodules satisfying braid relation from MCG presentation (2.1) in the paper.\n")
    print("The following bimodule is N(τ_A᛫τ_B᛫τ_A):")  
    X.show_short()
    print("The following bimodule is N(τ_B᛫τ_A᛫τ_B):")  
    Z.show_short()
    are_equal_smart_da(X,Z)
    print('\nType "yes" and press enter if you want more details about the bimodules and isomorphisms between them; type "no"  to proceed to the last relation.')
    input1= input()
    if input1=='yes':
        print("The following bimodule is N(τ_A᛫τ_B᛫τ_A):")  
        X.show()
        print("The following bimodule is N(τ_B᛫τ_A᛫τ_B):")  
        Z.show()
        are_equal_smart_da(X,Z,verbose=True)
        print("End of braid relation showcase. Now off to the last relation: (wait about 10 second)\n\n")  

    X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_,E_,E_,D_,C_,B_,A_)
    Z_=da_da_box_tensor_many_efficient_cancelations(E_,D_,C_,B_)
    Z=da_da_box_tensor_many_efficient_cancelations(Z_,Z_,Z_,Z_,Z_)
    print("We now check the last relation from MCG presentation (2.1) in the paper.\n")
    print("The following bimodule is N(τ_A᛫τ_B᛫τ_C᛫τ_D᛫τ_E᛫τ_E᛫τ_D᛫τ_C᛫τ_B᛫τ_A):")  
    X.show_short()
    print("The following bimodule is N((τ_E᛫τ_D᛫τ_C᛫τ_B)^5):")  
    Z.show_short()
    are_equal_smart_da(X,Z)
    print('\nType "yes" and press enter if you want more details about the bimodules and isomorphisms between them; type "no" and press enter otherwise.')
    input1= input()
    if input1=='yes':
        print("The following bimodule is N(τ_A᛫τ_B᛫τ_C᛫τ_D᛫τ_E᛫τ_E᛫τ_D᛫τ_C᛫τ_B᛫τ_A):")  
        X.show()
        print("The following bimodule is N((τ_E᛫τ_D᛫τ_C᛫τ_B)^5):")  
        Z.show()
        are_equal_smart_da(X,Z,verbose=True)

def showcase4():
    print("Here we perform Computation~4.2 from the paper, namely we compute Hochschild homology of the identity bimodule. For that we will use the bounded model for it, described in Section~2.3 of the paper.")
    print("This is the bimodule:")
    Y=da_da_box_tensor_many_no_cancelations(g2_ID)
    Y.show_short()
    print("This is the bounded model of the same bimodule:")
    Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded)
    Y.show_short()
    C=CH(Y)
    print("This is the Hochschild chain complex:")
    C.show_short()
    print("dim(HH(N(id))=" + str(homology_dim(C)))
    print('\nType "yes" and press enter if you want more details about the bimodule and its Hochschild homology; type "no" and press enter otherwise.')
    input1= input()
    if input1=='yes':
        print("This is the bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(g2_ID)
        Y.show()
        print("This is the bounded model of the same bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded)
        Y.show()
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show()
        print("\ndim(HH(N(id))=" + str(homology_dim(C)))

def showcase5():
    print("Here we perform Computation~4.3 from the paper. First we compute Hochschild homologies for all the five bimodules in question:\n")
    for Y in [A_,B_,C_,D_,E_]:
        print("This is the bimodule:")
        Y.show_short()
        print("This is the bounded model of the same bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        Y.show_short()
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show_short()
        print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    
    print('\nWe now provide more details, on demand:\nType "A" and press enter if you want more details about the bimodule τ_A and its Hochschild homology; \nType "B" and press enter if you want more details about the bimodule τ_B and its Hochschild homology; \nType "C" and press enter if you want more details about the bimodule τ_C and its Hochschild homology; \nType "D" and press enter if you want more details about the bimodule τ_D and its Hochschild homology; \nType "E" and press enter if you want more details about the bimodule τ_E and its Hochschild homology; \nType "no" and press enter otherwise.')
    input1= input()
    if input1 in ["A","B","C","D","E"]:
        print("This is the bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(eval(input1+'_'))
        Y.show()
        print("This is the bounded model of the same bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        Y.show()
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show()
        print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    
def showcase6():
    print("Here we perform Computation~4.3 from the paper. First we compute Hochschild homologies for all the five bimodules in question:\n")
    for Y in [A_inv,B_inv,C_inv,D_inv,E_inv]:
        print("This is the bimodule:")
        Y.show_short()
        print("This is the bounded model of the same bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        Y.show_short()
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show_short()
        print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    
    print('\nWe now provide more details, on demand:\nType "A" and press enter if you want more details about the bimodule τ_A^{-1} and its Hochschild homology; \nType "B" and press enter if you want more details about the bimodule τ_B^{-1} and its Hochschild homology; \nType "C" and press enter if you want more details about the bimodule τ_C^{-1} and its Hochschild homology; \nType "D" and press enter if you want more details about the bimodule τ_D^{-1} and its Hochschild homology; \nType "E" and press enter if you want more details about the bimodule τ_E^{-1} and its Hochschild homology; \nType "no" and press enter otherwise.')
    input1= input()
    if input1 in ["A","B","C","D","E"]:
        print("This is the bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(eval(input1+'_inv'))
        Y.show()
        print("This is the bounded model of the same bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        Y.show()
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show()
        print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    


def showcase7():
    print("Here we perform the last two computations from the proof of Theorem~1.2: (it may take about 4 minutes, after that we will output the results)\n")
    X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_)
    X=da_da_box_tensor_many_efficient_cancelations(X,X,X,X,X,X,X,X,X,X)
    Z=da_da_box_tensor_many_efficient_cancelations(A_,B_)
    Z=da_da_box_tensor_many_efficient_cancelations(Z,Z,Z,Z,Z,Z)
    for Y in [Z,X]:
        print("This is the bimodule:")
        Y.show_short()
        print("This is the bounded model of the same bimodule:")
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        Y.show_short()
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show_short()
        print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    print('\nWe now provide more details, on demand:\nType "yes" and press enter if you want more details; type "no" and press enter otherwise')
    input1= input()
    if input1=="yes":
        for Y in [Z,X]:
            print("This is the bimodule:")
            Y.show()
            Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
            C=CH(Y)
            print("This is the Hochschild chain complex:")
            C.show()
            print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    
def showcase8():
    print("Here we perform Computation~4.5 from the paper:")
    Y=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_)
    print("This is the bimodule:")
    Y.show_short()
    print("This is the bounded model of the same bimodule:")
    Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
    Y.show_short()
    C=CH(Y)
    print("This is the Hochschild chain complex:")
    C.show_short()
    print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    
    print('\nWe now provide more details, on demand:\nType "yes" and press enter if you want more details; type "no" and press enter otherwise')
    input1= input()
    if input1=="yes":
        Y=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_)
        print("This is the bimodule:")
        Y.show()
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show()
        print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')

def showcase9():
    print("Here we perform Computation~4.6 from the paper:")
    Y=da_da_box_tensor_many_efficient_cancelations(A_,A_,A_,A_,A_,B_,C_,D_,E_,E_,E_,E_,E_)
    print("This is the bimodule:")
    Y.show_short()
    print("This is the bounded model of the same bimodule:")
    Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
    Y.show_short()
    C=CH(Y)
    print("This is the Hochschild chain complex:")
    C.show_short()
    print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    
    print('\nWe now provide more details, on demand:\nType "yes" and press enter if you want more details; type "no" and press enter otherwise')
    input1= input()
    if input1=="yes":
        Y=da_da_box_tensor_many_efficient_cancelations(A_,A_,A_,A_,A_,B_,C_,D_,E_,E_,E_,E_,E_)
        print("This is the bimodule:")
        Y.show()
        Y=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,Y,g2_ID_bounded)
        C=CH(Y)
        print("This is the Hochschild chain complex:")
        C.show()
        print("Hochschild homology has dimension " + str(homology_dim(C))+'\n\n===========')
    
def showcase10():
    exec()

########## COMPUTATIONS ##########

# F2_THETA=composition(F2,THETA,A)
# F2_THETA.show()
# da_check_df_is_0(ID2,M_RHD,F2_THETA)



####################### genus two algebra computations 
# Bimodules: 
# g2_ID, g2_ID_bounded
# g2_K_RHD,g2_K_LHD
# g2_N_RHD,g2_N_LHD
# g2_T_RHD,g2_T_LHD
# g2_M_RHD,g2_M_LHD
# g2_L_RHD,g2_L_LHD


##### checking relations
# X=da_da_box_tensor_many_efficient_cancelations(A_,B_,C_,D_,E_,E_,D_,C_,B_,A_)
# Z_=da_da_box_tensor_many_efficient_cancelations(E_,D_,C_,B_)
# Z=da_da_box_tensor_many_efficient_cancelations(Z_,Z_,Z_,Z_,Z_)
# # X.show()
# print(are_equal_smart_da(X,Z))

##### computing HH 
# start_time = timeit.default_timer()
# X=da_da_box_tensor_many_efficient_cancelations(D_inv,C_inv,B_inv, A_inv)
# X=da_da_box_tensor_many_efficient_cancelations(X, X, X, X, X, X, X ,X ,X ,X)
# X=da_da_box_tensor_many_no_cancelations(g2_ID_bounded,X,g2_ID_bounded)
# X.show_short()
# print("\ndim(HH)=" + str(dimHH(X)))
# elapsed = timeit.default_timer() - start_time
# print elapsed
# HC=CH(X)
# HC.show()

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



func=sys.argv[1]
eval(func+'()')


