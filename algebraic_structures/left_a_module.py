# -*- coding: utf-8 -*- 
from .basics import check_idempotents_match_left_left, check_idempotents_match_right_left,check_idempotents_match_right_right
from .da_bimodule import Bunch_of_arrows
from random import shuffle

def left_a_in_mod_gen(left_a_arrow):
    return left_a_arrow[1]

def left_a_in_alg_tuple(left_a_arrow):
    return left_a_arrow[0]

def left_a_out_mod_gen(left_a_arrow):
    return left_a_arrow[2]


class Left_A_module(object):
    def __init__(self,gen_by_name,left_a_arrows,left_algebra,name,to_check=True):
        self.name=name
        self.gen_by_name=gen_by_name
        self.genset=self.gen_by_name.values()
        #differentials are represented by bunch of left_a_arrows with coefficients 1
        self.left_a_arrows=left_a_arrows
        self.left_a_arrows.delete_arrows_with_even_coeff()
        self.left_algebra=left_algebra

        if to_check==True:
            self.check()

    def check(self):
        ##Here we check that our " + self.name + " has all idempotents matching and d_squared=0:
        if not self.check_matching_of_idempotents_in_action():
            print( "\nSomething is wrong with idempotents!")
            raise NameError("Left_A_module " + self.name + " has problems with idempotents")
        
        d_squared=self.compute_d_squared()
        d_squared.delete_arrows_with_even_coeff()
        if d_squared:
            print( "d_squared is not 0! the terms that are not canceled:")
            d_squared.show()
            raise NameError("Left_A_module " + self.name + " doesn't satisfy d_squared=0.")

    def show(self):
        print( "==========")
        print( self.name + ':\n')
        print( 'Generators with their idempotents (' + str(len(self.genset))+ ' generators)')
        for gen in self.genset:
            print( str(gen.idem.left) + '___'  + str(gen))

        print( '\nActions:')
        for left_a_arrow in self.left_a_arrows:
            print( str(left_a_arrow))


    def check_matching_of_idempotents_in_action(self): 
        count_of_mismatches=0
        for left_a_arrow in self.left_a_arrows:

            # matching in last element of algebra and in gen
            if len(left_a_in_alg_tuple(left_a_arrow))!=0:
                if not check_idempotents_match_right_left(left_a_in_alg_tuple(left_a_arrow)[-1],left_a_in_mod_gen(left_a_arrow)): 
                    print( str(left_a_arrow) + '   idempotents are messed up in this left_a_arrow!1')
                    count_of_mismatches+=1
            
            # matching algebra elements between themselves, and with in generator
            for i in range(len(left_a_in_alg_tuple(left_a_arrow))):
                if (i==0):
                    #matching in algebra and gen
                    if not check_idempotents_match_left_left(left_a_out_mod_gen(left_a_arrow),left_a_in_alg_tuple(left_a_arrow)[0]): 
                        print( str(left_a_arrow) + '   idempotents are messed up in this left_a_arrow!1')
                        count_of_mismatches+=1  
                else:
                    #matching in algebras
                    if not check_idempotents_match_right_left(left_a_in_alg_tuple(left_a_arrow)[i-1],left_a_in_alg_tuple(left_a_arrow)[i]): 
                        print( left_a_arrow)
                        count_of_mismatches+=1

            # we don't check the following:
            # matchings between generators if differential is pure
            
        return count_of_mismatches==0

    def compute_d_squared(self):
        d_squared=Bunch_of_arrows()

        #contribution of double left_a_arrows
        for arrow1 in self.left_a_arrows:
            for arrow2 in self.left_a_arrows:
                if not left_a_out_mod_gen(arrow1)==left_a_in_mod_gen(arrow2): continue
                ar=(left_a_in_alg_tuple(arrow2) + left_a_in_alg_tuple(arrow1), left_a_in_mod_gen(arrow1)
                                ,left_a_out_mod_gen(arrow2))
                d_squared[ar]+=1


        #contribution of factorizing algebra elements, and differentials        
        for arrow in self.left_a_arrows:
            for index, a in enumerate(left_a_in_alg_tuple(arrow)):
                # factorizations
                for factorization in getattr(a,'factorizations', []):
                    new_tuple=left_a_in_alg_tuple(arrow)[:index] + factorization + left_a_in_alg_tuple(arrow)[index+1:]
                    ar=( new_tuple,left_a_in_mod_gen(arrow),
                        left_a_out_mod_gen(arrow))
                    d_squared[ar]+=1

                # predifferentials 
                for b in [algebra_diff_arrow[0] for algebra_diff_arrow in self.left_algebra.algebra_diff_arrows if algebra_diff_arrow[1]==a]:
                    new_tuple=left_a_in_alg_tuple(aa_arrow)[:index] + (b,) + left_a_in_alg_tuple(aa_arrow)[index+1:]
                    ar=( new_tuple,left_a_in_mod_gen(arrow),
                        left_a_out_mod_gen(arrow))
                    d_squared[ar]+=1

        return d_squared

def left_a_cancel_pure_differential(left_A_module_old,pure_differential):
    #z1--->z2
    z1=left_a_in_mod_gen(pure_differential)
    z2=left_a_out_mod_gen(pure_differential)

    #form generators
    generators_of_new_left_A=left_A_module_old.gen_by_name
    del generators_of_new_left_A[z1.name]
    del generators_of_new_left_A[z2.name]

    #form old differentials
    old_arrows_that_survive=[arrow for arrow in left_A_module_old.left_a_arrows if (left_a_out_mod_gen(arrow)!=z1 and left_a_out_mod_gen(arrow)!=z2  and left_a_in_mod_gen(arrow)!=z1 and left_a_in_mod_gen(arrow)!=z2 ) ]
    arrows_in_new_left_A=Bunch_of_arrows(old_arrows_that_survive)

    #form new differentials
    arrows_in_z2=[arrow for arrow in left_A_module_old.left_a_arrows if (left_a_out_mod_gen(arrow)==z2 and left_a_in_mod_gen(arrow)!=z1 and left_a_in_mod_gen(arrow)!=z2)]
    arrows_from_z1=[arrow for arrow in left_A_module_old.left_a_arrows if (left_a_in_mod_gen(arrow)==z1 and left_a_out_mod_gen(arrow)!=z2 and left_a_out_mod_gen(arrow)!=z1)]

    for arrow_in_z2 in arrows_in_z2:
        for arrow_from_z1 in arrows_from_z1:
            new_arrow=(left_a_in_alg_tuple(arrow_from_z1)+left_a_in_alg_tuple(arrow_in_z2),left_a_in_mod_gen(arrow_in_z2),
                        left_a_out_mod_gen(arrow_from_z1) )
            arrows_in_new_left_A[new_arrow]+=1


    arrows_in_new_left_A.delete_arrows_with_even_coeff()
    return Left_A_module(generators_of_new_left_A,arrows_in_new_left_A,left_A_module_old.left_algebra,name= left_A_module_old.name +'_red')

def left_a_randomly_cancel_until_possible(left_A):
    there_is_diff=0
    # arrs=left_A.left_a_arrows
    arrs=list(left_A.left_a_arrows)
    shuffle(arrs)
    for arrow in arrs:
        if (left_a_in_alg_tuple(arrow)==() and left_a_in_mod_gen(arrow)!=left_a_out_mod_gen(arrow)):
            # now we need to check that there are no other arrows between these two generators
            other_arrows=[ar for ar in arrs if (left_a_in_mod_gen(ar)==left_a_in_mod_gen(arrow) and left_a_out_mod_gen(ar)==left_a_out_mod_gen(arrow) and ar!=arrow)]
            if other_arrows: continue
            else:
                there_is_diff=1
                canceled_AA=left_a_cancel_pure_differential(left_A,arrow)
                return (left_a_randomly_cancel_until_possible(canceled_AA))

    if there_is_diff==0:
        return left_A


