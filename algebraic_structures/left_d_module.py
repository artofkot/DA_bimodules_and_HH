# -*- coding: utf-8 -*- 
from .chain_complex import ChainComplex
from .basics import AttrDict, debug
from .basics import check_idempotents_match_left_left, check_idempotents_match_right_left,check_idempotents_match_right_right
from .basics import Bunch_of_arrows
from random import shuffle
from .algebra import Generator, dg_algebra
from itertools import product

def left_d_in_mod_gen(left_d_arrow):
    return left_d_arrow[0]

def left_d_out_alg_gen(left_d_arrow):
    return left_d_arrow[1]

def left_d_out_mod_gen(left_d_arrow):
    return left_d_arrow[2]


class Left_D_module(object):
    def __init__(self,gen_by_name,left_d_arrows,left_algebra,name,to_check=True):
        self.name=name
        self.gen_by_name=gen_by_name
        self.genset=self.gen_by_name.values()
        #differentials are represented by bunch of left_d_arrows with coefficients 1
        self.left_d_arrows=left_d_arrows
        self.left_d_arrows.delete_arrows_with_even_coeff()
        self.left_algebra=left_algebra

        if to_check==True:
            self.check()

    def check(self):
        ##Here we check that our " + self.name + " has all idempotents matching and d_squared=0:
        if not self.check_matching_of_idempotents_in_action():
            print( "\nSomething is wrong with idempotents!")
            raise NameError("Left_D_module " + self.name + " has problems with idempotents")
        
        d_squared=self.compute_d_squared()
        d_squared.delete_arrows_with_even_coeff()
        if d_squared:
            print( "d_squared is not 0! the terms that are not canceled:")
            d_squared.show()
            raise NameError("Left_D_module " + self.name + " doesn't satisfy d_squared=0.")

    def show(self):
        print( "==========")
        print( self.name + ':\n')
        print( 'Generators with their idempotents (' + str(len(self.genset))+ ' generators)')
        for gen in self.genset:
            print( str(gen.idem.left) + '___'  + str(gen))

        print( '\n' + str(len(self.left_d_arrows)) + ' actions:')
        for left_d_arrow in self.left_d_arrows:
            print( str(left_d_arrow))


    def check_matching_of_idempotents_in_action(self): 
        count_of_mismatches=0
        for left_d_arrow in self.left_d_arrows:
            # matching out alg gen and out mod gen
            if not check_idempotents_match_right_left(left_d_out_alg_gen(left_d_arrow),left_d_out_mod_gen(left_d_arrow)): 
                print( str(left_d_arrow) + '   idempotents are messed up in this left_d_arrow!1')
                count_of_mismatches+=1
            
            # matching out alg gen and in mod gen
            if not check_idempotents_match_left_left(left_d_out_alg_gen(left_d_arrow),left_d_in_mod_gen(left_d_arrow)): 
                print( str(left_d_arrow) + '   idempotents are messed up in this left_d_arrow!1')
                count_of_mismatches+=1

            # we don't check the following:
            # matchings between generators if differential is pure
            
        return count_of_mismatches==0

    def compute_d_squared(self):
        d_squared=Bunch_of_arrows()

        #contribution of differential on algebra element
        for arrow in self.left_d_arrows:
            a = left_d_out_alg_gen(arrow)
            if a:
                for b in [algebra_diff_arrow[1] for algebra_diff_arrow in self.left_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                    ar=(left_d_in_mod_gen(arrow),
                        b,left_d_out_mod_gen(arrow))
                    d_squared[ar]+=1

        #contribution of double left_d_arrows
        for arrow1 in self.left_d_arrows:
            for arrow2 in self.left_d_arrows:
                if not left_d_out_mod_gen(arrow1)==left_d_in_mod_gen(arrow2): continue
                new_alg_gen=self.left_algebra.multiply(left_d_out_alg_gen(arrow1),left_d_out_alg_gen(arrow2))
                if not new_alg_gen: continue
                ar=(left_d_in_mod_gen(arrow1),
                    new_alg_gen,left_d_out_mod_gen(arrow2))
                d_squared[ar]+=1

        return d_squared

def left_d_cancel_pure_differential(Left_D_module_old,pure_differential):
    #z1--->z2
    z1=left_d_in_mod_gen(pure_differential)
    z2=left_d_out_mod_gen(pure_differential)

    #form generators
    generators_of_new_Left_D=Left_D_module_old.gen_by_name
    del generators_of_new_Left_D[z1.name]
    del generators_of_new_Left_D[z2.name]

    #form old differentials
    old_arrows_that_survive=[arrow for arrow in Left_D_module_old.left_d_arrows if (left_d_out_mod_gen(arrow)!=z1 and left_d_out_mod_gen(arrow)!=z2  and left_d_in_mod_gen(arrow)!=z1 and left_d_in_mod_gen(arrow)!=z2 ) ]
    arrows_in_new_Left_D=Bunch_of_arrows(old_arrows_that_survive)

    #form new differentials
    arrows_in_z2=[arrow for arrow in Left_D_module_old.left_d_arrows if (left_d_out_mod_gen(arrow)==z2 and left_d_in_mod_gen(arrow)!=z1 and left_d_in_mod_gen(arrow)!=z2)]
    arrows_from_z1=[arrow for arrow in Left_D_module_old.left_d_arrows if (left_d_in_mod_gen(arrow)==z1 and left_d_out_mod_gen(arrow)!=z2 and left_d_out_mod_gen(arrow)!=z1)]

    for arrow_in_z2 in arrows_in_z2:
        for arrow_from_z1 in arrows_from_z1:

            new_out_left_alg_gen=Left_D_module_old.left_algebra.multiply(left_d_out_alg_gen(arrow_in_z2),left_d_out_alg_gen(arrow_from_z1))

            if new_out_left_alg_gen:
                new_arrow=(left_d_in_mod_gen(arrow_in_z2),
                    new_out_left_alg_gen,left_d_out_mod_gen(arrow_from_z1))
                arrows_in_new_Left_D[new_arrow]+=1


    arrows_in_new_Left_D.delete_arrows_with_even_coeff()

    return Left_D_module(generators_of_new_Left_D,arrows_in_new_Left_D,Left_D_module_old.left_algebra,name= Left_D_module_old.name +'_red')


def left_d_randomly_cancel_until_possible(Left_D):
    there_is_diff=0
    arrs=list(Left_D.left_d_arrows)
    shuffle(arrs)
    for arrow in arrs:
        if (left_d_out_alg_gen(arrow)==1 and left_d_in_mod_gen(arrow)!=left_d_out_mod_gen(arrow)):
            # now we need to check that there are no other arrows between these two generators
            other_arrows=[ar for ar in arrs if (left_d_in_mod_gen(ar)==left_d_in_mod_gen(arrow) and left_d_out_mod_gen(ar)==left_d_out_mod_gen(arrow) and ar!=arrow)]
            if other_arrows: continue
            else:
                there_is_diff=1
                canceled_Left_D=left_d_cancel_pure_differential(Left_D,arrow)
                return (left_d_randomly_cancel_until_possible(canceled_Left_D))

    if there_is_diff==0:
        return Left_D

def morphism_space_for_left_D_structures(D1,D2):
    generators_of_mor_space=AttrDict({})
    A=D1.left_algebra
    for genD1 in D1.genset:
        for a in A.genset:
            if not a.idem.left==genD1.idem.left: continue
            for genD2 in D2.genset:
                if not a.idem.right==genD2.idem.left: continue
                f=(genD1,a,genD2)
                generators_of_mor_space[str(f)]=Generator(name=str(f), aux_info=f)
    arrows_in_mor_space=Bunch_of_arrows()
    for morphism_gen in generators_of_mor_space.values():
        f=morphism_gen.aux_info
        
        for arrow1 in [ar for ar in D1.left_d_arrows if left_d_out_mod_gen(ar)==left_d_in_mod_gen(f)]:
            a_new=A.multiply( left_d_out_alg_gen(arrow1),left_d_out_alg_gen(f) )
            if a_new:
                df=(left_d_in_mod_gen(arrow1),a_new,left_d_out_mod_gen(f))
                arrows_in_mor_space[(generators_of_mor_space[str(f)],generators_of_mor_space[str(df)])]+=1


        for arrow2 in [ar for ar in D2.left_d_arrows if left_d_in_mod_gen(ar)==left_d_out_mod_gen(f)]:

            a_new=A.multiply( left_d_out_alg_gen(f),left_d_out_alg_gen(arrow2) )
            if a_new:
                df=(left_d_in_mod_gen(f),a_new,left_d_out_mod_gen(arrow2))
                arrows_in_mor_space[(generators_of_mor_space[str(f)],generators_of_mor_space[str(df)])]+=1
    
    arrows_in_mor_space.delete_arrows_with_even_coeff()

    return ChainComplex(generators_of_mor_space,arrows_in_mor_space,name= 'Mor('+D1.name+','+D2.name +')',to_check=True)
