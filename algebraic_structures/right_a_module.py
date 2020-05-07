# -*- coding: utf-8 -*- 
from .basics import  check_idempotents_match_left_left, check_idempotents_match_right_left,check_idempotents_match_right_right
from .da_bimodule import Bunch_of_arrows

def right_a_in_mod_gen(right_a_arrow):
    return right_a_arrow[0]

def right_a_in_alg_tuple(right_a_arrow):
    return right_a_arrow[1]

def right_a_out_mod_gen(right_a_arrow):
    return right_a_arrow[2]


class Right_A_module(object):
    def __init__(self,gen_by_name,right_a_arrows,right_algebra,name,to_check=True):
        self.name=name
        self.gen_by_name=gen_by_name
        self.genset=self.gen_by_name.values()
        #differentials are represented by bunch of right_a_arrows with coefficients 1
        self.right_a_arrows=right_a_arrows
        self.right_a_arrows.delete_arrows_with_even_coeff()
        self.right_algebra=right_algebra

        if to_check==True:
            self.check()

    def check(self):
        ##Here we check that our " + self.name + " has all idempotents matching and d_squared=0:
        if not self.check_matching_of_idempotents_in_action():
            print( "\nSomething is wrong with idempotents!")
            raise NameError("Right_A_module " + self.name + " has problems with idempotents")
        
        d_squared=self.compute_d_squared()
        d_squared.delete_arrows_with_even_coeff()
        if d_squared:
            print( "d_squared is not 0! the terms that are not canceled:")
            d_squared.show()
            raise NameError("Right_A_module " + self.name + " doesn't satisfy d_squared=0.")

    def show(self):
        print( "==========")
        print( self.name + ':\n')
        print( 'Generators with their idempotents (' + str(len(self.genset))+ ' generators)')
        for gen in self.genset:
            print( str(gen) + '___' + str(gen.idem.right))

        print( '\nActions:')
        for right_a_arrow in self.right_a_arrows:
            print( str(right_a_arrow))

    # def show_for_tex(self):
    #     print( '=============')
    #     print( self.name + ':\n')
    #     print( str(len(self.genset)) + ' generators with their idempotents: \n')
    #     for gen in self.genset:
    #         # _{i_2}{(t_{12})}_{i_2}
    #         print( '${(' + gen.tex_name + ')}_{' + gen.idem.right.tex_name +'}$,')

    #     print( '\nActions:\n')
    #     for generator1 in self.genset:
    #         for generator2 in self.genset:
    #             arrows=[arrow for arrow in self.right_a_arrows if (right_a_in_mod_gen(arrow)==generator1 and right_a_out_mod_gen(arrow)==generator2)]
    #             if len(arrows)!=0:
    #                 for arrow in arrows:
    #                     print( '${}\\otimes {} \\longrightarrow {} $, '.format(getattr(right_a_in_mod_gen(arrow),'tex_name',1),getattr(right_a_in_mod_gen(arrow),'tex_name',1)))
       


    def check_matching_of_idempotents_in_action(self): 
        count_of_mismatches=0
        for right_a_arrow in self.right_a_arrows:
            # matching in last element of algebra and out gen
            if len(right_a_in_alg_tuple(right_a_arrow))!=0:
                if not check_idempotents_match_right_right(right_a_in_alg_tuple(right_a_arrow)[-1],right_a_out_mod_gen(right_a_arrow)): 
                    print( str(right_a_arrow) + '   idempotents are messed up in this right_a_arrow!1')
                    count_of_mismatches+=1
            
            # matching algebra elements between themselves, and with in generator
            for i in range(len(right_a_in_alg_tuple(right_a_arrow))):
                if (i==0):
                    #matching in algebra and gen
                    if not check_idempotents_match_right_left(right_a_in_mod_gen(right_a_arrow),right_a_in_alg_tuple(right_a_arrow)[0]): 
                        print( str(right_a_arrow) + '   idempotents are messed up in this right_a_arrow!1')
                        count_of_mismatches+=1  
                else:
                    #matching in algebras
                    if not check_idempotents_match_right_left(right_a_in_alg_tuple(right_a_arrow)[i-1],right_a_in_alg_tuple(right_a_arrow)[i]): 
                        print( right_a_arrow)
                        count_of_mismatches+=1
            # we don't check the following:
            # matchings between generators if differential is pure

        return count_of_mismatches==0

    def compute_d_squared(self):
        d_squared=Bunch_of_arrows()

        #contribution of double right_a_arrows
        for arrow1 in self.right_a_arrows:
            for arrow2 in self.right_a_arrows:
                if not right_a_out_mod_gen(arrow1)==right_a_in_mod_gen(arrow2): continue
                ar=(right_a_in_mod_gen(arrow1), right_a_in_alg_tuple(arrow1) + right_a_in_alg_tuple(arrow2)
                                ,right_a_out_mod_gen(arrow2))
                d_squared[ar]+=1


        #contribution of factorizing algebra elements , and differentials   
        for arrow in self.right_a_arrows:
            for index, a in enumerate(right_a_in_alg_tuple(arrow)):
                # factorizations
                for factorization in getattr(a,'factorizations', []):
                    new_tuple=right_a_in_alg_tuple(arrow)[:index] + factorization + right_a_in_alg_tuple(arrow)[index+1:]
                    ar=(right_a_in_mod_gen(arrow), new_tuple,
                        right_a_out_mod_gen(arrow))
                    d_squared[ar]+=1

                # predifferentials 
                for b in [algebra_diff_arrow[0] for algebra_diff_arrow in self.right_algebra.algebra_diff_arrows if algebra_diff_arrow[1]==a]:
                    new_tuple=right_a_in_alg_tuple(aa_arrow)[:index] + (b,) + right_a_in_alg_tuple(aa_arrow)[index+1:]
                    ar=(right_a_in_mod_gen(arrow), new_tuple,
                        right_a_out_mod_gen(arrow))
                    d_squared[ar]+=1

        return d_squared



