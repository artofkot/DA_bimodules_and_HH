# -*- coding: utf-8 -*- 
from .basics import check_idempotents_match_left_left, check_idempotents_match_right_left,check_idempotents_match_right_right
from .basics import Bunch_of_arrows

def right_d_in_mod_gen(right_d_arrow):
    return right_d_arrow[0]

def right_d_out_alg_gen(right_d_arrow):
    return right_d_arrow[2]

def right_d_out_mod_gen(right_d_arrow):
    return right_d_arrow[1]


class Right_D_module(object):
    def __init__(self,gen_by_name,right_d_arrows,right_algebra,name,to_check=True):
        self.name=name
        self.gen_by_name=gen_by_name
        self.genset=self.gen_by_name.values()
        #differentials are represented by bunch of right_d_arrows with coefficients 1
        self.right_d_arrows=right_d_arrows
        self.right_d_arrows.delete_arrows_with_even_coeff()
        self.right_algebra=right_algebra

        if to_check==True:
            self.check()

    def check(self):
        ##Here we check that our " + self.name + " has all idempotents matching and d_squared=0:
        if not self.check_matching_of_idempotents_in_action():
            print( "\nSomething is wrong with idempotents!")
            raise NameError("Right_D_module " + self.name + " has problems with idempotents")
        
        d_squared=self.compute_d_squared()
        d_squared.delete_arrows_with_even_coeff()
        if d_squared:
            print( "d_squared is not 0! the terms that are not canceled:")
            d_squared.show()
            raise NameError("Right_D_module " + self.name + " doesn't satisfy d_squared=0.")

    def show(self):
        print( "==========")
        print( self.name + ':\n')
        print( 'Generators with their idempotents (' + str(len(self.genset))+ ' generators)')
        for gen in self.genset:
            print( str(gen) + '___'  + str(gen.idem.right))

        print( '\nActions:')
        for right_d_arrow in self.right_d_arrows:
            print( str(right_d_arrow))


    def check_matching_of_idempotents_in_action(self): 
        count_of_mismatches=0
        for right_d_arrow in self.right_d_arrows:
            # matching out alg gen and out mod gen
            if not check_idempotents_match_right_right(right_d_out_alg_gen(right_d_arrow),right_d_in_mod_gen(right_d_arrow)): 
                print( str(right_d_arrow) + '   idempotents are messed up in this right_d_arrow!1')
                count_of_mismatches+=1
            
            # matching out alg gen and in mod gen
            if not check_idempotents_match_right_left(right_d_out_mod_gen(right_d_arrow),right_d_out_alg_gen(right_d_arrow)): 
                print( str(right_d_arrow) + '   idempotents are messed up in this right_d_arrow!1')
                count_of_mismatches+=1

            # matching out alg gen and in mod gen
            if not check_idempotents_match_left_left(right_d_out_mod_gen(right_d_arrow),right_d_in_mod_gen(right_d_arrow)): 
                print( str(right_d_arrow) + '   idempotents are messed up in this right_d_arrow!1')
                count_of_mismatches+=1

            # we don't check the following:
            # matchings between generators if differential is pure
            
        return count_of_mismatches==0

    def compute_d_squared(self):
        d_squared=Bunch_of_arrows()

        #contribution of differential on algebra element
        for arrow in self.right_d_arrows:
            a = right_d_out_alg_gen(arrow)
            if a:
                for b in [algebra_diff_arrow[1] for algebra_diff_arrow in self.right_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                    ar=(right_d_in_mod_gen(arrow),
                        right_d_out_mod_gen(arrow),b)
                    d_squared[ar]+=1

        #contribution of double right_d_arrows
        for arrow1 in self.right_d_arrows:
            for arrow2 in self.right_d_arrows:
                if not right_d_out_mod_gen(arrow1)==right_d_in_mod_gen(arrow2): continue
                new_alg_gen=self.right_algebra.multiply(right_d_out_alg_gen(arrow2),right_d_out_alg_gen(arrow1))
                if not new_alg_gen: continue
                ar=(right_d_in_mod_gen(arrow1),
                    right_d_out_mod_gen(arrow2),new_alg_gen)
                d_squared[ar]+=1

        return d_squared



