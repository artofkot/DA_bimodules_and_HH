# -*- coding: utf-8 -*- 
from .basics import ( debug,
    in_red,check_idempotents_match_left_left, 
    check_idempotents_match_right_left,
    check_idempotents_match_right_right)
from .da_bimodule import Bunch_of_arrows
from random import shuffle

# def dd_arrow_to_str(tuplee): #works only for DD arrows
    # if len(tuplee)==4:
    #     return str(tuplee[0]) +'⊗'+ str(tuplee[1]) + "---->" + str(tuplee[2]) +'⊗'+ str(tuplee[3])  

def dd_in_mod_gen(dd_arrow):
    return dd_arrow[0]

def dd_out_left_alg_gen(dd_arrow):
    return dd_arrow[1]

def dd_out_mod_gen(dd_arrow):
    return dd_arrow[2]

def dd_out_right_alg_gen(dd_arrow):
    return dd_arrow[3]


class DD_bimodule(object):
    def __init__(self,gen_by_name,dd_arrows,left_algebra,right_algebra,name,to_check=True):
        self.name=name
        self.gen_by_name=gen_by_name
        self.genset=self.gen_by_name.values()
        #differentials are represented by bunch of dd_arrows with coefficients 1
        self.dd_arrows=dd_arrows
        self.dd_arrows.delete_arrows_with_even_coeff()
        self.left_algebra=left_algebra
        self.right_algebra=right_algebra

        if to_check==True:
            self.check()

    def check(self):
        ##Here we check that our " + self.name + " has all idempotents matching and d_squared=0:
        if not self.check_matching_of_idempotents_in_action():
            print ("\nSomething is wrong with idempotents!")
            raise NameError("DD_bimodule " + self.name + " has problems with idempotents")
        
        d_squared=self.compute_d_squared()
        d_squared.delete_arrows_with_even_coeff()
        if d_squared:
            print ("d_squared is not 0! the terms that are not canceled:")
            d_squared.show()
            raise NameError("DD_bimodule " + self.name + " doesn't satisfy d_squared=0.")

    def show(self):
        print ("==========")
        print (self.name + ':\n')
        print ('{} generators with their idempotents:'.format(str(len(self.genset))))
        for gen in self.genset:
            print( str(gen.idem.left) + '___' + str(gen) + '___' + str(gen.idem.right))

        print( '\n{} actions:'.format(len(list(self.dd_arrows))))
        for dd_arrow in self.dd_arrows:
            print( str(dd_arrow))

    def show_for_tex(self):
        print( '=============')
        print( self.name + ':\n')
        print( str(len(self.genset)) + ' generators with their idempotents: \n')
        for gen in self.genset:
            # _{i_2}{(t_{12})}_{i_2}
            print( '$_{' + gen.idem.left.tex_name + '}{(' + gen.tex_name + ')}_{' + gen.idem.right.tex_name +'}$,')

        print( '\nActions:\n')
        for generator1 in self.genset:
            for generator2 in self.genset:
                arrows=[arrow for arrow in self.dd_arrows if (dd_in_mod_gen(arrow)==generator1 and dd_out_mod_gen(arrow)==generator2)]
                if len(arrows)!=0:
                    for arrow in arrows:
                        print( '${}\\rightarrow {} \\otimes {} \\otimes {} $, '.format(getattr(dd_in_mod_gen(arrow),'tex_name',1),getattr(dd_out_left_alg_gen(arrow),'tex_name',1),getattr(dd_out_mod_gen(arrow),'tex_name',1),getattr(dd_out_right_alg_gen(arrow),'tex_name',1)))
                
            


    def check_matching_of_idempotents_in_action(self): 
        count_of_mismatches=0
        for dd_arrow in self.dd_arrows:

            # matching out left algebra and out gen
            if not check_idempotents_match_right_left( dd_out_left_alg_gen(dd_arrow),dd_out_mod_gen(dd_arrow) ): 
                print( str(dd_arrow) + '   idempotents are messed up in this dd_arrow!1')
                count_of_mismatches+=1

            # matching out left algebra and in gen
            if not check_idempotents_match_left_left( dd_in_mod_gen(dd_arrow),dd_out_left_alg_gen(dd_arrow) ): 
                print( str(dd_arrow) + '   idempotents are messed up in this dd_arrow!1')
                count_of_mismatches+=1

            # matching out right algebra and out gen
            if not check_idempotents_match_right_left(dd_out_mod_gen(dd_arrow),dd_out_right_alg_gen(dd_arrow)): 
                print( str(dd_arrow) + '   idempotents are messed up in this dd_arrow!1')
                count_of_mismatches+=1

            # matching out right algebra and in gen
            if not check_idempotents_match_right_right(dd_in_mod_gen(dd_arrow),dd_out_right_alg_gen(dd_arrow)): 
                print( str(dd_arrow) + '   idempotents are messed up in this dd_arrow!1')
                count_of_mismatches+=1


            # we don't check the following:
            # matchings between generators if differential is pure
        return count_of_mismatches==0

    def compute_d_squared(self):
        d_squared=Bunch_of_arrows()
        
        #contribution of differential on both of outgoing algebra elements
        for arrow in self.dd_arrows:
            # left
            a = dd_out_left_alg_gen(arrow)
            if a:
                for b in [algebra_diff_arrow[1] for algebra_diff_arrow in self.left_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                    ar=(dd_in_mod_gen(arrow),b,dd_out_mod_gen(arrow),dd_out_right_alg_gen(arrow))
                    d_squared[ar]+=1

            # right
            a = dd_out_right_alg_gen(arrow)
            if a:
                for b in [algebra_diff_arrow[1] for algebra_diff_arrow in self.right_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                    ar=(dd_in_mod_gen(arrow),dd_out_left_alg_gen(arrow),dd_out_mod_gen(arrow),b)
                    d_squared[ar]+=1

        #contribution of double dd_arrows
        for dd_arrow1 in self.dd_arrows:
            for dd_arrow2 in self.dd_arrows:
                # if str(dd_arrow1)=="(p01_c1, 1, p01, c1)" and str(dd_arrow2)=="(p01, p01, l1, 1)":
                #     debug('Yes')
                if not dd_out_mod_gen(dd_arrow1)==dd_in_mod_gen(dd_arrow2): continue
                a1a2_left=self.left_algebra.multiply(dd_out_left_alg_gen(dd_arrow1),dd_out_left_alg_gen(dd_arrow2))
                a1a2_right=self.right_algebra.multiply(dd_out_right_alg_gen(dd_arrow2),dd_out_right_alg_gen(dd_arrow1))
                if a1a2_left and a1a2_right:
                    ar=(dd_in_mod_gen(dd_arrow1),a1a2_left,dd_out_mod_gen(dd_arrow2),a1a2_right)
                    d_squared[ar]+=1
                    # if str(ar)=="(p01_c1, p01, l1, c1)": 
                    #     debug(str(dd_arrow1))
                    #     debug(str(dd_arrow2))
        return d_squared


# bijection = list of pairs (tuples)
# set of bijections= list of lists of pairs
# possible_images_of_generators = dictionary with keys from genset1, and values lists of gens in genset2

from .da_bimodule import find_bijections

def are_equal_smart_dd(DD1,DD2):
    print( 'number of generators in bimodules we compare are {} and {}'.format(str(len(DD1.genset)),str(len(DD2.genset))))
    print( 'number of actions in bimodules we compare are {} and {}'.format(str(len(DD1.dd_arrows)),str(len(DD2.dd_arrows))))
    
    if len(DD1.genset)!=len(DD2.genset): 
        print( 'incompatible sizes')
        return False
    
    # next we compute all the possible module homos w.r.t. idempotent ring
    possible_images_of_generators={}
    for gen1 in DD1.genset:
        possible_images_of_generators[gen1]=[]
        for gen2 in DD2.genset:
            if gen1.idem.left==gen2.idem.left and gen1.idem.right==gen2.idem.right: 
                possible_images_of_generators[gen1].append(gen2)
    
    # now, let us check if any of module homos satisfy df=0
    module_homos=find_bijections(DD1.genset,possible_images_of_generators,DD2.genset)
    for bijection in module_homos:
        f=Bunch_of_arrows()
        for pair in bijection:
            f[(pair[0],
                1,pair[1],1)]+=1
        if dd_check_df_is_0(DD1,DD2,f): return True
    print( 'there is no bijective isomorphism')
    return False


# morphism is represented by bunch of arrows, where all coefficients are 1!
# and also I assume that there are no differentials in algebra!
def dd_compute_df(DD1,DD2,f):
    df=Bunch_of_arrows()

    #contribution of double arrows, with morphism on the second place
    for arrow1 in DD1.dd_arrows:
        for arrow2 in [arrow2 for arrow2 in f if dd_out_mod_gen(arrow1)==dd_in_mod_gen(arrow2)]:
            left_a=DD1.left_algebra.multiply(dd_out_left_alg_gen(arrow1),dd_out_left_alg_gen(arrow2))
            right_a=DD1.right_algebra.multiply(dd_out_right_alg_gen(arrow2),dd_out_right_alg_gen(arrow1))
            if left_a and right_a:
                
                new_arrow_in_df=(dd_in_mod_gen(arrow1),
                    left_a,dd_out_mod_gen(arrow2),right_a)
                
                df[new_arrow_in_df]+=1

    #contribution of double arrows, with morphism on the first place
    for arrow1 in f:
        for arrow2 in [arrow2 for arrow2 in DD2.dd_arrows if dd_out_mod_gen(arrow1)==dd_in_mod_gen(arrow2)]:
            # if str(arrow1)='(x1, 1, x1⊗x1, 1)' and str(arrow2)=
            left_a=DD2.left_algebra.multiply(dd_out_left_alg_gen(arrow1),dd_out_left_alg_gen(arrow2))
            right_a=DD2.right_algebra.multiply(dd_out_right_alg_gen(arrow2),dd_out_right_alg_gen(arrow1))
            
            if left_a and right_a:
                new_arrow_in_df=(dd_in_mod_gen(arrow1),
                    left_a,dd_out_mod_gen(arrow2),right_a)


                df[new_arrow_in_df]+=1


    #contribution of differentials      
    for arrow in f:
        # contribution of differential on left algebra out element
        a = dd_out_left_alg_gen(arrow)
        if a:
            for b in [algebra_diff_arrow[1] for algebra_diff_arrow in DD1.left_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                new_arrow_in_df=(dd_in_mod_gen(arrow),
                            b,dd_out_mod_gen(arrow),dd_out_right_alg_gen(arrow) )
                df[new_arrow_in_df]+=1

        # contribution of differential on right algebra out element
        a = dd_out_right_alg_gen(arrow)
        if a:
            for b in [algebra_diff_arrow[1] for algebra_diff_arrow in DD2.right_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                new_arrow_in_df=(dd_in_mod_gen(arrow),
                    dd_out_left_alg_gen(arrow),dd_out_mod_gen(arrow),b)
                df[new_arrow_in_df]+=1

    return df

def dd_check_df_is_0(DD1,DD2,f):
    df=dd_compute_df(DD1,DD2,f)
    df.delete_arrows_with_even_coeff()
    # print( '\n{\nHere are all the arrows in df')
    # df.show()
    # print( '}')
    for arrow in df:
        if df[arrow] % 2 != 0:
            return False
    print( 'Everything is ok, we found isomorphism f such that df=0!')
    return True

def dd_cancel_pure_differential(DDbimodule_old,pure_differential):
    #z1--->z2
    z1=dd_in_mod_gen(pure_differential)
    z2=dd_out_mod_gen(pure_differential)

    #form generators
    generators_of_new_DD=DDbimodule_old.gen_by_name
    del generators_of_new_DD[z1.name]
    del generators_of_new_DD[z2.name]

    #form old differentials
    old_arrows_that_survive=[arrow for arrow in DDbimodule_old.dd_arrows if (dd_out_mod_gen(arrow)!=z1 and dd_out_mod_gen(arrow)!=z2  and dd_in_mod_gen(arrow)!=z1 and dd_in_mod_gen(arrow)!=z2 ) ]
    arrows_in_new_DD=Bunch_of_arrows(old_arrows_that_survive)

    #form new differentials
    arrows_in_z2=[arrow for arrow in DDbimodule_old.dd_arrows if (dd_out_mod_gen(arrow)==z2 and dd_in_mod_gen(arrow)!=z1 and dd_in_mod_gen(arrow)!=z2)]
    arrows_from_z1=[arrow for arrow in DDbimodule_old.dd_arrows if (dd_in_mod_gen(arrow)==z1 and dd_out_mod_gen(arrow)!=z2 and dd_out_mod_gen(arrow)!=z1)]

    for arrow_in_z2 in arrows_in_z2:
        for arrow_from_z1 in arrows_from_z1:

            new_out_left_alg_gen=DDbimodule_old.left_algebra.multiply(dd_out_left_alg_gen(arrow_in_z2),dd_out_left_alg_gen(arrow_from_z1))
            new_out_right_alg_gen=DDbimodule_old.right_algebra.multiply(dd_out_right_alg_gen(arrow_from_z1),dd_out_right_alg_gen(arrow_in_z2))

            if new_out_left_alg_gen and new_out_right_alg_gen:
                # debug(arrow_in_z2)
                # debug(arrow_from_z1)
                # print( '\n')
                new_arrow=(dd_in_mod_gen(arrow_in_z2),
                    new_out_left_alg_gen,dd_out_mod_gen(arrow_from_z1),new_out_right_alg_gen )
                arrows_in_new_DD[new_arrow]+=1


    arrows_in_new_DD.delete_arrows_with_even_coeff()

    return DD_bimodule(generators_of_new_DD,arrows_in_new_DD,DDbimodule_old.left_algebra,DDbimodule_old.right_algebra,name= DDbimodule_old.name +'_red')


def dd_randomly_cancel_until_possible(DD):
    there_is_diff=0
    arrs=list(DD.dd_arrows)
    shuffle(arrs)
    for arrow in arrs:
        if (dd_out_left_alg_gen(arrow)==1 and dd_out_right_alg_gen(arrow)==1 and dd_in_mod_gen(arrow)!=dd_out_mod_gen(arrow)):
            # now we need to check that there are no other arrows between these two generators
            other_arrows=[ar for ar in arrs if (dd_in_mod_gen(ar)==dd_in_mod_gen(arrow) and dd_out_mod_gen(ar)==dd_out_mod_gen(arrow) and ar!=arrow)]
            if other_arrows: continue
            else:
                there_is_diff=1
                canceled_DD=dd_cancel_pure_differential(DD,arrow)
                return (dd_randomly_cancel_until_possible(canceled_DD))

    if there_is_diff==0:
        return DD




