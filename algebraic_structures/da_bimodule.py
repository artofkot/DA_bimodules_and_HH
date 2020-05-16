# -*- coding: utf-8 -*- 
from .basics import *
from random import shuffle
from .basics import Bunch_of_arrows


def da_arrow_to_str(tuplee): #works only for DA arrows
    if len(tuplee)==4:
        return str(tuplee[0]) +'⊗'+ str(tuplee[1]) + "---->" + str(tuplee[2]) +'⊗'+ str(tuplee[3])  

def da_in_mod_gen(da_arrow):
    return da_arrow[0]

def da_in_alg_tuple(da_arrow):
    return da_arrow[1]

def da_out_alg_gen(da_arrow):
    return da_arrow[2]

def da_out_mod_gen(da_arrow):
    return da_arrow[3]

# conventions: D side is left, A side is right
class DA_bimodule(object):
    def __init__(self,gen_by_name,da_arrows,left_algebra,right_algebra,name,to_check=True,number_of_times_reduced=0):
        self.name=name
        self.number_of_times_reduced=number_of_times_reduced
        self.gen_by_name=gen_by_name
        self.genset=list(self.gen_by_name.values())
        #differentials are represented by bunch of arrows with coefficients 1
        self.da_arrows=da_arrows
        self.da_arrows.delete_arrows_with_even_coeff()
        self.left_algebra=left_algebra
        self.right_algebra=right_algebra

        if to_check==True:
            self.check()

    def check(self):
        ##Here we check that our " + self.name + " has all idempotents matching and d_squared=0:
        if not self.check_matching_of_idempotents_in_action():
            print ("\nSomething is wrong with idempotents!")
            raise NameError("DA_bimodule " + self.name + " has problems with idempotents")
        
        d_squared=self.compute_d_squared()
        d_squared.delete_arrows_with_even_coeff()

        if d_squared:
            print ("d_squared is not 0! the terms that are not canceled:")
            d_squared.show()
            raise NameError("DA_bimodule " + self.name + " doesn't satisfy d_squared=0.")


    def show(self):
        suffix=""
        if self.number_of_times_reduced>0:
            suffix="_canceled_"+str(self.number_of_times_reduced)+"_times"
        print ('Bimodule name: '+ self.name+suffix)
        print (str(len(self.genset)) + ' generators with their idempotents:')
        for gen in self.genset:
            print (str(gen.idem.left) + '___' + str(gen) + '___' + str(gen.idem.right) )

        print ('\n' + str(len(self.da_arrows))+ ' actions')
        for generator1 in self.genset:
            for generator2 in self.genset:
                arrows=[arrow for arrow in self.da_arrows if (da_in_mod_gen(arrow)==generator1 and da_out_mod_gen(arrow)==generator2)]
                if len(arrows)!=0:
                    print (generator1)
                    print ("↓", end=" ")
                    for ind, arrow in enumerate(arrows):
                        print(str(da_out_alg_gen(arrow)) + '⊗' + str(da_in_alg_tuple(arrow)),end=" ")
                        if ind+1!=len(arrows): print ('+', end=" ")
                        else: print(" ")
                    print (str(generator2) + '\n')
                
        # for arrow in self.da_arrows:
        #     print da_arrow_to_str(arrow)



    def show_for_tex(self):
        print ('=============')
        print ('Bimodule name: '+ self.name + ':')
        print ('\n\\vspace{0.2cm' +'}\n')
        print (str(len(self.genset)) + ' generators with their idempotents:')
        for gen in self.genset:
            # _{i_2}{(t_{12})}_{i_2}
            print ('$_{' + str(gen.idem.left) + '}{(' + str(gen) + ')}_{' + str(gen.idem.right) +'}$,')
        print ('\n\\vspace{0.2cm'+'}\n')

        print ('Actions:')
        for generator1 in self.genset:
            for generator2 in self.genset:
                arrows=[arrow for arrow in self.da_arrows if (da_in_mod_gen(arrow)==generator1 and da_out_mod_gen(arrow)==generator2)]
                if len(arrows)!=0:
                    for arrow in arrows:
                        print ('$'+str(da_in_mod_gen(arrow)) +'\\otimes ' + str(da_in_alg_tuple(arrow)) + '\\rightarrow ' + str(da_out_alg_gen(arrow)) +'\\otimes '+ str(da_out_mod_gen(arrow)) + '$, ')
                
                
        # for arrow in self.da_arrows:
        #     print da_arrow_to_str(arrow)

    def show_short(self):
        suffix=""
        if self.number_of_times_reduced>0:
            suffix="_canceled_"+str(self.number_of_times_reduced)+"_times"
        print ('Bimodule name: '+ self.name+suffix)
        print ('It has ' + str(len(self.genset))+ ' generators.')
        # for gen in self.genset:
        #     print (str(gen.idem.left) + '___' + str(gen) + '___' + str(gen.idem.right))
        print ('It has ' + str(len(self.da_arrows))+ ' actions.\n')


    def check_matching_of_idempotents_in_action(self): 
        count_of_mismatches=0
        for da_arrow in self.da_arrows:

            # matching out algebra and out gen
            if not check_idempotents_match_right_left( da_out_alg_gen(da_arrow),da_out_mod_gen(da_arrow) ): 
                print (da_arrow )
                count_of_mismatches+=1

            # matching out algebra and in gen
            if not check_idempotents_match_left_left( da_in_mod_gen(da_arrow),da_out_alg_gen(da_arrow) ): 
                print (da_arrow_to_str(da_arrow) + '   idempotents are messed up in this da_arrow!1')
                count_of_mismatches+=1

            # matching in first element of algebra and in gen
            if len(da_in_alg_tuple(da_arrow))!=0:
                if not check_idempotents_match_right_left( da_in_mod_gen(da_arrow),da_in_alg_tuple(da_arrow)[0] ): 
                    print (da_arrow_to_str(da_arrow) + '   idempotents are messed up in this da_arrow!2')
                    count_of_mismatches+=1

            if len(da_in_alg_tuple(da_arrow))!=0:
                if not check_idempotents_match_right_right( da_in_alg_tuple(da_arrow)[-1],da_out_mod_gen(da_arrow) ): 
                    print (da_arrow_to_str(da_arrow) + '   idempotents are messed up in this da_arrow!2')
                    count_of_mismatches+=1
            
            for i in range(len(da_in_alg_tuple(da_arrow))):
                if (i==0):
                    #matching in algebra and gen
                    if not check_idempotents_match_right_left(da_in_mod_gen(da_arrow),da_in_alg_tuple(da_arrow)[0]): 
                        print (da_arrow_to_str(da_arrow) + '   idempotents are messed up in this da_arrow!3')
                        count_of_mismatches+=1
                else:
                    #matching in algebras
                    if not check_idempotents_match_right_left(da_in_alg_tuple(da_arrow)[i-1],da_in_alg_tuple(da_arrow)[i]): 
                        print (da_arrow)
                        count_of_mismatches+=1

            # we don't check the following:
            # matchings between generators if tuple is empty
        return count_of_mismatches==0

    def compute_d_squared(self):
        d_squared=Bunch_of_arrows()
        #contribution of double arrows
        for arrow1 in self.da_arrows:
            for arrow2 in self.da_arrows:
                if not da_out_mod_gen(arrow1)==da_in_mod_gen(arrow2): continue
                a1a2=self.left_algebra.multiply(da_out_alg_gen(arrow1),da_out_alg_gen(arrow2))
                if a1a2:
                    ar=(da_in_mod_gen(arrow1), da_in_alg_tuple(arrow1) + da_in_alg_tuple(arrow2),a1a2,da_out_mod_gen(arrow2))
                    d_squared[ar]+=1


        for arrow in self.da_arrows:
            # contribution of differential on an algebra out element
            a = da_out_alg_gen(arrow)
            if a:
                for b in [algebra_diff_arrow[1] for algebra_diff_arrow in self.left_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                    ar=(da_in_mod_gen(arrow), da_in_alg_tuple(arrow),
                            b,da_out_mod_gen(arrow))
                    d_squared[ar]+=1


            # contribution of factorizing algebra elements in a tuple , and differentials       
            for index, a in enumerate(da_in_alg_tuple(arrow)):
                # factorization
                for factorization in getattr(a,'factorizations', []):
                    new_tuple=da_in_alg_tuple(arrow)[:index] + factorization + da_in_alg_tuple(arrow)[index+1:]
                    ar=(da_in_mod_gen(arrow), new_tuple,
                        da_out_alg_gen(arrow),da_out_mod_gen(arrow))
                    d_squared[ar]+=1

                # predifferentials on the right
                for b in [algebra_diff_arrow[0] for algebra_diff_arrow in self.right_algebra.algebra_diff_arrows if algebra_diff_arrow[1]==a]:
                    new_tuple=da_in_alg_tuple(arrow)[:index] + (b,) + da_in_alg_tuple(arrow)[index+1:]
                    ar=(da_in_mod_gen(arrow), new_tuple,
                        da_out_alg_gen(arrow),da_out_mod_gen(arrow))
                    d_squared[ar]+=1



        return d_squared

def cancel_pure_differential(DAbimodule_old,pure_differential):
    #z1--->z2
    z1=da_in_mod_gen(pure_differential)
    z2=da_out_mod_gen(pure_differential)

    #form generators
    generators_of_new_DA=DAbimodule_old.gen_by_name
    del generators_of_new_DA[z1.name]
    del generators_of_new_DA[z2.name]

    #form old differentials
    old_arrows_that_survive=[arrow for arrow in DAbimodule_old.da_arrows if (da_out_mod_gen(arrow)!=z1 and da_out_mod_gen(arrow)!=z2  and da_in_mod_gen(arrow)!=z1 and da_in_mod_gen(arrow)!=z2 ) ]
    arrows_in_new_DA=Bunch_of_arrows(old_arrows_that_survive)

    #form new differentials
    arrows_in_z2=[arrow for arrow in DAbimodule_old.da_arrows if (da_out_mod_gen(arrow)==z2 and da_in_mod_gen(arrow)!=z1 and da_in_mod_gen(arrow)!=z2)]
    arrows_from_z1=[arrow for arrow in DAbimodule_old.da_arrows if (da_in_mod_gen(arrow)==z1 and da_out_mod_gen(arrow)!=z2 and da_out_mod_gen(arrow)!=z1)]

    for arrow_in_z2 in arrows_in_z2:
        for arrow_from_z1 in arrows_from_z1:
            new_out_alg_gen=DAbimodule_old.left_algebra.multiply(da_out_alg_gen(arrow_in_z2),da_out_alg_gen(arrow_from_z1))

            if new_out_alg_gen:
                new_arrow=(da_in_mod_gen(arrow_in_z2), da_in_alg_tuple(arrow_in_z2)+da_in_alg_tuple(arrow_from_z1),
                    new_out_alg_gen,da_out_mod_gen(arrow_from_z1) )
                arrows_in_new_DA[new_arrow]+=1


    arrows_in_new_DA.delete_arrows_with_even_coeff()

    DADA=DA_bimodule(generators_of_new_DA,arrows_in_new_DA,DAbimodule_old.left_algebra,DAbimodule_old.right_algebra,name= DAbimodule_old.name,number_of_times_reduced= DAbimodule_old.number_of_times_reduced+1)
    # DADA.show_short()
    return DADA

def da_randomly_cancel_until_possible(DA1):
    there_is_diff=0
    # arrs=DA1.da_arrows
    DA1.da_arrows.delete_arrows_with_even_coeff()
    arrs=list(DA1.da_arrows)
    shuffle(arrs)
    for arrow in arrs:
        if (da_in_alg_tuple(arrow)==() and da_out_alg_gen(arrow)==1 and da_in_mod_gen(arrow)!=da_out_mod_gen(arrow)):
            
            # now we need to check that there are no other arrows between these two generators
            other_arrows=[ar for ar in arrs if (da_in_mod_gen(ar)==da_in_mod_gen(arrow) and da_out_mod_gen(ar)==da_out_mod_gen(arrow) and ar!=arrow)]
            if other_arrows: continue
            else:
                there_is_diff=1
                canceled_DA=cancel_pure_differential(DA1,arrow)
                return (da_randomly_cancel_until_possible(canceled_DA))

    if there_is_diff==0:
        return DA1

def da_cancel_until_possible(DA1):
    there_is_diff=0
    # arrs=DA1.da_arrows
    DA1.da_arrows.delete_arrows_with_even_coeff()
    arrs=sorted(list(DA1.da_arrows),key=str)
    for arrow in arrs:
        if (da_in_alg_tuple(arrow)==() and da_out_alg_gen(arrow)==1 and da_in_mod_gen(arrow)!=da_out_mod_gen(arrow)):
            
            # now we need to check that there are no other arrows between these two generators
            other_arrows=[ar for ar in arrs if (da_in_mod_gen(ar)==da_in_mod_gen(arrow) and da_out_mod_gen(ar)==da_out_mod_gen(arrow) and ar!=arrow)]
            if other_arrows: continue
            else:
                there_is_diff=1
                canceled_DA=cancel_pure_differential(DA1,arrow)
                return (da_cancel_until_possible(canceled_DA))

    if there_is_diff==0:
        return DA1

# bijection = list of pairs (tuples)
# set of bijections= list of lists of pairs
# possible_images_of_generators = dictionary with keys from genset1, and values lists of gens in genset2
def find_bijections(genset1,possible_images_of_generators,genset2):
    # print genset1

    if len(genset1)!=len(genset2): raise NameError("tried to find bijections between different size sets")
    if len(genset1)==1:
        gen1=genset1[0]
        gen2=genset2[0]
        if gen2 in possible_images_of_generators[gen1]: 
            all_bijections=[[(gen1,gen2)]]
        else: all_bijections=[]
    else:
        all_bijections=[]
        gen1=genset1[0]
        # if len(genset1)==10:
        #     print gen1
        for gen2 in possible_images_of_generators[gen1]:
            if not gen2 in genset2:continue
            # print possible_images_of_generators[gen1]
            # print str(gen1) +' '+ str(gen2)
            a1=list(genset1)
            a2=list(genset2)
            a1.remove(gen1)
            a2.remove(gen2)
            bijections_in_the_rest=find_bijections(a1,possible_images_of_generators,a2)
            bijections_with_fixed_first_pair=[bij+[(gen1,gen2)] for bij in bijections_in_the_rest]
            all_bijections=all_bijections+bijections_with_fixed_first_pair
    return all_bijections
    
def are_equal_smart_da(DA1,DA2,verbose=False):
    print ('Number of generators in bimodules we compare are {} and {}'.format(str(len(DA1.genset)),str(len(DA2.genset))))
    print ('Number of actions in bimodules we compare are {} and {}'.format(str(len(DA1.da_arrows)),str(len(DA2.da_arrows))))
    
    if len(DA1.genset)!=len(DA2.genset): return False
    
    # next we compute all the possible module homos w.r.t. idempotent ring
    possible_images_of_generators={}
    for gen1 in DA1.genset:
        possible_images_of_generators[gen1]=[]
        for gen2 in DA2.genset:
            if gen1.idem.left==gen2.idem.left and gen1.idem.right==gen2.idem.right: 
                possible_images_of_generators[gen1].append(gen2)
    # now, let us check any of module homos satisfy df=0
    
    # print find_bijections(['x','z'],{'x':['y','t'],'z':['y','t']},['y','t'])
    module_homos=find_bijections(DA1.genset,possible_images_of_generators,DA2.genset)
    
    for bijection in module_homos:
        f=Bunch_of_arrows()
        for pair in bijection:
            f[(pair[0],(),
                1,pair[1])]+=1
        if da_check_df_is_0(DA1,DA2,f): 
            if verbose: 
                print('The following are DA actions in the isomorphism f:')
                f.show()
            return True
    print ('there is no bijective isomorphism')
    return False


# morphism is represented by bunch of arrows, where all coefficients are 1!
# and also I assume tha there are no differentials in algebra!
def da_compute_df(DA1,DA2,f):
    df=Bunch_of_arrows()

    #contribution of double arrows, with morphism on the second place
    for arrow1 in DA1.da_arrows:
        for arrow2 in [arrow2 for arrow2 in f if da_out_mod_gen(arrow1)==da_in_mod_gen(arrow2)]:
        # for arrow2 in f:
        #     if not da_out_mod_gen(arrow1)==da_in_mod_gen(arrow2): continue
            a1a2=DA1.left_algebra.multiply(da_out_alg_gen(arrow1),da_out_alg_gen(arrow2))
            if a1a2:
                new_arrow_in_df=(da_in_mod_gen(arrow1), da_in_alg_tuple(arrow1) + da_in_alg_tuple(arrow2),a1a2,da_out_mod_gen(arrow2))
                df[new_arrow_in_df]+=1

    #contribution of double arrows, with morphism on the first place
    for arrow1 in f:
        for arrow2 in [arrow2 for arrow2 in DA2.da_arrows if da_out_mod_gen(arrow1)==da_in_mod_gen(arrow2)]:
        # for arrow2 in DA2.da_arrows:
        #     if not da_out_mod_gen(arrow1)==da_in_mod_gen(arrow2): continue
            a1a2=DA1.left_algebra.multiply(da_out_alg_gen(arrow1),da_out_alg_gen(arrow2))
            if a1a2:
                new_arrow_in_df=(da_in_mod_gen(arrow1), da_in_alg_tuple(arrow1) + da_in_alg_tuple(arrow2),
                            a1a2,da_out_mod_gen(arrow2))
                df[new_arrow_in_df]+=1


    #contribution of factorizing algebra elements, as well as predifferentiating       
    for arrow in f:
        # contribution of differential on an algebra out element
        a = da_out_alg_gen(arrow)
        if a:
            for b in [algebra_diff_arrow[1] for algebra_diff_arrow in DA1.left_algebra.algebra_diff_arrows if algebra_diff_arrow[0]==a]:
                new_arrow_in_df=(da_in_mod_gen(arrow), da_in_alg_tuple(arrow),
                            b,da_out_mod_gen(arrow))
                df[new_arrow_in_df]+=1


        for index, a in enumerate(da_in_alg_tuple(arrow)):
            # factorizing
            for factorization in getattr(a,'factorizations', []):
                new_tuple=da_in_alg_tuple(arrow)[:index] + factorization + da_in_alg_tuple(arrow)[index+1:]
                new_arrow_in_df=(da_in_mod_gen(arrow), new_tuple,
                    da_out_alg_gen(arrow),da_out_mod_gen(arrow))
                # add arrow to df
                df[new_arrow_in_df]+=1

            # predifferential
            for b in [algebra_diff_arrow[0] for algebra_diff_arrow in DA1.right_algebra.algebra_diff_arrows if algebra_diff_arrow[1]==a]:
                new_tuple=da_in_alg_tuple(arrow)[:index] + (b,) + da_in_alg_tuple(arrow)[index+1:]
                new_arrow_in_df=(da_in_mod_gen(arrow), new_tuple,
                        da_out_alg_gen(arrow),da_out_mod_gen(arrow))
                df[new_arrow_in_df]+=1
    return df

def da_check_df_is_0(DA1,DA2,f):
    df=da_compute_df(DA1,DA2,f)
    df.delete_arrows_with_even_coeff()
    # print '\n{\nHere are all the arrows in df'
    # df.show()
    # print '}'
    for arrow in df:
        if df[arrow] % 2 != 0:
            # print "Something is wrong with df=0 for this one, " + str(len(df)) + "arrows in df"
            return False
    print ('Everything is ok, we found isomorphism f such that df=0!')
    return True

def composition(f,g,A):
    f_g=Bunch_of_arrows()

    #contribution of double arrows
    for arrow1 in f:
        for arrow2 in g:
            if not da_out_mod_gen(arrow1)==da_in_mod_gen(arrow2): continue
            a1a2=A.multiply(da_out_alg_gen(arrow1),da_out_alg_gen(arrow2))
            if a1a2:
                new_arrow_in_df=(da_in_mod_gen(arrow1), da_in_alg_tuple(arrow1) + da_in_alg_tuple(arrow2),
                            a1a2,da_out_mod_gen(arrow2))
                f_g[new_arrow_in_df]+=1

    f_g.delete_arrows_with_even_coeff()
    return f_g


