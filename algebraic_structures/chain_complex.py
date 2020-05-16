# -*- coding: utf-8 -*- 
from .basics import Bunch_of_arrows

class ChainComplex(object):
    def __init__(self,gen_by_name,arrows,name,to_check=True,number_of_times_reduced=0):
        self.name=name
        self.number_of_times_reduced=number_of_times_reduced
        self.gen_by_name=gen_by_name
        self.genset=self.gen_by_name.values()
        #differentials are represented by bunch of arrows with coefficients 1
        self.arrows=arrows
        self.arrows.delete_arrows_with_even_coeff()

        if to_check==True:
            dd=self.compute_dd()
            dd.delete_arrows_with_even_coeff()
            if dd:
                # self.show()
                # print '\ndd is wrong:'
                # dd.show()
                raise NameError("Chain complex " + self.name + " doesn't satisfy dd=0 !!!")

    def compute_dd(self):
        dd=Bunch_of_arrows()
        #contribution of double arrows
        for arrow1 in self.arrows:
            for arrow2 in self.arrows:
                if not arrow1[1]==arrow2[0]: continue
                else:  
                    dd[(arrow1[0],arrow2[1])]+=1
        return dd

    def show(self):
        suffix=""
        if self.number_of_times_reduced>0:
            suffix="_canceled_"+str(self.number_of_times_reduced)+"_times"
        print ('Chain complex name: '+self.name+suffix)
        print ("{} generators:".format(len(self.genset)))
        for gen in self.genset:
            print (str(gen) + " in grading " + str(gen.Z2grading))
        print ("\n{} differentials:".format(len(self.arrows)))
        self.arrows.show()

    def show_short(self):
        suffix=""
        if self.number_of_times_reduced>0:
            suffix="_canceled_"+str(self.number_of_times_reduced)+"_times"
        print ('Chain complex name: '+ self.name+suffix)
        print ('It has ' + str(len(self.genset))+ ' generators.')
        print ('It has ' + str(len(self.arrows))+ ' differentials.\n')
    
    def show_gradings(self):
        a=0
        b=0
        for gen in self.genset:
            if gen.Z2grading==0: a+=1
            else: b+=1
        print (str(a)+ ' generators in grading 0')
        print (str(b)+ ' generators in grading 1')

def cancel_differential(C,d):
    z1=d[0]
    z2=d[1]
    
    new_generators_by_name=C.gen_by_name.copy()
    del new_generators_by_name[z1.name]
    del new_generators_by_name[z2.name]
    
    old_arrows_that_survive=[arrow for arrow in C.arrows if (arrow[0]!=z1 and arrow[1]!=z2  and arrow[1]!=z1 and arrow[0]!=z2 ) ]
    new_arrows=Bunch_of_arrows(old_arrows_that_survive)
    
    arrows_in_z2=[arrow for arrow in C.arrows if (arrow[1]==z2 and arrow[0]!=z1 and arrow[0]!=z2)]
    arrows_from_z1=[arrow for arrow in C.arrows if (arrow[0]==z1 and arrow[1]!=z1 and arrow[1]!=z2)]

    for arrow_in_z2 in arrows_in_z2:
        for arrow_from_z1 in arrows_from_z1:
            new_arrows[(arrow_in_z2[0],arrow_from_z1[1])]+=1

    new_arrows.delete_arrows_with_even_coeff()
    C2=ChainComplex(new_generators_by_name,new_arrows, C.name, number_of_times_reduced=C.number_of_times_reduced+1)
    return C2

def homology_vector_space(C):
    C.arrows.delete_arrows_with_even_coeff()
    there_is_diff=0
    for arrow in C.arrows:
        if arrow[0]==arrow[1]: continue
        there_is_diff=1
        canceled_C=cancel_differential(C,arrow)
        return (homology_vector_space(canceled_C))

    if there_is_diff==0:
        # print ("\nGenerators of H(" + C.name + ") are:")
        # print (C.genset)
        # print ('\nHomology has dimension ' + str(len(C.genset)))
        C.name='Homology vector space after all cancellations'
        return C

def homology_dim(C):
    C.arrows.delete_arrows_with_even_coeff()
    there_is_diff=0
    for arrow in C.arrows:
        if arrow[0]==arrow[1]: continue
        there_is_diff=1
        canceled_C=cancel_differential(C,arrow)
        return (homology_dim(canceled_C))

    if there_is_diff==0:
        # print ("\nGenerators of H(" + C.name + ") are:")
        # print (C.genset)
        # print ('\nHomology has dimension ' + str(len(C.genset)))
        return len(C.genset)
        

