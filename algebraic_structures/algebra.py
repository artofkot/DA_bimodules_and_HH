# -*- coding: utf-8 -*- 
from .basics import AttrDict, Bunch_of_arrows, Generator, debug
# from sets import Set

# I will assume that algebra is vector space over F2, and does not have differentials, 
# which happens in the first symmetric product
class dg_algebra(object):
    def __init__(self,gen_by_name,
                idem_by_name,
                multiplication_table,
                name,
                algebra_diff_arrows=Bunch_of_arrows([]),
                a_inf_actions=Bunch_of_arrows([]),
                to_check=True):
        self.algebra_diff_arrows=algebra_diff_arrows
        self.name=name
        self.gen_by_name=gen_by_name #dictionary of generators, where keys are their names
        self.genset=self.gen_by_name.values() #list of generators
        self.idem_by_name=idem_by_name #dictionary with idempotents, where keys are their names
        self.idemset=idem_by_name.values() #list of idempotents
        self.multiplication_table=multiplication_table
        self.a_inf_actions=a_inf_actions

        # adding factorizations of elements which are not idempotents into non-idempotent product
        for gen1 in set(self.genset) - set(self.idemset):
            for gen2 in set(self.genset) - set(self.idemset):
                gen3=self.multiplication_table.get((gen1,gen2),None)
                if gen3:
                    gen3.add_factorizations((gen1,gen2))
                    
        for gen in self.genset:
            for idem1 in self.idemset:
                for idem2 in self.idemset:
                    if (self.multiplication_table.get((gen,idem2))==gen) and (self.multiplication_table.get((idem1,gen))==gen):
                        gen.add_idems(idem1,idem2)
        if to_check==True:
            self.check()

    def check(self):
        # Here we check that our satisfies associativity, and Leibnitz rule (i.e. A_infty relations hold):
        
        #### Associativity we don't check, it takes too much time cause its hard to violate it:
        #### all my algebras are associative, and Bohua's algebras also
        # for gen1 in self.genset:
        #     for gen2 in self.genset:
        #         for gen3 in self.genset:
        #             if self.multiply(self.multiply(gen1,gen2), gen3)!=self.multiply(gen1, self.multiply(gen2,gen3)):
        #                 print "Associativity in" + self.name + " is broken"
        #                 raise NameError("Associativity in" + self.name + " is broken")

        # Leibnitz rule
        d_squared=Bunch_of_arrows([])
        for gen1 in self.genset:
            for gen2 in self.genset:
                gen3=self.multiply(gen1,gen2)
                if gen3: 
                    for a in self.algebra_diff_arrows:
                        if a[0]==gen3: d_squared[(gen1,gen2,a[1])]+=1

        for a in self.algebra_diff_arrows:
            for gen in self.genset:
                if self.multiply(gen,a[1]):  d_squared[(gen,a[0],self.multiply(gen,a[1]))]+=1
                if self.multiply(a[1],gen):  d_squared[(a[0],gen,self.multiply(a[1],gen))]+=1

        d_squared.delete_arrows_with_even_coeff()
        if d_squared:
            print("d_squared is not 0 in" + self.name + ". The terms that are not canceled:")
            d_squared.show()
            raise NameError("algebra " + self.name + " doesn't satisfy d_squared=0.")


    def mult2(self, gen1, gen2):
        if gen1==1: return gen2
        if gen2==1: return gen1
        if (gen1,gen2) in self.multiplication_table:
            return self.multiplication_table[(gen1,gen2)]
        else: return 0

    def multiply(self, *args):
        x=args[0]
        for i in range(len(args)-1):
            x=self.mult2(x,args[i+1])
        return x

    def show(self):
        print ('========\n'+ self.name)
        print('({} generators in this algebra)'.format(len(self.genset)) ) 
        for gen in self.genset:
            print(str(gen.idem.left) + ' ' + str(gen) + ' ' + str(gen.idem.right))
        print('Multiplications in algebra:')
        for k in self.multiplication_table:
            print(str(k[0]) + '*' + str(k[1]) + '=' + str(self.multiplication_table[k]))
