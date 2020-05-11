# -*- coding: utf-8 -*- 
from algebraic_structures.da_bimodule import  Bunch_of_arrows
from algebraic_structures.da_bimodule  import da_in_mod_gen, da_out_mod_gen, da_in_alg_tuple, da_out_alg_gen
from .chain_complex import ChainComplex, homology_dim
from .basics import is_cyclic, AttrDict
# import time, sys, random


def is_bounded(DA1):
    graph={}
    for gen in DA1.genset:
        graph[gen]=()
        for arrow in DA1.da_arrows:
            if da_in_mod_gen(arrow)==gen:
                graph[gen]=graph[gen]+(da_out_mod_gen(arrow),)
    return not is_cyclic(graph)


def CH(DA1):
    def looping(container_of_final_diffs, current_gen, current_alg_tuple, first_gen, DA1,depth):
        # print 'depth=' + str(depth)
        # print 'first_gen=' + str(first_gen)
        # print 'current_gen=' + str(current_gen)
        # print 'current_alg_tuple' + str(current_alg_tuple)
        try:
            ending_arrows=[ar for ar in DA1.da_arrows if (da_in_mod_gen(ar)==current_gen and da_in_alg_tuple(ar)==current_alg_tuple and da_out_alg_gen(ar)==1  )  ]
            for ar in ending_arrows:
                container_of_final_diffs[(first_gen,da_out_mod_gen(ar))]+=1
        except: pass

        for ind in range(len(current_alg_tuple)+1):
            possible_next_arrows=[ar for ar in DA1.da_arrows if (da_in_mod_gen(ar)==current_gen and da_in_alg_tuple(ar)==current_alg_tuple[:ind] and da_out_alg_gen(ar)!=1 )]

            for arrow in possible_next_arrows:
                looping(container_of_final_diffs,
                        current_gen=da_out_mod_gen(arrow),
                        current_alg_tuple=current_alg_tuple[ind:]+(da_out_alg_gen(arrow),),
                        first_gen=first_gen,
                        DA1=DA1,depth=depth+1)


    # для убыстрения закоменть след две строки, если уверен что бимодули ограниченные.
    if not is_bounded(DA1):
        raise NameError("Bimodule " + DA1.name + " is not bounded, so can't take HH of it!")
    
    generators=[gen for gen in DA1.genset if gen.idem.left==gen.idem.right]
    generators_of_chain_complex_by_name=AttrDict({})
    for gen in generators:
        generators_of_chain_complex_by_name[gen.name]=gen


    differentials=Bunch_of_arrows()

    for first_gen in generators:
        first_arrows=[ar for ar in DA1.da_arrows if (da_in_mod_gen(ar)==first_gen and da_in_alg_tuple(ar)==() )]
        for first_arrow in first_arrows:
            if da_out_alg_gen(first_arrow)==1:
                # finding pure differentials
                differentials[(first_gen,da_out_mod_gen(first_arrow))]+=1
            else: 
                # finding all other differentials
                looping(container_of_final_diffs=differentials,
                    current_gen=da_out_mod_gen(first_arrow),
                    current_alg_tuple=(da_out_alg_gen(first_arrow),),
                    first_gen=first_gen,
                    DA1=DA1,depth=1)

    return ChainComplex(generators_of_chain_complex_by_name, differentials, "CH("+DA1.name+")")

def dimHH(DA1):
    C=CH(DA1)
    print("This is the Hochschild chain complex:")
    C.show_short()
    return homology_dim(C)








