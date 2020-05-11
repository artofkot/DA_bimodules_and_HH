# -*- coding: utf-8 -*- 
import sys
# sys.path.append('../../')
sys.path.append('../')

from algebraic_structures.specific_algebras_modules_bimodules import torus_A

from input_DA_bimodules import (
    ID1,g2_M_RHD,g2_ID,g2_T_LHD) 

from algebraic_structures.basics import (
    AttrDict, Generator, Bunch_of_arrows,
    in_red)

from algebraic_structures.dd_bimodule import (
    DD_bimodule,are_equal_smart_dd)

from algebraic_structures.da_bimodule import (
    da_randomly_cancel_until_possible,
    are_equal_smart_da,
    DA_bimodule)

from algebraic_structures.aa_bimodule import (
    AA_bimodule, a_AA_a,
    aa_randomly_cancel_until_possible)

from algebraic_structures.dualization import (
    dualization_of_DD, dualization_of_AA)

from algebraic_structures.tensoring import (
 dd_aa_box_tensor_product,
 aa_da_box_tensor_product,
 da_dd_box_tensor_product,
 da_da_box_tensor_product,
 da_da_box_tensor_many_efficient_cancelations )

from .algebras_from_Bohua_program import (
    B1,A_pmc1_str3,
    B2,
    B3,
    B4,
    B5,
    B6,)

###### genus=2 

###### genus=2. First, for split pointed matched circle, which we denote by pmc1
def init_DD_pmc1(B1,A_pmc1_str3):
    gen_by_name=AttrDict({
            "x0": Generator("x0"),
            "x1": Generator("x1"),
            "x2": Generator("x2"),
            "x3": Generator("x3"),
            })

    gen_by_name.x0.add_idems(B1.idem_by_name.i0,A_pmc1_str3.idem_by_name['|(1, 3),(4, 6),(5, 7)|'] )
    gen_by_name.x1.add_idems(B1.idem_by_name.i1,A_pmc1_str3.idem_by_name['|(0, 2),(4, 6),(5, 7)|'])
    gen_by_name.x2.add_idems(B1.idem_by_name.i2,A_pmc1_str3.idem_by_name['|(0, 2),(1, 3),(5, 7)|'] )
    gen_by_name.x3.add_idems(B1.idem_by_name.i3,A_pmc1_str3.idem_by_name['|(0, 2),(1, 3),(4, 6)|'])
    


    dd_arrows=Bunch_of_arrows([
                    (gen_by_name.x0,
        B1.gen_by_name['r1'], 
        gen_by_name.x1, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(0->1)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r3'], 
        gen_by_name.x1, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(2->3)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r123'], 
        gen_by_name.x1, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(0->3)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r2'], 
        gen_by_name.x0, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(1->2)|']),

                    (gen_by_name.x2,
        B1.gen_by_name['r5'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(4->5)|']),

                    (gen_by_name.x2,
        B1.gen_by_name['r7'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(6->7)|']),

                    (gen_by_name.x2,
        B1.gen_by_name['r567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(4->7)|']),

                    (gen_by_name.x3,
        B1.gen_by_name['r6'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(5->6)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r1234'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(0->4)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r34'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(2->4)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r3456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(2->6)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r123456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(0->6)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r4'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(3->4)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(3->6)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r234'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(1->4)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r23456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(1->6)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r45'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(3->5)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r4567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(3->7)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r2345'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(1->5)|']),

                    (gen_by_name.x1,
        B1.gen_by_name['r234567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(1->7)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r345'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(2->5)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r34567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(2->7)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r12345'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(0->5)|']),

                    (gen_by_name.x0,
        B1.gen_by_name['r1234567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(0->7)|']),

        ])
    
    return DD_bimodule(gen_by_name,dd_arrows,B1,A_pmc1_str3,name="K_pmc1")
def init_DD_pmc1_M_RHD(A_pmc1_str1,A_pmc1_str3):
    gen_by_name=AttrDict({
            "x0": Generator("x0"),
            "x1": Generator("x1"),
            "x2": Generator("x2"),
            "x3": Generator("x3"),
            "r": Generator("r"),
            })

    gen_by_name.x0.add_idems(A_pmc1_str1.idem_by_name.i0,A_pmc1_str3.idem_by_name['|(1, 3),(4, 6),(5, 7)|'] )
    gen_by_name.x1.add_idems(A_pmc1_str1.idem_by_name.i1,A_pmc1_str3.idem_by_name['|(0, 2),(4, 6),(5, 7)|'])
    gen_by_name.x2.add_idems(A_pmc1_str1.idem_by_name.i2,A_pmc1_str3.idem_by_name['|(0, 2),(1, 3),(5, 7)|'] )
    gen_by_name.x3.add_idems(A_pmc1_str1.idem_by_name.i3,A_pmc1_str3.idem_by_name['|(0, 2),(1, 3),(4, 6)|'])
    gen_by_name.r.add_idems(A_pmc1_str1.idem_by_name.i1,A_pmc1_str3.idem_by_name['|(1, 3),(4, 6),(5, 7)|'])
        

    dd_arrows=Bunch_of_arrows([
                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r1'], 
        gen_by_name.x1, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(0->1)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r123'], 
        gen_by_name.x1, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(0->3)|']),

                    (gen_by_name.x2,
        A_pmc1_str1.gen_by_name['r5'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(4->5)|']),

                    (gen_by_name.x2,
        A_pmc1_str1.gen_by_name['r7'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(6->7)|']),

                    (gen_by_name.x2,
        A_pmc1_str1.gen_by_name['r567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(4->7)|']),

                    (gen_by_name.x3,
        A_pmc1_str1.gen_by_name['r6'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(1, 3),(5->6)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r1234'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(0->4)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r34'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(5, 7),(1->4),(2->3)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r3456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(5, 7),(1->6),(2->3)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r123456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(0->6)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r4'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(3->4)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(3->6)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r234'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(1->4)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r23456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(0, 2),(5, 7),(1->6)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r45'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(3->5)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r4567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(3->7)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r2345'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(1->5)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r234567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(0, 2),(4, 6),(1->7)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r345'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(4, 6),(1->5),(2->3)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r34567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(4, 6),(1->7),(2->3)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r12345'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(0->5)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r1234567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(0->7)|']),

                    (gen_by_name.r,
        A_pmc1_str1.gen_by_name['r2'], 
        gen_by_name.x0, 
        1),
                    (gen_by_name.r,
        1, 
        gen_by_name.x1, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(2->3)|']),

                    (gen_by_name.x1,
        A_pmc1_str1.gen_by_name['r23'], 
        gen_by_name.r, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(1->2)|']),

                    (gen_by_name.x0,
        A_pmc1_str1.gen_by_name['r3'], 
        gen_by_name.r, 
        A_pmc1_str3.gen_by_name['|(4, 6),(5, 7),(1->3)|']),

                    (gen_by_name.r,
        A_pmc1_str1.gen_by_name['r45'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(2->5)|']),

                    (gen_by_name.r,
        A_pmc1_str1.gen_by_name['r4567'], 
        gen_by_name.x3, 
        A_pmc1_str3.gen_by_name['|(1, 3),(4, 6),(2->7)|']),

                    (gen_by_name.r,
        A_pmc1_str1.gen_by_name['r4'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(2->4)|']),

                    (gen_by_name.r,
        A_pmc1_str1.gen_by_name['r456'], 
        gen_by_name.x2, 
        A_pmc1_str3.gen_by_name['|(1, 3),(5, 7),(2->6)|']),

        ])
    
    return DD_bimodule(gen_by_name,dd_arrows,A_pmc1_str1,A_pmc1_str3,name="DD_pmc1_M_RHD")
K_pmc1=init_DD_pmc1(B1,A_pmc1_str3)
DA_pmc1=g2_ID

####### experimenting with M_RHD
# g2_M_RHD=g2_M_RHD
# # g2_M_RHD.show()
# DD_pmc1_M_RHD=init_DD_pmc1_M_RHD(B1,A_pmc1_str3)
# mb_dd=da_dd_box_tensor_product(g2_M_RHD,K_pmc1)

# YES, THEY ARE EQUAL!


### DA1 arcslide
def init_DA1(B1,B2):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "r":Generator("r")
                })
    
    gen_by_name.x0.add_idems(B1.idem_by_name.i0,B2.idem_by_name['|(0, 2)|'])
    gen_by_name.x1.add_idems(B1.idem_by_name.i1,B2.idem_by_name['|(1, 6)|'])
    gen_by_name.x2.add_idems(B1.idem_by_name.i2,B2.idem_by_name['|(3, 5)|'])
    gen_by_name.x3.add_idems(B1.idem_by_name.i3,B2.idem_by_name['|(4, 7)|'])
    gen_by_name.r.add_idems(B1.idem_by_name.i2,B2.idem_by_name['|(1, 6)|'])

    arrows=Bunch_of_arrows([
        # short-near-chords:
        (              gen_by_name.x0,(B2.gen_by_name['|(0->1)|'],),
                B1.gen_by_name.r1,gen_by_name.x1),

        (              gen_by_name.x1,(B2.gen_by_name['|(1->2)|'],),
                B1.gen_by_name.r2,gen_by_name.x0),

        (              gen_by_name.x1,(),
                B1.gen_by_name.r4,gen_by_name.r),

        (              gen_by_name.x2,(B2.gen_by_name['|(3->4)|'],),
                B1.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(B2.gen_by_name['|(4->5)|'],),
                B1.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(B2.gen_by_name['|(5->6)|'],),
                1,gen_by_name.r),

        (              gen_by_name.r,(B2.gen_by_name['|(6->7)|'],),
                B1.gen_by_name.r7,gen_by_name.x3),

        # the rest, kind of follows from short ones, but not exactly

        (              gen_by_name.x2,(B2.gen_by_name['|(5->7)|'],),
                B1.gen_by_name.r7,gen_by_name.x3),

        (              gen_by_name.x2,(B2.gen_by_name['|(3->5)|'],),
                B1.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(B2.gen_by_name['|(4->7)|'],),
                B1.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(B2.gen_by_name['|(3->7)|'],),
                B1.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x0,(B2.gen_by_name['|(0->2)|'],),
                B1.gen_by_name.r12,gen_by_name.x0),

        #### next one might be basic choice
        (              gen_by_name.x3,(B2.gen_by_name['|(4->6)|'],),
                B1.gen_by_name.r6,gen_by_name.r),

        (              gen_by_name.x2,(B2.gen_by_name['|(3->6)|'],),
                B1.gen_by_name.r56,gen_by_name.r),

        (              gen_by_name.x0,(B2.gen_by_name['|(2->3)|'],B2.gen_by_name['|(5->6)|']),
                B1.gen_by_name.r3,gen_by_name.x1),

        (              gen_by_name.x0,(B2.gen_by_name['|(2->3)|'],),
                B1.gen_by_name.r34,gen_by_name.x2),

        (              gen_by_name.x0,(B2.gen_by_name['|(0->3)|'],B2.gen_by_name['|(5->6)|']),
                B1.gen_by_name.r123,gen_by_name.x1),

        (              gen_by_name.x1,(B2.gen_by_name['|(1->3)|'],B2.gen_by_name['|(5->6)|']),
                B1.gen_by_name.r23,gen_by_name.x1),

        (              gen_by_name.x1,(B2.gen_by_name['|(1->3)|'],),
                B1.gen_by_name.r234,gen_by_name.x2),

        (              gen_by_name.x0,(B2.gen_by_name['|(2->4)|'],),
                B1.gen_by_name.r345,gen_by_name.x3),

        (              gen_by_name.x0,(B2.gen_by_name['|(2->7)|'],),
                B1.gen_by_name.r34567,gen_by_name.x3),

        (              gen_by_name.x0,(B2.gen_by_name['|(2->6)|'],),
                B1.gen_by_name.r3456,gen_by_name.r),

        (              gen_by_name.x0,(B2.gen_by_name['|(2->5)|'],),
                B1.gen_by_name.r3456,gen_by_name.x2),

        (              gen_by_name.x1,(B2.gen_by_name['|(1->4)|'],),
                B1.gen_by_name.r2345,gen_by_name.x3),

        (              gen_by_name.x1,(B2.gen_by_name['|(1->7)|'],),
                B1.gen_by_name.r234567,gen_by_name.x3),

        (              gen_by_name.x0,(B2.gen_by_name['|(0->7)|'],),
                B1.gen_by_name.r1234567,gen_by_name.x3),

        (              gen_by_name.x0,(B2.gen_by_name['|(0->3)|'],),
                B1.gen_by_name.r1234,gen_by_name.x2),

        (              gen_by_name.x0,(B2.gen_by_name['|(0->5)|'],),
                B1.gen_by_name.r123456,gen_by_name.x2),

        (              gen_by_name.x0,(B2.gen_by_name['|(0->4)|'],),
                B1.gen_by_name.r12345,gen_by_name.x3),

        (              gen_by_name.x0,(B2.gen_by_name['|(0->6)|'],),
                B1.gen_by_name.r123456,gen_by_name.r),

        (              gen_by_name.x1,(B2.gen_by_name['|(1->6)|'],),
                B1.gen_by_name.r23456,gen_by_name.r),

        (              gen_by_name.x1,(B2.gen_by_name['|(1->5)|'],),
                B1.gen_by_name.r23456,gen_by_name.x2),

        ])

    return DA_bimodule(gen_by_name,arrows,B1,B2,name="N(η)")
DA1=init_DA1(B1,B2)
# DA1.show_for_tex()

### DA2 arcslide, over B2, B3
def init_DA2(B2,B3):
    gen_by_name=AttrDict({
                "y0": Generator("y0"),
                "y1": Generator("y1"),
                "y2": Generator("y2"),
                "y3": Generator("y3"),
                "r":Generator("r")
                })
    
    gen_by_name.y0.add_idems(B2.idem_by_name['|(3, 5)|'],B3.idem_by_name['|(2, 4)|'])
    gen_by_name.y1.add_idems(B2.idem_by_name['|(0, 2)|'],B3.idem_by_name['|(0, 3)|'])
    gen_by_name.y2.add_idems(B2.idem_by_name['|(1, 6)|'],B3.idem_by_name['|(1, 6)|'])
    gen_by_name.y3.add_idems(B2.idem_by_name['|(4, 7)|'],B3.idem_by_name['|(5, 7)|'])
    gen_by_name.r.add_idems(B2.idem_by_name['|(1, 6)|'],B3.idem_by_name['|(2, 4)|'])

    arrows=Bunch_of_arrows([
        # short-near-chords:
        (              gen_by_name.y1,(B3.gen_by_name['|(0->1)|'],),
            B2.gen_by_name['|(0->1)|'],gen_by_name.y2),
        (              gen_by_name.y2,(B3.gen_by_name['|(1->2)|'],),
            1,gen_by_name.r),
        (              gen_by_name.r,(B3.gen_by_name['|(2->3)|'],),
            B2.gen_by_name['|(1->2)|'],gen_by_name.y1),
        (              gen_by_name.y1,(B3.gen_by_name['|(3->4)|'],),
            B2.gen_by_name['|(2->3)|'],gen_by_name.y0),
        (              gen_by_name.y0,(B3.gen_by_name['|(4->5)|'],),
            B2.gen_by_name['|(3->4)|'],gen_by_name.y3),
        (              gen_by_name.y0,(),
            B2.gen_by_name['|(5->6)|'],gen_by_name.r),
        (              gen_by_name.y2,(B3.gen_by_name['|(6->7)|'],),
            B2.gen_by_name['|(6->7)|'],gen_by_name.y3),

        # the rest, kind of follows from short ones, but not exactly
        # the game changer for DA bimodule is usually periodic domain, which
        # contains σ and the domain next to it. See below

        (              gen_by_name.y1,(B3.gen_by_name['|(3->5)|'],),
            B2.gen_by_name['|(2->4)|'],gen_by_name.y3),

        (              gen_by_name.y2,(B3.gen_by_name['|(1->3)|'],),
            B2.gen_by_name['|(1->2)|'],gen_by_name.y1),

        (              gen_by_name.r,(B3.gen_by_name['|(2->4)|'],),
            B2.gen_by_name['|(1->3)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B3.gen_by_name['|(0->2)|'],),
            B2.gen_by_name['|(0->1)|'],gen_by_name.r),

        (              gen_by_name.y1,(B3.gen_by_name['|(0->3)|'],),
            B2.gen_by_name['|(0->2)|'],gen_by_name.y1),

        (              gen_by_name.r,(B3.gen_by_name['|(2->5)|'],),
            B2.gen_by_name['|(1->4)|'],gen_by_name.y3),

        (              gen_by_name.y2,(B3.gen_by_name['|(1->4)|'],),
            B2.gen_by_name['|(1->3)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B3.gen_by_name['|(0->4)|'],),
            B2.gen_by_name['|(0->3)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B3.gen_by_name['|(1->5)|'],),
            B2.gen_by_name['|(1->4)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B3.gen_by_name['|(0->5)|'],),
            B2.gen_by_name['|(0->4)|'],gen_by_name.y3),
        
        # next one is game changer 
        (              gen_by_name.y3,(B3.gen_by_name['|(5->6)|'],),
            B2.gen_by_name['|(4->6)|'],gen_by_name.y2),

        (              gen_by_name.y3,(B3.gen_by_name['|(5->7)|'],),
            B2.gen_by_name['|(4->7)|'],gen_by_name.y3),

        (              gen_by_name.y2,(B3.gen_by_name['|(1->6)|'],),
            B2.gen_by_name['|(1->6)|'],gen_by_name.y2),

        (              gen_by_name.y1,(B3.gen_by_name['|(3->6)|'],),
            B2.gen_by_name['|(2->6)|'],gen_by_name.y2),

        (              gen_by_name.y0,(B3.gen_by_name['|(4->6)|'],),
            B2.gen_by_name['|(3->6)|'],gen_by_name.y2),

        (              gen_by_name.r,(B3.gen_by_name['|(2->6)|'],),
            B2.gen_by_name['|(1->6)|'],gen_by_name.y2),

        (              gen_by_name.y0,(B3.gen_by_name['|(4->6)|'],B3.gen_by_name['|(1->2)|']),
            B2.gen_by_name['|(3->5)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B3.gen_by_name['|(5->6)|'],B3.gen_by_name['|(1->2)|']),
            B2.gen_by_name['|(4->5)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B3.gen_by_name['|(1->7)|'],),
            B2.gen_by_name['|(1->7)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B3.gen_by_name['|(0->7)|'],),
            B2.gen_by_name['|(0->7)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B3.gen_by_name['|(0->6)|'],),
            B2.gen_by_name['|(0->6)|'],gen_by_name.y2),

        (              gen_by_name.y1,(B3.gen_by_name['|(3->7)|'],),
            B2.gen_by_name['|(2->7)|'],gen_by_name.y3),

        (              gen_by_name.y0,(B3.gen_by_name['|(4->7)|'],),
            B2.gen_by_name['|(3->7)|'],gen_by_name.y3),

        (              gen_by_name.r,(B3.gen_by_name['|(2->7)|'],),
            B2.gen_by_name['|(1->7)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B3.gen_by_name['|(3->6)|'],B3.gen_by_name['|(1->2)|']),
            B2.gen_by_name['|(2->5)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B3.gen_by_name['|(0->6)|'],B3.gen_by_name['|(1->2)|']),
            B2.gen_by_name['|(0->5)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B3.gen_by_name['|(1->6)|'],B3.gen_by_name['|(1->2)|']),
            B2.gen_by_name['|(1->5)|'],gen_by_name.y0),

        (              gen_by_name.r,(B3.gen_by_name['|(2->6)|'],B3.gen_by_name['|(1->2)|']),
            B2.gen_by_name['|(1->5)|'],gen_by_name.y0),

        ])

    return DA_bimodule(gen_by_name,arrows,B2,B3,name="N(μ_1)")
DA2=init_DA2(B2,B3)
# DA2.show_for_tex()

### DA3 arcslide, over B3, B4
def init_DA3(B3,B4):
    gen_by_name=AttrDict({
                "y0": Generator("y0"),
                "y1": Generator("y1"),
                "y2": Generator("y2"),
                "y3": Generator("y3"),
                "r":Generator("r")
                })
    
    gen_by_name.y0.add_idems(B3.idem_by_name['|(5, 7)|'],B4.idem_by_name['|(2, 7)|'])
    gen_by_name.y1.add_idems(B3.idem_by_name['|(1, 6)|'],B4.idem_by_name['|(1, 6)|'])
    gen_by_name.y2.add_idems(B3.idem_by_name['|(2, 4)|'],B4.idem_by_name['|(3, 5)|'])
    gen_by_name.y3.add_idems(B3.idem_by_name['|(0, 3)|'],B4.idem_by_name['|(0, 4)|'])
    gen_by_name.r.add_idems(B3.idem_by_name['|(1, 6)|'],B4.idem_by_name['|(2, 7)|'])

    arrows=Bunch_of_arrows([
        # short-near-chords:
        (              gen_by_name.y3,(B4.gen_by_name['|(0->1)|'],),
            B3.gen_by_name['|(0->1)|'],gen_by_name.y1),

        (              gen_by_name.y1,(B4.gen_by_name['|(1->2)|'],),
            1,gen_by_name.r),

        (              gen_by_name.r,(B4.gen_by_name['|(2->3)|'],),
            B3.gen_by_name['|(1->2)|'],gen_by_name.y2),

        (              gen_by_name.y2,(B4.gen_by_name['|(3->4)|'],),
            B3.gen_by_name['|(2->3)|'],gen_by_name.y3),

        (              gen_by_name.y3,(B4.gen_by_name['|(4->5)|'],),
            B3.gen_by_name['|(3->4)|'],gen_by_name.y2),

        (              gen_by_name.y0,(),
            B3.gen_by_name['|(5->6)|'],gen_by_name.r),

        (              gen_by_name.y1,(B4.gen_by_name['|(6->7)|'],),
            B3.gen_by_name['|(6->7)|'],gen_by_name.y0),


        # the rest, kind of follows from short ones, but not exactly
        # the game changer for DA bimodule is usually periodic domain, which
        # contains σ and the domain next to it. See below

        (              gen_by_name.y2,(B4.gen_by_name['|(3->5)|'],),
            B3.gen_by_name['|(2->4)|'],gen_by_name.y2),

        (              gen_by_name.r,(B4.gen_by_name['|(2->4)|'],),
            B3.gen_by_name['|(1->3)|'],gen_by_name.y3),

        (              gen_by_name.y3,(B4.gen_by_name['|(0->2)|'],),
            B3.gen_by_name['|(0->1)|'],gen_by_name.r),

        (              gen_by_name.y1,(B4.gen_by_name['|(1->3)|'],),
            B3.gen_by_name['|(1->2)|'],gen_by_name.y2),

        (              gen_by_name.y1,(B4.gen_by_name['|(1->4)|'],),
            B3.gen_by_name['|(1->3)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B4.gen_by_name['|(1->5)|'],),
            B3.gen_by_name['|(1->4)|'],gen_by_name.y2),

        (              gen_by_name.r,(B4.gen_by_name['|(2->5)|'],),
            B3.gen_by_name['|(1->4)|'],gen_by_name.y2),

        (              gen_by_name.y3,(B4.gen_by_name['|(0->4)|'],),
            B3.gen_by_name['|(0->3)|'],gen_by_name.y3),

        (              gen_by_name.y3,(B4.gen_by_name['|(0->3)|'],),
            B3.gen_by_name['|(0->2)|'],gen_by_name.y2),

        (              gen_by_name.y3,(B4.gen_by_name['|(0->5)|'],),
            B3.gen_by_name['|(0->4)|'],gen_by_name.y2),

        # game changer
        (              gen_by_name.y2,(B4.gen_by_name['|(5->6)|'],),
            B3.gen_by_name['|(4->6)|'],gen_by_name.y1),

        (              gen_by_name.y2,(B4.gen_by_name['|(5->6)|'],B4.gen_by_name['|(1->2)|']),
            B3.gen_by_name['|(4->5)|'],gen_by_name.y0),

        (              gen_by_name.r,(B4.gen_by_name['|(2->6)|'],B4.gen_by_name['|(1->2)|']),
            B3.gen_by_name['|(1->5)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B4.gen_by_name['|(1->6)|'],),
            B3.gen_by_name['|(1->6)|'],gen_by_name.y1),

        (              gen_by_name.y1,(B4.gen_by_name['|(1->6)|'],B4.gen_by_name['|(1->2)|']),
            B3.gen_by_name['|(1->5)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B4.gen_by_name['|(3->6)|'],B4.gen_by_name['|(1->2)|']),
            B3.gen_by_name['|(2->5)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B4.gen_by_name['|(0->6)|'],B4.gen_by_name['|(1->2)|']),
            B3.gen_by_name['|(0->5)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B4.gen_by_name['|(4->6)|'],),
            B3.gen_by_name['|(3->6)|'],gen_by_name.y1),

        (              gen_by_name.r,(B4.gen_by_name['|(2->6)|'],),
            B3.gen_by_name['|(1->6)|'],gen_by_name.y1),

        (              gen_by_name.y2,(B4.gen_by_name['|(5->7)|'],),
            B3.gen_by_name['|(4->7)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B4.gen_by_name['|(3->6)|'],),
            B3.gen_by_name['|(2->6)|'],gen_by_name.y1),

        (              gen_by_name.y3,(B4.gen_by_name['|(0->6)|'],),
            B3.gen_by_name['|(0->6)|'],gen_by_name.y1),

        (              gen_by_name.y3,(B4.gen_by_name['|(4->6)|'],B4.gen_by_name['|(1->2)|']),
            B3.gen_by_name['|(3->5)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B4.gen_by_name['|(0->7)|'],),
            B3.gen_by_name['|(0->7)|'],gen_by_name.y0),

         (              gen_by_name.y3,(B4.gen_by_name['|(4->7)|'],),
            B3.gen_by_name['|(3->7)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B4.gen_by_name['|(1->7)|'],),
            B3.gen_by_name['|(1->7)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B4.gen_by_name['|(3->7)|'],),
            B3.gen_by_name['|(2->7)|'],gen_by_name.y0),

        (              gen_by_name.r,(B4.gen_by_name['|(2->7)|'],),
            B3.gen_by_name['|(1->7)|'],gen_by_name.y0),

        ])

    return DA_bimodule(gen_by_name,arrows,B3,B4,name="N(μ_2)")
DA3=init_DA3(B3,B4)
# DA3.show_for_tex()

### DA4 arcslide, over B4, B5
def init_DA4(B4,B5):
    gen_by_name=AttrDict({
                "y0": Generator("y0"),
                "y1": Generator("y1"),
                "y2": Generator("y2"),
                "y3": Generator("y3"),
                "r":Generator("r")
                })
    
    gen_by_name.y0.add_idems(B4.idem_by_name['|(3, 5)|'],B5.idem_by_name['|(2, 4)|'])
    gen_by_name.y1.add_idems(B4.idem_by_name['|(1, 6)|'],B5.idem_by_name['|(1, 6)|'])
    gen_by_name.y2.add_idems(B4.idem_by_name['|(2, 7)|'],B5.idem_by_name['|(3, 7)|'])
    gen_by_name.y3.add_idems(B4.idem_by_name['|(0, 4)|'],B5.idem_by_name['|(0, 5)|'])
    gen_by_name.r.add_idems(B4.idem_by_name['|(1, 6)|'],B5.idem_by_name['|(2, 4)|'])

    arrows=Bunch_of_arrows([
        # short-near-chords:
        (              gen_by_name.y3,(B5.gen_by_name['|(0->1)|'],),
            B4.gen_by_name['|(0->1)|'],gen_by_name.y1),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->2)|'],),
            1,gen_by_name.r),

        (              gen_by_name.r,(B5.gen_by_name['|(2->3)|'],),
            B4.gen_by_name['|(1->2)|'],gen_by_name.y2),

        (              gen_by_name.y2,(B5.gen_by_name['|(3->4)|'],),
            B4.gen_by_name['|(2->3)|'],gen_by_name.y0),

        (              gen_by_name.y0,(B5.gen_by_name['|(4->5)|'],),
            B4.gen_by_name['|(3->4)|'],gen_by_name.y3),

        (              gen_by_name.y0,(),
            B4.gen_by_name['|(5->6)|'],gen_by_name.r),

        (              gen_by_name.y1,(B5.gen_by_name['|(6->7)|'],),
            B4.gen_by_name['|(6->7)|'],gen_by_name.y2),

        # the rest, kind of follows from short ones, but not exactly
        # the game changer for DA bimodule is usually periodic domain, which
        # contains σ and the domain next to it. See below

        (              gen_by_name.r,(B5.gen_by_name['|(2->4)|'],),
            B4.gen_by_name['|(1->3)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B5.gen_by_name['|(0->2)|'],),
            B4.gen_by_name['|(0->1)|'],gen_by_name.r),

        (              gen_by_name.y2,(B5.gen_by_name['|(3->5)|'],),
            B4.gen_by_name['|(2->4)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->3)|'],),
            B4.gen_by_name['|(1->2)|'],gen_by_name.y2),

        (              gen_by_name.y3,(B5.gen_by_name['|(0->4)|'],),
            B4.gen_by_name['|(0->3)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->4)|'],),
            B4.gen_by_name['|(1->3)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B5.gen_by_name['|(0->3)|'],),
            B4.gen_by_name['|(0->2)|'],gen_by_name.y2),

        (              gen_by_name.r,(B5.gen_by_name['|(2->5)|'],),
            B4.gen_by_name['|(1->4)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->4)|'],),
            B4.gen_by_name['|(1->3)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->5)|'],),
            B4.gen_by_name['|(1->4)|'],gen_by_name.y3),

        (              gen_by_name.y3,(B5.gen_by_name['|(0->5)|'],),
            B4.gen_by_name['|(0->4)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->4)|'],),
            B4.gen_by_name['|(1->3)|'],gen_by_name.y0),


        # game changer

        (              gen_by_name.y3,(B5.gen_by_name['|(5->6)|'],),
            B4.gen_by_name['|(4->6)|'],gen_by_name.y1),

        (              gen_by_name.y3,(B5.gen_by_name['|(5->7)|'],),
            B4.gen_by_name['|(4->7)|'],gen_by_name.y2),

        (              gen_by_name.y2,(B5.gen_by_name['|(3->6)|'],),
            B4.gen_by_name['|(2->6)|'],gen_by_name.y1),

        (              gen_by_name.y0,(B5.gen_by_name['|(4->6)|'],),
            B4.gen_by_name['|(3->6)|'],gen_by_name.y1),

        (              gen_by_name.y3,(B5.gen_by_name['|(0->6)|'],),
            B4.gen_by_name['|(0->6)|'],gen_by_name.y1),

        (              gen_by_name.r,(B5.gen_by_name['|(2->6)|'],),
            B4.gen_by_name['|(1->6)|'],gen_by_name.y1),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->6)|'],),
            B4.gen_by_name['|(1->6)|'],gen_by_name.y1),

        (              gen_by_name.y3,(B5.gen_by_name['|(5->6)|'],B5.gen_by_name['|(1->2)|']),
            B4.gen_by_name['|(4->5)|'],gen_by_name.y0),

        (              gen_by_name.y0,(B5.gen_by_name['|(4->7)|'],),
            B4.gen_by_name['|(3->7)|'],gen_by_name.y2),

        (              gen_by_name.y2,(B5.gen_by_name['|(3->7)|'],),
            B4.gen_by_name['|(2->7)|'],gen_by_name.y2),

        (              gen_by_name.y3,(B5.gen_by_name['|(0->7)|'],),
            B4.gen_by_name['|(0->7)|'],gen_by_name.y2),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->7)|'],),
            B4.gen_by_name['|(1->7)|'],gen_by_name.y2),

        (              gen_by_name.r,(B5.gen_by_name['|(2->7)|'],),
            B4.gen_by_name['|(1->7)|'],gen_by_name.y2),

        (              gen_by_name.y0,(B5.gen_by_name['|(4->6)|'],B5.gen_by_name['|(1->2)|']),
            B4.gen_by_name['|(3->5)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B5.gen_by_name['|(1->6)|'],B5.gen_by_name['|(1->2)|']),
            B4.gen_by_name['|(1->5)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B5.gen_by_name['|(3->6)|'],B5.gen_by_name['|(1->2)|']),
            B4.gen_by_name['|(2->5)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B5.gen_by_name['|(0->6)|'],B5.gen_by_name['|(1->2)|']),
            B4.gen_by_name['|(0->5)|'],gen_by_name.y0),

        (              gen_by_name.r,(B5.gen_by_name['|(2->6)|'],B5.gen_by_name['|(1->2)|']),
            B4.gen_by_name['|(1->5)|'],gen_by_name.y0),






        ])

    return DA_bimodule(gen_by_name,arrows,B4,B5,name="N(μ_3)")
DA4=init_DA4(B4,B5)
# DA4.show_for_tex()

### DA5 arcslide, over B5, B6=B2
def init_DA5(B5,B6):
    gen_by_name=AttrDict({
                "y0": Generator("y0"),
                "y1": Generator("y1"),
                "y2": Generator("y2"),
                "y3": Generator("y3"),
                "r":Generator("r")
                })
    
    gen_by_name.y0.add_idems(B5.idem_by_name['|(0, 5)|'],B6.idem_by_name['|(0, 2)|'])
    gen_by_name.y1.add_idems(B5.idem_by_name['|(1, 6)|'],B6.idem_by_name['|(1, 6)|'])
    gen_by_name.y2.add_idems(B5.idem_by_name['|(2, 4)|'],B6.idem_by_name['|(3, 5)|'])
    gen_by_name.y3.add_idems(B5.idem_by_name['|(3, 7)|'],B6.idem_by_name['|(4, 7)|'])
    gen_by_name.r.add_idems(B5.idem_by_name['|(1, 6)|'],B6.idem_by_name['|(0, 2)|'])

    arrows=Bunch_of_arrows([
        # short-near-chords:
        (              gen_by_name.y0,(B6.gen_by_name['|(0->1)|'],),
            B5.gen_by_name['|(0->1)|'],gen_by_name.y1),

        (              gen_by_name.y1,(B6.gen_by_name['|(1->2)|'],),
            1,gen_by_name.r),

        (              gen_by_name.r,(B6.gen_by_name['|(2->3)|'],),
            B5.gen_by_name['|(1->2)|'],gen_by_name.y2),

        (              gen_by_name.y2,(B6.gen_by_name['|(3->4)|'],),
            B5.gen_by_name['|(2->3)|'],gen_by_name.y3),

        (              gen_by_name.y3,(B6.gen_by_name['|(4->5)|'],),
            B5.gen_by_name['|(3->4)|'],gen_by_name.y2),

        (              gen_by_name.y0,(),
            B5.gen_by_name['|(5->6)|'],gen_by_name.r),

        (              gen_by_name.y1,(B6.gen_by_name['|(6->7)|'],),
            B5.gen_by_name['|(6->7)|'],gen_by_name.y3),


        # the rest, kind of follows from short ones, but not exactly
        # the game changer for DA bimodule is usually periodic domain, which
        # contains σ and the domain next to it. See below

        (              gen_by_name.y1,(B6.gen_by_name['|(1->3)|'],),
            B5.gen_by_name['|(1->2)|'],gen_by_name.y2),

        (              gen_by_name.r,(B6.gen_by_name['|(2->4)|'],),
            B5.gen_by_name['|(1->3)|'],gen_by_name.y3),

        (              gen_by_name.y2,(B6.gen_by_name['|(3->5)|'],),
            B5.gen_by_name['|(2->4)|'],gen_by_name.y2),

        (              gen_by_name.y0,(B6.gen_by_name['|(0->2)|'],),
            B5.gen_by_name['|(0->1)|'],gen_by_name.r),

        (              gen_by_name.y0,(B6.gen_by_name['|(0->4)|'],),
            B5.gen_by_name['|(0->3)|'],gen_by_name.y3),

        (              gen_by_name.r,(B6.gen_by_name['|(2->5)|'],),
            B5.gen_by_name['|(1->4)|'],gen_by_name.y2),

        (              gen_by_name.y0,(B6.gen_by_name['|(0->3)|'],),
            B5.gen_by_name['|(0->2)|'],gen_by_name.y2),

        (              gen_by_name.y1,(B6.gen_by_name['|(1->5)|'],),
            B5.gen_by_name['|(1->4)|'],gen_by_name.y2),

        (              gen_by_name.y1,(B6.gen_by_name['|(1->4)|'],),
            B5.gen_by_name['|(1->3)|'],gen_by_name.y3),

        (              gen_by_name.y0,(B6.gen_by_name['|(0->5)|'],),
            B5.gen_by_name['|(0->4)|'],gen_by_name.y2),


        # game changer

        (              gen_by_name.y2,(B6.gen_by_name['|(5->6)|'],),
            B5.gen_by_name['|(4->6)|'],gen_by_name.y1),

        (              gen_by_name.y2,(B6.gen_by_name['|(5->7)|'],),
            B5.gen_by_name['|(4->7)|'],gen_by_name.y3),

        (              gen_by_name.y2,(B6.gen_by_name['|(3->6)|'],),
            B5.gen_by_name['|(2->6)|'],gen_by_name.y1),

        (              gen_by_name.r,(B6.gen_by_name['|(2->6)|'],),
            B5.gen_by_name['|(1->6)|'],gen_by_name.y1),

        (              gen_by_name.y3,(B6.gen_by_name['|(4->6)|'],),
            B5.gen_by_name['|(3->6)|'],gen_by_name.y1),

        (              gen_by_name.y0,(B6.gen_by_name['|(0->6)|'],),
            B5.gen_by_name['|(0->6)|'],gen_by_name.y1),

        (              gen_by_name.y1,(B6.gen_by_name['|(1->6)|'],),
            B5.gen_by_name['|(1->6)|'],gen_by_name.y1),

        (              gen_by_name.y2,(B6.gen_by_name['|(3->7)|'],),
            B5.gen_by_name['|(2->7)|'],gen_by_name.y3),

        (              gen_by_name.r,(B6.gen_by_name['|(2->7)|'],),
            B5.gen_by_name['|(1->7)|'],gen_by_name.y3),

        (              gen_by_name.y0,(B6.gen_by_name['|(0->7)|'],),
            B5.gen_by_name['|(0->7)|'],gen_by_name.y3),

        (              gen_by_name.y3,(B6.gen_by_name['|(4->7)|'],),
            B5.gen_by_name['|(3->7)|'],gen_by_name.y3),

        (              gen_by_name.y1,(B6.gen_by_name['|(1->7)|'],),
            B5.gen_by_name['|(1->7)|'],gen_by_name.y3),

        (              gen_by_name.y2,(B6.gen_by_name['|(5->6)|'],B6.gen_by_name['|(1->2)|']),
            B5.gen_by_name['|(4->5)|'],gen_by_name.y0),

        (              gen_by_name.y2,(B6.gen_by_name['|(3->6)|'],B6.gen_by_name['|(1->2)|']),
            B5.gen_by_name['|(2->5)|'],gen_by_name.y0),

        (              gen_by_name.y3,(B6.gen_by_name['|(4->6)|'],B6.gen_by_name['|(1->2)|']),
            B5.gen_by_name['|(3->5)|'],gen_by_name.y0),

        (              gen_by_name.y1,(B6.gen_by_name['|(1->6)|'],B6.gen_by_name['|(1->2)|']),
            B5.gen_by_name['|(1->5)|'],gen_by_name.y0),

        (              gen_by_name.r,(B6.gen_by_name['|(2->6)|'],B6.gen_by_name['|(1->2)|']),
            B5.gen_by_name['|(1->5)|'],gen_by_name.y0),

        (              gen_by_name.y0,(B6.gen_by_name['|(0->6)|'],B6.gen_by_name['|(1->2)|']),
            B5.gen_by_name['|(0->5)|'],gen_by_name.y0),


        ])

    return DA_bimodule(gen_by_name,arrows,B5,B6,name="N(μ_4)")
DA5=init_DA5(B5,B6)
# DA5.show_for_tex()

### DA6 arcslide, over B6=B2, B1
def init_DA6(B6,B1):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "r":Generator("r")
                })
    
    gen_by_name.x0.add_idems(B6.idem_by_name['|(0, 2)|'],B1.idem_by_name.i0)
    gen_by_name.x1.add_idems(B6.idem_by_name['|(1, 6)|'],B1.idem_by_name.i1)
    gen_by_name.x2.add_idems(B6.idem_by_name['|(3, 5)|'],B1.idem_by_name.i2)
    gen_by_name.x3.add_idems(B6.idem_by_name['|(4, 7)|'],B1.idem_by_name.i3)
    gen_by_name.r.add_idems(B6.idem_by_name['|(3, 5)|'],B1.idem_by_name.i1)

    arrows=Bunch_of_arrows([
        # short-near-chords:
        (              gen_by_name.x0,(B1.gen_by_name.r1,),
            B6.gen_by_name['|(0->1)|'],gen_by_name.x1),

        (              gen_by_name.x1,(B1.gen_by_name.r2,),
            B6.gen_by_name['|(1->2)|'],gen_by_name.x0),

        (              gen_by_name.x0,(B1.gen_by_name.r3,),
            B6.gen_by_name['|(2->3)|'],gen_by_name.r),

        (              gen_by_name.r,(B1.gen_by_name.r4,),
            1,gen_by_name.x2),

        (              gen_by_name.x2,(B1.gen_by_name.r5,),
            B6.gen_by_name['|(3->4)|'],gen_by_name.x3),

        (              gen_by_name.x3,(B1.gen_by_name.r6,),
            B6.gen_by_name['|(4->5)|'],gen_by_name.x2),

        (              gen_by_name.r,(),
            B6.gen_by_name['|(5->6)|'],gen_by_name.x1),

        # the rest, kind of follows from short ones, but not exactly
        # the game changer for DA bimodule is usually periodic domain, which
        # contains σ and the domain next to it. See below

        (              gen_by_name.x0,(B1.gen_by_name.r34,),
            B6.gen_by_name['|(2->3)|'],gen_by_name.x2),

        (              gen_by_name.x2,(B1.gen_by_name.r56,),
            B6.gen_by_name['|(3->5)|'],gen_by_name.x2),

        (              gen_by_name.x0,(B1.gen_by_name.r12,),
            B6.gen_by_name['|(0->2)|'],gen_by_name.x0),

        (              gen_by_name.r,(B1.gen_by_name.r45,),
            B6.gen_by_name['|(3->4)|'],gen_by_name.x3),

        (              gen_by_name.x1,(B1.gen_by_name.r23,),
            B6.gen_by_name['|(1->3)|'],gen_by_name.r),

        (              gen_by_name.x0,(B1.gen_by_name.r345,),
            B6.gen_by_name['|(2->4)|'],gen_by_name.x3),

        (              gen_by_name.x0,(B1.gen_by_name.r1234,),
            B6.gen_by_name['|(0->3)|'],gen_by_name.x2),

        (              gen_by_name.x0,(B1.gen_by_name.r3456,),
            B6.gen_by_name['|(2->5)|'],gen_by_name.x2),

        (              gen_by_name.r,(B1.gen_by_name.r456,),
            B6.gen_by_name['|(3->5)|'],gen_by_name.x2),

        (              gen_by_name.x1,(B1.gen_by_name.r2345,),
            B6.gen_by_name['|(1->4)|'],gen_by_name.x3),

        (              gen_by_name.x1,(B1.gen_by_name.r234,),
            B6.gen_by_name['|(1->3)|'],gen_by_name.x2),

        (              gen_by_name.x1,(B1.gen_by_name.r23456,),
            B6.gen_by_name['|(1->5)|'],gen_by_name.x2),

        (              gen_by_name.x0,(B1.gen_by_name.r12345,),
            B6.gen_by_name['|(0->4)|'],gen_by_name.x3),

        (              gen_by_name.x0,(B1.gen_by_name.r123456,),
            B6.gen_by_name['|(0->5)|'],gen_by_name.x2),

        (              gen_by_name.x0,(B1.gen_by_name.r123,),
            B6.gen_by_name['|(0->3)|'],gen_by_name.r),

        # game changer

        (              gen_by_name.x2,(B1.gen_by_name.r7,),
            B6.gen_by_name['|(5->7)|'],gen_by_name.x3),

        (              gen_by_name.x1,(B1.gen_by_name.r4,B1.gen_by_name.r7),
            B6.gen_by_name['|(6->7)|'],gen_by_name.x3),

        (              gen_by_name.x3,(B1.gen_by_name.r67,),
            B6.gen_by_name['|(4->7)|'],gen_by_name.x3),

        (              gen_by_name.x2,(B1.gen_by_name.r567,),
            B6.gen_by_name['|(3->7)|'],gen_by_name.x3),

        (              gen_by_name.x0,(B1.gen_by_name.r34567,),
            B6.gen_by_name['|(2->7)|'],gen_by_name.x3),

        (              gen_by_name.x1,(B1.gen_by_name.r234567,),
            B6.gen_by_name['|(1->7)|'],gen_by_name.x3),

        (              gen_by_name.x0,(B1.gen_by_name.r1234567,),
            B6.gen_by_name['|(0->7)|'],gen_by_name.x3),

        (              gen_by_name.r,(B1.gen_by_name.r4567,),
            B6.gen_by_name['|(3->7)|'],gen_by_name.x3),

        ])

    return DA_bimodule(gen_by_name,arrows,B6,B1,name="N(η^{-1})")
DA6=init_DA6(B6,B1)
# DA6.show_for_tex()

### FINAL COMPUTATION - WORKS!
WOW_DA=da_da_box_tensor_many_efficient_cancelations(DA1,DA2,DA3,DA4,DA5,DA6)
g2_T_LHD_red=da_randomly_cancel_until_possible(g2_T_LHD)

print("Below we compute the tensor product of 6 arc-slide bimodules, and verify that it is homotopy equivalent to N(τ_C^{-1}).\n")
print("We first describe all the 6 arc-slide bimodules:")
for X in [DA1,DA2,DA3,DA4,DA5,DA6]:
    X.show_short()
print("Next we describe the tensor product of 6 bimodules")  
WOW_DA.show_short()
print("Next we describe the bimodule N(τ_C^{-1})")  
g2_T_LHD.show_short()
print("Next we describe the bimodule N(τ_C^{-1}), but simplified with the one possible cancellation")
g2_T_LHD_red.show_short()
print("Next we check if the two bimodules in question are isomorphic:")
are_equal_smart_da(WOW_DA,g2_T_LHD_red)
print('\nType "yes" and press enter if you want more details about the bimodules and the isomorphism between them; type "no" and press enter otherwise.')
input1= input()
if input1=='yes':
    print("We first describe all the 6 arc-slide bimodules.")
    for X in [DA1,DA2,DA3,DA4,DA5,DA6]:
        X.show()
    print("Next we describe the tensor product of 6 bimodules (already maximally canceled)")  
    WOW_DA.show()
    print("Next we describe the bimodule N(τ_C^{-1}), simplified with the one possible cancellation")
    g2_T_LHD_red.show()
    print("Next we check if the two bimodules in question are isomorphic:")
    are_equal_smart_da(WOW_DA,g2_T_LHD_red,verbose=True)
    



############# genus=1 playing around with AA(id) - not important
# def init_g1_DD_id(torus_A,torus_A2):
#     gen_by_name=AttrDict({ 
#             "(i0⊗i1)": Generator("(i0⊗i1)"),
#             "(i1⊗i0)": Generator("(i1⊗i0)"),
#             })

#     gen_by_name['(i0⊗i1)'].add_idems(torus_A.idem_by_name.i0,torus_A2.idem_by_name.i1 )
#     gen_by_name['(i1⊗i0)'].add_idems(torus_A.idem_by_name.i1,torus_A2.idem_by_name.i0)

#     dd_arrows=Bunch_of_arrows([
#                     (gen_by_name['(i0⊗i1)'],
#         torus_A.gen_by_name['r1'], 
#         gen_by_name['(i1⊗i0)'], 
#         torus_A2.gen_by_name['r1']),

#                     (gen_by_name['(i0⊗i1)'],
#         torus_A.gen_by_name['r3'], 
#         gen_by_name['(i1⊗i0)'], 
#         torus_A2.gen_by_name['r3']),

#                     (gen_by_name['(i0⊗i1)'],
#         torus_A.gen_by_name['r123'], 
#         gen_by_name['(i1⊗i0)'], 
#         torus_A2.gen_by_name['r123']),

#                     (gen_by_name['(i1⊗i0)'],
#         torus_A.gen_by_name['r2'], 
#         gen_by_name['(i0⊗i1)'], 
#         torus_A2.gen_by_name['r2']),

#         ])
    
#     return DD_bimodule(gen_by_name,dd_arrows,torus_A,torus_A2,name="g1_DD_id")
# def init_g1_AA_id(torus_A,torus_A2):
#     gen_by_name=AttrDict({
#             "y": Generator("y"),
#             "x": Generator("x"),
#             "w1": Generator("w1"),
#             "w2": Generator("w2"),
#             "z1": Generator("z1"),
#             "z2": Generator("z2"),
#             })

#     gen_by_name.y.add_idems(torus_A.idem_by_name.i0,torus_A2.idem_by_name.i1 )
#     gen_by_name.x.add_idems(torus_A.idem_by_name.i1,torus_A2.idem_by_name.i0)
#     gen_by_name.z1.add_idems(torus_A.idem_by_name.i1,torus_A2.idem_by_name.i1 )
#     gen_by_name.z2.add_idems(torus_A.idem_by_name.i1,torus_A2.idem_by_name.i1)
#     gen_by_name.w1.add_idems(torus_A.idem_by_name.i0,torus_A2.idem_by_name.i0 )
#     gen_by_name.w2.add_idems(torus_A.idem_by_name.i0,torus_A2.idem_by_name.i0)

#     dd_arrows=Bunch_of_arrows([
#             ((torus_A.gen_by_name['r2'],),
#             gen_by_name.y,
#             (torus_A2.gen_by_name['r2'],),
#             gen_by_name.x),

#             ((torus_A.gen_by_name['r2'],),
#             gen_by_name.y,
#             (torus_A2.gen_by_name['r23'],),
#             gen_by_name.z2),

#             ((torus_A.gen_by_name['r12'],),
#             gen_by_name.y,
#             (torus_A2.gen_by_name['r2'],),
#             gen_by_name.w2),

#             ((),
#             gen_by_name.w1,
#             (torus_A2.gen_by_name['r1'],),
#             gen_by_name.y),

#             ((torus_A.gen_by_name['r3'],),
#             gen_by_name.z1,
#             (),
#             gen_by_name.y),

#             ((torus_A2.gen_by_name['r2'],),
#             gen_by_name.w1,
#             (torus_A2.gen_by_name['r12'],),
#             gen_by_name.x),

#             ((torus_A.gen_by_name['r23'],),
#             gen_by_name.z1,
#             (torus_A.gen_by_name['r2'],),
#             gen_by_name.x),

#             ((),
#             gen_by_name.x,
#             (torus_A2.gen_by_name['r3'],),
#             gen_by_name.z2),

#             ((torus_A.gen_by_name['r1'],),
#             gen_by_name.x,
#             (),
#             gen_by_name.w2),

#             ((),
#             gen_by_name.w1,
#             (),
#             gen_by_name.w2),

#             ((torus_A2.gen_by_name['r12'],),
#             gen_by_name.w1,
#             (torus_A2.gen_by_name['r12'],),
#             gen_by_name.w2),

#             ((),
#             gen_by_name.z1,
#             (),
#             gen_by_name.z2),

#             ((torus_A2.gen_by_name['r23'],),
#             gen_by_name.z1,
#             (torus_A2.gen_by_name['r23'],),
#             gen_by_name.z2),

#             ((torus_A2.gen_by_name['r123'],),
#             gen_by_name.z1,
#             (torus_A2.gen_by_name['r2'],),
#             gen_by_name.w2),

#             ((torus_A2.gen_by_name['r2'],),
#             gen_by_name.w1,
#             (torus_A2.gen_by_name['r123'],),
#             gen_by_name.z2),

#             ((torus_A2.gen_by_name['r2'],),
#             gen_by_name.w1,
#             (torus_A2.gen_by_name['r3'],torus_A2.gen_by_name['r2'],torus_A2.gen_by_name['r1']),
#             gen_by_name.z2),


#         ])
    
#     return AA_bimodule(gen_by_name,dd_arrows,torus_A,torus_A2,name="g1_AA_id")
# g1_AA_id=init_g1_AA_id(torus_A,torus_A)
# g1_DD_id=init_g1_DD_id(torus_A,torus_A)
# g1_DA_id=ID1
# dual_g1_DD_id=dualization_of_DD(g1_DD_id)
# a_A_a=a_A_a(torus_A)
# dual_a_A_a=dualization_of_AA(a_A_a)

# da=dd_aa_box_tensor_product(dual_g1_DD_id,a_A_a)
# aa=aa_da_box_tensor_product(dual_a_A_a,da)
# aa.show()
# aa=aa_randomly_cancel_until_possible(aa)

# mb_da_id=dd_aa_box_tensor_product(g1_DD_id,aa)
# mb_da_id=da_randomly_cancel_until_possible(mb_da_id)
# mb_da_id.show()
##### it is DA(id), but it has too many generators


