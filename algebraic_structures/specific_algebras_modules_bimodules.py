# -*- coding: utf-8 -*- 
from .basics import AttrDict, Bunch_of_arrows, Generator, debug
# from sets import Set
from .algebra import *
from .aa_bimodule import *
from .dd_bimodule import *

########## Algebras

def init_torus_algebra():
    gen_by_name=AttrDict({
                "r1": Generator("r1"),
                "r2": Generator("r2"),
                "r3": Generator("r3"),
                "r12": Generator("r12"),
                "r23": Generator("r23"),
                "r123": Generator("r123"),
                "i0": Generator("i0"),
                "i1": Generator("i1"),
                })

    # now this stuff is done when algebra is initiated
    # gen_by_name.r12.add_factorizations((gen_by_name.r1,gen_by_name.r2))
    # gen_by_name.r23.add_factorizations((gen_by_name.r2,gen_by_name.r3))
    # gen_by_name.r123.add_factorizations((gen_by_name.r1,gen_by_name.r23),(gen_by_name.r12,gen_by_name.r3))

    idem_by_name=AttrDict({
                    "i0": gen_by_name.i0,
                    "i1": gen_by_name.i1,
                    })
    multiplication_table={(gen_by_name.r1,gen_by_name.r2):gen_by_name.r12,
                        (gen_by_name.r2,gen_by_name.r3):gen_by_name.r23,
                        (gen_by_name.r1,gen_by_name.r23):gen_by_name.r123,
                        (gen_by_name.r12,gen_by_name.r3):gen_by_name.r123,
                        (gen_by_name.i0,gen_by_name.r1):gen_by_name.r1,
                        (gen_by_name.i0,gen_by_name.r3):gen_by_name.r3,
                        (gen_by_name.i0,gen_by_name.r12):gen_by_name.r12,
                        (gen_by_name.i0,gen_by_name.r123):gen_by_name.r123,
                        (gen_by_name.i0,gen_by_name.i0):gen_by_name.i0,
                        (gen_by_name.i1,gen_by_name.i1):gen_by_name.i1,
                        (gen_by_name.i1,gen_by_name.r23):gen_by_name.r23,
                        (gen_by_name.i1,gen_by_name.r2):gen_by_name.r2,
                        (gen_by_name.r1,gen_by_name.i1):gen_by_name.r1,
                        (gen_by_name.r2,gen_by_name.i0):gen_by_name.r2,
                        (gen_by_name.r3,gen_by_name.i1):gen_by_name.r3,
                        (gen_by_name.r12,gen_by_name.i0):gen_by_name.r12,
                        (gen_by_name.r23,gen_by_name.i1):gen_by_name.r23,
                        (gen_by_name.r123,gen_by_name.i1):gen_by_name.r123,
                                    }

    return dg_algebra(gen_by_name=gen_by_name,idem_by_name=idem_by_name,multiplication_table=multiplication_table,name='torus_A')
torus_A=init_torus_algebra()
def init_torus_algebra_2bp():
    gen_by_name=AttrDict({
                "r1": Generator("r1"),
                "r2": Generator("r2"),
                "r12": Generator("r12"),
                "ks1": Generator("ks1"),
                "ks2": Generator("ks2"),
                "ks12": Generator("ks12"),
                "i0": Generator("i0"),
                "i1": Generator("i1"),
                "i2": Generator("i2"),
                })

    # now this stuff is done when algebra is initiated
    # gen_by_name.r12.add_factorizations((gen_by_name.r1,gen_by_name.r2))
    # gen_by_name.r23.add_factorizations((gen_by_name.r2,gen_by_name.r3))
    # gen_by_name.r123.add_factorizations((gen_by_name.r1,gen_by_name.r23),(gen_by_name.r12,gen_by_name.r3))

    idem_by_name=AttrDict({
                    "i0": gen_by_name.i0,
                    "i1": gen_by_name.i1,
                    "i2": gen_by_name.i2,
                    })
    multiplication_table={
                        (gen_by_name.i0,gen_by_name.r1):gen_by_name.r1,
                        (gen_by_name.i0,gen_by_name.r12):gen_by_name.r12,
                        (gen_by_name.i0,gen_by_name.i0):gen_by_name.i0,
                        (gen_by_name.i0,gen_by_name.ks1):gen_by_name.ks1,
                        (gen_by_name.i0,gen_by_name.ks12):gen_by_name.ks12,
                        (gen_by_name.i1,gen_by_name.i1):gen_by_name.i1,
                        (gen_by_name.i1,gen_by_name.r2):gen_by_name.r2,
                        (gen_by_name.i1,gen_by_name.ks2):gen_by_name.ks2,
                        (gen_by_name.i2,gen_by_name.i2):gen_by_name.i2,
                        (gen_by_name.r2,gen_by_name.i2):gen_by_name.r2,
                        (gen_by_name.r12,gen_by_name.i2):gen_by_name.r12,
                        (gen_by_name.ks12,gen_by_name.i2):gen_by_name.ks12,
                        (gen_by_name.ks2,gen_by_name.i2):gen_by_name.ks2,
                        (gen_by_name.r1,gen_by_name.i1):gen_by_name.r1,
                        (gen_by_name.ks1,gen_by_name.i1):gen_by_name.ks1,
                        (gen_by_name.r1,gen_by_name.r2):gen_by_name.r12,
                        (gen_by_name.ks1,gen_by_name.ks2):gen_by_name.ks12,
                                    }
    return dg_algebra(gen_by_name=gen_by_name,idem_by_name=idem_by_name,multiplication_table=multiplication_table,name='torus_A')
torus_2bp_A=init_torus_algebra_2bp()
def init_genus2_algebra():
    gen_by_name=AttrDict({
                "r1": Generator("r1"),
                "r2": Generator("r2"),
                "r3": Generator("r3"),
                "r12": Generator("r12"),
                "r23": Generator("r23"),
                "r123": Generator("r123"),

                "i0": Generator("i0"),
                "i1": Generator("i1"),

                "r5": Generator("r5"),
                "r6": Generator("r6"),
                "r7": Generator("r7"),
                "r56": Generator("r56"),
                "r67": Generator("r67"),
                "r567": Generator("r567"),

                "i2": Generator("i2"),
                "i3": Generator("i3"),

                "r1234": Generator("r1234"),
                "r12345": Generator("r12345"),
                "r123456": Generator("r123456"),
                "r1234567": Generator("r1234567"),
                "r234": Generator("r234"),
                "r2345": Generator("r2345"),
                "r23456": Generator("r23456"),
                "r234567": Generator("r234567"),
                "r34": Generator("r34"),
                "r345": Generator("r345"),
                "r3456": Generator("r3456"),
                "r34567": Generator("r34567"),
                "r4": Generator("r4"),
                "r45": Generator("r45"),
                "r456": Generator("r456"),
                "r4567": Generator("r4567")
                })

    # now this stuff is done when algebra is initiated
    # gen_by_name.r12.add_factorizations((gen_by_name.r1,gen_by_name.r2))
    # gen_by_name.r23.add_factorizations((gen_by_name.r2,gen_by_name.r3))
    # gen_by_name.r34.add_factorizations((gen_by_name.r3,gen_by_name.r4))

    # gen_by_name.r45.add_factorizations((gen_by_name.r4,gen_by_name.r5))
    # gen_by_name.r56.add_factorizations((gen_by_name.r5,gen_by_name.r6))
    # gen_by_name.r67.add_factorizations((gen_by_name.r6,gen_by_name.r7))

    # gen_by_name.r123.add_factorizations((gen_by_name.r1,gen_by_name.r23),(gen_by_name.r12,gen_by_name.r3))
    # gen_by_name.r234.add_factorizations((gen_by_name.r2,gen_by_name.r34),(gen_by_name.r23,gen_by_name.r4))
    # gen_by_name.r345.add_factorizations((gen_by_name.r3,gen_by_name.r45),(gen_by_name.r34,gen_by_name.r5))
    # gen_by_name.r456.add_factorizations((gen_by_name.r4,gen_by_name.r56),(gen_by_name.r45,gen_by_name.r6))
    # gen_by_name.r567.add_factorizations((gen_by_name.r5,gen_by_name.r67),(gen_by_name.r56,gen_by_name.r7))

    # gen_by_name.r1234.add_factorizations((gen_by_name.r1,gen_by_name.r234),(gen_by_name.r12,gen_by_name.r34),
    #                                     (gen_by_name.r123,gen_by_name.r4))
    # gen_by_name.r2345.add_factorizations((gen_by_name.r2,gen_by_name.r345),(gen_by_name.r23,gen_by_name.r45),
    #                                     (gen_by_name.r234,gen_by_name.r5))
    # gen_by_name.r3456.add_factorizations((gen_by_name.r3,gen_by_name.r456),(gen_by_name.r34,gen_by_name.r56),
    #                                     (gen_by_name.r345,gen_by_name.r6))
    # gen_by_name.r4567.add_factorizations((gen_by_name.r4,gen_by_name.r567),(gen_by_name.r45,gen_by_name.r67),
    #                                     (gen_by_name.r456,gen_by_name.r7))

    # gen_by_name.r12345.add_factorizations((gen_by_name.r1,gen_by_name.r2345),(gen_by_name.r12,gen_by_name.r345),
    #                                     (gen_by_name.r123,gen_by_name.r45),(gen_by_name.r1234,gen_by_name.r5))
    # gen_by_name.r23456.add_factorizations((gen_by_name.r2,gen_by_name.r3456),(gen_by_name.r23,gen_by_name.r456),
    #                                     (gen_by_name.r234,gen_by_name.r56),(gen_by_name.r2345,gen_by_name.r6))
    # gen_by_name.r34567.add_factorizations((gen_by_name.r3,gen_by_name.r4567),(gen_by_name.r34,gen_by_name.r567),
    #                                     (gen_by_name.r345,gen_by_name.r67),(gen_by_name.r3456,gen_by_name.r7))   


    # gen_by_name.r123456.add_factorizations((gen_by_name.r1,gen_by_name.r23456),(gen_by_name.r12,gen_by_name.r3456),
    #                                     (gen_by_name.r123,gen_by_name.r456),(gen_by_name.r1234,gen_by_name.r56),
    #                                     (gen_by_name.r12345,gen_by_name.r6))
    # gen_by_name.r234567.add_factorizations((gen_by_name.r2,gen_by_name.r34567),(gen_by_name.r23,gen_by_name.r4567),
    #                                     (gen_by_name.r234,gen_by_name.r567),(gen_by_name.r2345,gen_by_name.r67),
    #                                     (gen_by_name.r23456,gen_by_name.r7)) 

    # gen_by_name.r1234567.add_factorizations((gen_by_name.r1,gen_by_name.r234567),(gen_by_name.r12,gen_by_name.r34567),
    #                                     (gen_by_name.r123,gen_by_name.r4567),(gen_by_name.r1234,gen_by_name.r567),
    #                                     (gen_by_name.r12345,gen_by_name.r67),(gen_by_name.r123456,gen_by_name.r7)) 


    idem_by_name=AttrDict({
                    "i0": gen_by_name.i0,
                    "i1": gen_by_name.i1,
                    "i2": gen_by_name.i2,
                    "i3": gen_by_name.i3,
                    })
    multiplication_table={(gen_by_name.r1,gen_by_name.r2):gen_by_name.r12, #first torus algebra part
                        (gen_by_name.r2,gen_by_name.r3):gen_by_name.r23,
                        (gen_by_name.r1,gen_by_name.r23):gen_by_name.r123,
                        (gen_by_name.r12,gen_by_name.r3):gen_by_name.r123,
                        (gen_by_name.i0,gen_by_name.r1):gen_by_name.r1,
                        (gen_by_name.i0,gen_by_name.r3):gen_by_name.r3,
                        (gen_by_name.i0,gen_by_name.r12):gen_by_name.r12,
                        (gen_by_name.i0,gen_by_name.r123):gen_by_name.r123,
                        (gen_by_name.i0,gen_by_name.i0):gen_by_name.i0,
                        (gen_by_name.i1,gen_by_name.i1):gen_by_name.i1,
                        (gen_by_name.i1,gen_by_name.r23):gen_by_name.r23,
                        (gen_by_name.i1,gen_by_name.r2):gen_by_name.r2,
                        (gen_by_name.r1,gen_by_name.i1):gen_by_name.r1,
                        (gen_by_name.r2,gen_by_name.i0):gen_by_name.r2,
                        (gen_by_name.r3,gen_by_name.i1):gen_by_name.r3,
                        (gen_by_name.r12,gen_by_name.i0):gen_by_name.r12,
                        (gen_by_name.r23,gen_by_name.i1):gen_by_name.r23,
                        (gen_by_name.r123,gen_by_name.i1):gen_by_name.r123,

                        (gen_by_name.r5,gen_by_name.r6):gen_by_name.r56, #second torus algebra part
                        (gen_by_name.r6,gen_by_name.r7):gen_by_name.r67,
                        (gen_by_name.r5,gen_by_name.r67):gen_by_name.r567,
                        (gen_by_name.r56,gen_by_name.r7):gen_by_name.r567,
                        (gen_by_name.i2,gen_by_name.r5):gen_by_name.r5,
                        (gen_by_name.i2,gen_by_name.r7):gen_by_name.r7,
                        (gen_by_name.i2,gen_by_name.r56):gen_by_name.r56,
                        (gen_by_name.i2,gen_by_name.r567):gen_by_name.r567,
                        (gen_by_name.i2,gen_by_name.i2):gen_by_name.i2,
                        (gen_by_name.i3,gen_by_name.i3):gen_by_name.i3,
                        (gen_by_name.i3,gen_by_name.r67):gen_by_name.r67,
                        (gen_by_name.i3,gen_by_name.r6):gen_by_name.r6,
                        (gen_by_name.r5,gen_by_name.i3):gen_by_name.r5,
                        (gen_by_name.r6,gen_by_name.i2):gen_by_name.r6,
                        (gen_by_name.r7,gen_by_name.i3):gen_by_name.r7,
                        (gen_by_name.r56,gen_by_name.i2):gen_by_name.r56,
                        (gen_by_name.r67,gen_by_name.i3):gen_by_name.r67,
                        (gen_by_name.r567,gen_by_name.i3):gen_by_name.r567,

                        (gen_by_name.i0,gen_by_name.r1234):gen_by_name.r1234, #rest of idempotent multiplications
                        (gen_by_name.i0,gen_by_name.r12345):gen_by_name.r12345,
                        (gen_by_name.i0,gen_by_name.r123456):gen_by_name.r123456,
                        (gen_by_name.i0,gen_by_name.r1234567):gen_by_name.r1234567,
                        (gen_by_name.i0,gen_by_name.r34):gen_by_name.r34,
                        (gen_by_name.i0,gen_by_name.r345):gen_by_name.r345,
                        (gen_by_name.i0,gen_by_name.r3456):gen_by_name.r3456,
                        (gen_by_name.i0,gen_by_name.r34567):gen_by_name.r34567,

                        (gen_by_name.i1,gen_by_name.r234):gen_by_name.r234,
                        (gen_by_name.i1,gen_by_name.r2345):gen_by_name.r2345,
                        (gen_by_name.i1,gen_by_name.r23456):gen_by_name.r23456,
                        (gen_by_name.i1,gen_by_name.r234567):gen_by_name.r234567,
                        (gen_by_name.i1,gen_by_name.r4):gen_by_name.r4,
                        (gen_by_name.i1,gen_by_name.r45):gen_by_name.r45,
                        (gen_by_name.i1,gen_by_name.r456):gen_by_name.r456,
                        (gen_by_name.i1,gen_by_name.r4567):gen_by_name.r4567,    

                        (gen_by_name.r1234,gen_by_name.i2):gen_by_name.r1234,
                        (gen_by_name.r12345,gen_by_name.i3):gen_by_name.r12345,
                        (gen_by_name.r123456,gen_by_name.i2):gen_by_name.r123456,
                        (gen_by_name.r1234567,gen_by_name.i3):gen_by_name.r1234567,
                        (gen_by_name.r34,gen_by_name.i2):gen_by_name.r34,
                        (gen_by_name.r345,gen_by_name.i3):gen_by_name.r345,
                        (gen_by_name.r3456,gen_by_name.i2):gen_by_name.r3456,
                        (gen_by_name.r34567,gen_by_name.i3):gen_by_name.r34567,

                        (gen_by_name.r234,gen_by_name.i2):gen_by_name.r234,
                        (gen_by_name.r2345,gen_by_name.i3):gen_by_name.r2345,
                        (gen_by_name.r23456,gen_by_name.i2):gen_by_name.r23456,
                        (gen_by_name.r234567,gen_by_name.i3):gen_by_name.r234567,
                        (gen_by_name.r4,gen_by_name.i2):gen_by_name.r4,
                        (gen_by_name.r45,gen_by_name.i3):gen_by_name.r45,
                        (gen_by_name.r456,gen_by_name.i2):gen_by_name.r456,
                        (gen_by_name.r4567,gen_by_name.i3):gen_by_name.r4567,  

                        (gen_by_name.r4,gen_by_name.r5):gen_by_name.r45, #actual elements (havinig r4 inside) factorizations
                        (gen_by_name.r3,gen_by_name.r4):gen_by_name.r34, 
                        (gen_by_name.r1,gen_by_name.r23):gen_by_name.r123,
                        (gen_by_name.r12,gen_by_name.r3):gen_by_name.r123,
                        (gen_by_name.r2,gen_by_name.r34):gen_by_name.r234,
                        (gen_by_name.r23,gen_by_name.r4):gen_by_name.r234,
                        (gen_by_name.r3,gen_by_name.r45):gen_by_name.r345,
                        (gen_by_name.r34,gen_by_name.r5):gen_by_name.r345,
                        (gen_by_name.r4,gen_by_name.r56):gen_by_name.r456,
                        (gen_by_name.r45,gen_by_name.r6):gen_by_name.r456,
                        (gen_by_name.r5,gen_by_name.r67):gen_by_name.r567,
                        (gen_by_name.r56,gen_by_name.r7):gen_by_name.r567,

  
                        (gen_by_name.r1,gen_by_name.r234):gen_by_name.r1234,
                        (gen_by_name.r12,gen_by_name.r34):gen_by_name.r1234,
                        (gen_by_name.r123,gen_by_name.r4):gen_by_name.r1234,

                        (gen_by_name.r2,gen_by_name.r345):gen_by_name.r2345,
                        (gen_by_name.r23,gen_by_name.r45):gen_by_name.r2345,
                        (gen_by_name.r234,gen_by_name.r5):gen_by_name.r2345,

                        (gen_by_name.r3,gen_by_name.r456):gen_by_name.r3456,
                        (gen_by_name.r34,gen_by_name.r56):gen_by_name.r3456,
                        (gen_by_name.r345,gen_by_name.r6):gen_by_name.r3456,

                        (gen_by_name.r4,gen_by_name.r567):gen_by_name.r4567,
                        (gen_by_name.r45,gen_by_name.r67):gen_by_name.r4567,
                        (gen_by_name.r456,gen_by_name.r7):gen_by_name.r4567,

                        (gen_by_name.r1,gen_by_name.r2345):gen_by_name.r12345,
                        (gen_by_name.r12,gen_by_name.r345):gen_by_name.r12345,
                        (gen_by_name.r123,gen_by_name.r45):gen_by_name.r12345,
                        (gen_by_name.r1234,gen_by_name.r5):gen_by_name.r12345,

                        (gen_by_name.r2,gen_by_name.r3456):gen_by_name.r23456,
                        (gen_by_name.r23,gen_by_name.r456):gen_by_name.r23456,
                        (gen_by_name.r234,gen_by_name.r56):gen_by_name.r23456,
                        (gen_by_name.r2345,gen_by_name.r6):gen_by_name.r23456,

                        (gen_by_name.r3,gen_by_name.r4567):gen_by_name.r34567,
                        (gen_by_name.r34,gen_by_name.r567):gen_by_name.r34567,
                        (gen_by_name.r345,gen_by_name.r67):gen_by_name.r34567,
                        (gen_by_name.r3456,gen_by_name.r7):gen_by_name.r34567,   

                        (gen_by_name.r1,gen_by_name.r23456):gen_by_name.r123456,
                        (gen_by_name.r12,gen_by_name.r3456):gen_by_name.r123456,
                        (gen_by_name.r123,gen_by_name.r456):gen_by_name.r123456,
                        (gen_by_name.r1234,gen_by_name.r56):gen_by_name.r123456,
                        (gen_by_name.r12345,gen_by_name.r6):gen_by_name.r123456,

                        (gen_by_name.r2,gen_by_name.r34567):gen_by_name.r234567,
                        (gen_by_name.r23,gen_by_name.r4567):gen_by_name.r234567,
                        (gen_by_name.r234,gen_by_name.r567):gen_by_name.r234567,
                        (gen_by_name.r2345,gen_by_name.r67):gen_by_name.r234567,
                        (gen_by_name.r23456,gen_by_name.r7):gen_by_name.r234567, 

                        (gen_by_name.r1,gen_by_name.r234567):gen_by_name.r1234567,
                        (gen_by_name.r12,gen_by_name.r34567):gen_by_name.r1234567,
                        (gen_by_name.r123,gen_by_name.r4567):gen_by_name.r1234567,
                        (gen_by_name.r1234,gen_by_name.r567):gen_by_name.r1234567,
                        (gen_by_name.r12345,gen_by_name.r67):gen_by_name.r1234567,
                        (gen_by_name.r123456,gen_by_name.r7):gen_by_name.r1234567
                
                                    }

    return dg_algebra(gen_by_name=gen_by_name,idem_by_name=idem_by_name,multiplication_table=multiplication_table,name='g2_A')
g2_A=init_genus2_algebra()
def init_pillowcase_algebra():
    gen_by_name=AttrDict({
                "i0": Generator("i0","i_0"),
                "i1": Generator("i1","i_1"),
                "i2": Generator("i2","i_2"),
                "j0": Generator("j0","j_0"),
                "j1": Generator("j1","j_1"),
                "j2": Generator("j2","j_2"),

                "r0": Generator("r0","r_0","r'_0"),
                "r2": Generator("r2","r_2","r'_2"),

                "et1": Generator("et1","\\eta_1","\\eta'_1"),
                "et2": Generator("et2","\\eta_2","\\eta'_2"),
                "et3": Generator("et3","\\eta_3","\\eta'_3"),
                "et12": Generator("et12","\\eta_{12}","\\eta'_{21}"),
                "et23": Generator("et23","\\eta_{23}","\\eta'_{32}"),
                "et123": Generator("et123","\\eta_{123}","\\eta'_{321}"),

                "ks1": Generator("ks1","\\xi_1","\\xi'_1"),
                "ks2": Generator("ks2","\\xi_2","\\xi'_2"),
                "ks3": Generator("ks3","\\xi_3","\\xi'_3"),
                "ks12": Generator("ks12","\\xi_{12}","\\xi'_{21}"),
                "ks23": Generator("ks23","\\xi_{23}","\\xi'_{32}"),
                "ks123": Generator("ks123","\\xi_{123}","\\xi'_{321}"),
                })


    # now this stuff is done when algebra is initiated
    # gen_by_name.et12.add_factorizations((gen_by_name.et1,gen_by_name.et2))
    # gen_by_name.et23.add_factorizations((gen_by_name.et2,gen_by_name.et3))
    # gen_by_name.et123.add_factorizations((gen_by_name.et1,gen_by_name.et23),(gen_by_name.et12,gen_by_name.et3))
    # gen_by_name.ks12.add_factorizations((gen_by_name.ks1,gen_by_name.ks2))
    # gen_by_name.ks23.add_factorizations((gen_by_name.ks2,gen_by_name.ks3))
    # gen_by_name.ks123.add_factorizations((gen_by_name.ks1,gen_by_name.ks23),(gen_by_name.ks12,gen_by_name.ks3))

    idem_by_name=AttrDict({
                    "i0": gen_by_name.i0,
                    "i1": gen_by_name.i1,
                    "i2": gen_by_name.i2,
                    "j0": gen_by_name.j0,
                    "j1": gen_by_name.j1,
                    "j2": gen_by_name.j2
                    })

    multiplication_table={(gen_by_name.et1,gen_by_name.et2):gen_by_name.et12,
                        (gen_by_name.et2,gen_by_name.et3):gen_by_name.et23,
                        (gen_by_name.et1,gen_by_name.et23):gen_by_name.et123,
                        (gen_by_name.et12,gen_by_name.et3):gen_by_name.et123,

                        (gen_by_name.ks1,gen_by_name.ks2):gen_by_name.ks12,
                        (gen_by_name.ks2,gen_by_name.ks3):gen_by_name.ks23,
                        (gen_by_name.ks1,gen_by_name.ks23):gen_by_name.ks123,
                        (gen_by_name.ks12,gen_by_name.ks3):gen_by_name.ks123,

                        (gen_by_name.i0,gen_by_name.r0):gen_by_name.r0,
                        (gen_by_name.i0,gen_by_name.i0):gen_by_name.i0,
                        (gen_by_name.i0,gen_by_name.et1):gen_by_name.et1,
                        (gen_by_name.i0,gen_by_name.et12):gen_by_name.et12,
                        (gen_by_name.i0,gen_by_name.et123):gen_by_name.et123,

                        (gen_by_name.r0,gen_by_name.j0):gen_by_name.r0,
                        (gen_by_name.j0,gen_by_name.j0):gen_by_name.j0,
                        (gen_by_name.et3,gen_by_name.j0):gen_by_name.et3,
                        (gen_by_name.et23,gen_by_name.j0):gen_by_name.et23,
                        (gen_by_name.et123,gen_by_name.j0):gen_by_name.et123,

                        (gen_by_name.r2,gen_by_name.j2):gen_by_name.r2,
                        (gen_by_name.j2,gen_by_name.j2):gen_by_name.j2,
                        (gen_by_name.j2,gen_by_name.ks3):gen_by_name.ks3,
                        (gen_by_name.ks2,gen_by_name.j2):gen_by_name.ks2,
                        (gen_by_name.ks12,gen_by_name.j2):gen_by_name.ks12,

                        (gen_by_name.i2,gen_by_name.r2):gen_by_name.r2,
                        (gen_by_name.i2,gen_by_name.i2):gen_by_name.i2,
                        (gen_by_name.i2,gen_by_name.ks2):gen_by_name.ks2,
                        (gen_by_name.i2,gen_by_name.ks23):gen_by_name.ks23,
                        (gen_by_name.ks1,gen_by_name.i2):gen_by_name.ks1,

                        (gen_by_name.i1,gen_by_name.i1):gen_by_name.i1,
                        (gen_by_name.i1,gen_by_name.ks1):gen_by_name.ks1,
                        (gen_by_name.i1,gen_by_name.ks12):gen_by_name.ks12,
                        (gen_by_name.i1,gen_by_name.ks123):gen_by_name.ks123,
                        (gen_by_name.i1,gen_by_name.et2):gen_by_name.et2,
                        (gen_by_name.i1,gen_by_name.et23):gen_by_name.et23,
                        (gen_by_name.et1,gen_by_name.i1):gen_by_name.et1,

                        (gen_by_name.j1,gen_by_name.j1):gen_by_name.j1,
                        (gen_by_name.j1,gen_by_name.et3):gen_by_name.et3,
                        (gen_by_name.et12,gen_by_name.j1):gen_by_name.et12,
                        (gen_by_name.et2,gen_by_name.j1):gen_by_name.et2,
                        (gen_by_name.ks3,gen_by_name.j1):gen_by_name.ks3,
                        (gen_by_name.ks23,gen_by_name.j1):gen_by_name.ks23,
                        (gen_by_name.ks123,gen_by_name.j1):gen_by_name.ks123,
                                    }

    return dg_algebra(gen_by_name=gen_by_name,idem_by_name=idem_by_name,multiplication_table=multiplication_table,name='pil_A')
pil_A=init_pillowcase_algebra()
def init_B_r():
    gen_by_name=AttrDict({
                "l0": Generator("l0","l_0"),
                "l1": Generator("l1","l_1"),
                "p01": Generator("p01","p_\{01\}"),
                "q10": Generator("q10","p_\{10\}"),
                "c0": Generator("c0","c_0"),
                "c1": Generator("c1","c_1"),
                })


    idem_by_name=AttrDict({
                    "l0": gen_by_name.l0,
                    "l1": gen_by_name.l1,
                    })

    multiplication_table={(gen_by_name.l0,gen_by_name.l0):gen_by_name.l0,
                        (gen_by_name.l0,gen_by_name.p01):gen_by_name.p01,
                        (gen_by_name.l0,gen_by_name.c0):gen_by_name.c0,
                        (gen_by_name.c0,gen_by_name.l0):gen_by_name.c0,
                        (gen_by_name.l1,gen_by_name.l1):gen_by_name.l1,
                        (gen_by_name.l1,gen_by_name.q10):gen_by_name.q10,
                        (gen_by_name.l1,gen_by_name.c1):gen_by_name.c1,
                        (gen_by_name.c1,gen_by_name.l1):gen_by_name.c1,
                        (gen_by_name.p01,gen_by_name.q10):gen_by_name.c0,
                        (gen_by_name.p01,gen_by_name.l1):gen_by_name.p01,
                        (gen_by_name.q10,gen_by_name.p01):gen_by_name.c1,
                        (gen_by_name.q10,gen_by_name.l0):gen_by_name.q10,
                                    }

    return dg_algebra(gen_by_name=gen_by_name,idem_by_name=idem_by_name,multiplication_table=multiplication_table,name='B_r')
B_r=init_B_r()

def init_Apq():
    gen_by_name=AttrDict({
                "a": Generator("a","a"),
                "b": Generator("b","b"),
                "c": Generator("c","c"),
                "d": Generator("d","d"),
                "p2": Generator("p2","p2"),
                "p3": Generator("p3","p3"),
                "p4": Generator("p4","p4"),
                "p32": Generator("p32","p32"),
                "p43": Generator("p43","p43"),
                "p432": Generator("p432","p432"),
                "q2": Generator("q2","q2"),
                "q3": Generator("q3","q3"),
                "q4": Generator("q4","q4"),
                "q23": Generator("q23","q23"),
                "q34": Generator("q34","q34"),
                "q234": Generator("q234","q234"),
                })


    idem_by_name=AttrDict({
                    "a": gen_by_name.a,
                    "b": gen_by_name.b,
                    "c": gen_by_name.c,
                    "d": gen_by_name.d,
                    })

    multiplication_table={(gen_by_name.a,gen_by_name.a):gen_by_name.a,
                        (gen_by_name.b,gen_by_name.b):gen_by_name.b,
                        (gen_by_name.c,gen_by_name.c):gen_by_name.c,
                        (gen_by_name.d,gen_by_name.d):gen_by_name.d,
                        (gen_by_name.a,gen_by_name.q2):gen_by_name.q2,
                        (gen_by_name.q2,gen_by_name.b):gen_by_name.q2,
                        (gen_by_name.b,gen_by_name.q3):gen_by_name.q3,
                        (gen_by_name.q3,gen_by_name.c):gen_by_name.q3,
                        (gen_by_name.c,gen_by_name.q4):gen_by_name.q4,
                        (gen_by_name.q4,gen_by_name.d):gen_by_name.q4,
                        (gen_by_name.a,gen_by_name.q23):gen_by_name.q23,
                        (gen_by_name.q23,gen_by_name.c):gen_by_name.q23,
                        (gen_by_name.b,gen_by_name.q34):gen_by_name.q34,
                        (gen_by_name.q34,gen_by_name.d):gen_by_name.q34,
                        (gen_by_name.a,gen_by_name.q234):gen_by_name.q234,
                        (gen_by_name.q234,gen_by_name.d):gen_by_name.q234,
                        (gen_by_name.b,gen_by_name.p2):gen_by_name.p2,
                        (gen_by_name.p2,gen_by_name.a):gen_by_name.p2,
                        (gen_by_name.c,gen_by_name.p3):gen_by_name.p3,
                        (gen_by_name.p3,gen_by_name.b):gen_by_name.p3,
                        (gen_by_name.d,gen_by_name.p4):gen_by_name.p4,
                        (gen_by_name.p4,gen_by_name.c):gen_by_name.p4,
                        (gen_by_name.c,gen_by_name.p32):gen_by_name.p32,
                        (gen_by_name.p32,gen_by_name.a):gen_by_name.p32,
                        (gen_by_name.d,gen_by_name.p43):gen_by_name.p43,
                        (gen_by_name.p43,gen_by_name.b):gen_by_name.p43,
                        (gen_by_name.d,gen_by_name.p432):gen_by_name.p432,
                        (gen_by_name.p432,gen_by_name.a):gen_by_name.p432,
                        (gen_by_name.q2,gen_by_name.q3):gen_by_name.q23,
                        (gen_by_name.q3,gen_by_name.q4):gen_by_name.q34,
                        (gen_by_name.q23,gen_by_name.q4):gen_by_name.q234,
                        (gen_by_name.q2,gen_by_name.q34):gen_by_name.q234,
                        (gen_by_name.p4,gen_by_name.p3):gen_by_name.p43,
                        (gen_by_name.p3,gen_by_name.p2):gen_by_name.p32,
                        (gen_by_name.p4,gen_by_name.p32):gen_by_name.p432,
                        (gen_by_name.p43,gen_by_name.p2):gen_by_name.p432,
                                    }

    return dg_algebra(gen_by_name=gen_by_name,idem_by_name=idem_by_name,multiplication_table=multiplication_table,name='B_r')
Apq=init_Apq()


########### bimodules 

# dual(A^bar^A) DD-bimodule from Mor paper, for algebra with no differential
def init_DD_bar_dual(directed_A):
    vertices=directed_A.idemset
    edges=directed_A.genset
    #the following works only for graphs without cycles
    def paths(idem1,idem2):
        paths_to_add=set()
        if idem1==idem2: 
            paths_to_add.add((idem1,))
            return paths_to_add
        
        outgoing_edges=[edge for edge in (set(edges)-set(vertices)) if edge.idem.left==idem1]
        for edge in outgoing_edges:
            if edge.idem.right==idem2:
                paths_to_add.add((edge,))
            else:
                auxilary_paths=paths(edge.idem.right,idem2)
                for a_path in auxilary_paths:  
                    paths_to_add.add((edge,)+a_path)
        return paths_to_add
    def reverse_str(path):
        if path[0].idem.left==path[0].idem.right: out = '('+str(path[0])+"'" +')'
        else:
            out='('
            for edge in reversed(list(path)):
                out=out+'-'+str(edge)+','
            out=out[:-1]
            out=out+')'
        return out
    def reverse_tex_str(path):
        if path[0].idem.left==path[0].idem.right: out = '('+str(path[0].tex_name)+"'" +')'
        else:
            out='('
            for edge in reversed(list(path)):
                out=out+str(edge.reverse_tex_name)+','
            out=out[:-1]
            out=out+')'
        return out

    all_paths_in_pillow_graph=set()
    for idem1 in directed_A.idemset:
        for idem2 in directed_A.idemset:
            all_paths_in_pillow_graph=all_paths_in_pillow_graph.union(paths(idem1,idem2))

    gen_by_name=AttrDict({})
    for path in all_paths_in_pillow_graph:
        gen_by_name['b' + reverse_str(path)]=Generator('b' + reverse_str(path),'b' + reverse_tex_str(path))

        idemleft=path[-1].idem.right
        idemright=path[0].idem.left
        gen_by_name['b' + reverse_str(path)].add_idems(idemleft,idemright) #because A^bar^A reverses idempotents, and dual reverses one more time

    dd_arrows=Bunch_of_arrows([])
    for path in all_paths_in_pillow_graph:
        in_path_gen=gen_by_name['b'+reverse_str(path)]
        for chord in (set(directed_A.genset)-set(directed_A.idemset)):
            rightaction_out=(chord.idem.right==path[0].idem.left)
            if rightaction_out:
                if path[0] in (directed_A.idemset):
                    out_path=(chord,)
                else:
                    out_path=(chord,)+path
                out_path_gen=gen_by_name['b'+reverse_str(out_path)]
                dd_arrows[(in_path_gen,1,out_path_gen,chord)]+=1

            leftaction_out=(chord.idem.left==path[-1].idem.right)
            if leftaction_out:
                if path[0] in (directed_A.idemset):
                    out_path=(chord,)
                else:
                    out_path=path+(chord,)
                out_path_gen=gen_by_name['b'+reverse_str(out_path)]
                dd_arrows[(in_path_gen,chord,out_path_gen,1)]+=1

        for ind, chord_in_path in enumerate(path):
            for factorization in getattr(chord_in_path,'factorizations', []):
                out_path=path[:ind]+factorization+path[ind+1:]
                out_path_gen=gen_by_name['b'+reverse_str(out_path)]
                dd_arrows[(in_path_gen,1,out_path_gen,1)]+=1
    
    return DD_bimodule(gen_by_name,dd_arrows,directed_A,directed_A,name="DD_bar_dual")
DD_bar_dual=init_DD_bar_dual(pil_A)
DD_bar_r_dual=dd_randomly_cancel_until_possible(DD_bar_dual)

