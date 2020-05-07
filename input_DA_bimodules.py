# -*- coding: utf-8 -*- 
import sys
# sys.path.append('../')

from algebraic_structures.algebra import (
    AttrDict, Generator, )
from algebraic_structures.specific_algebras_modules_bimodules import torus_A, g2_A
from algebraic_structures.da_bimodule import (
    Bunch_of_arrows, DA_bimodule)

from algebraic_structures.aa_bimodule import (
    AA_bimodule)

from algebraic_structures.dd_bimodule import (
    DD_bimodule)

# from evarist.models import (solution_filters, events,
#                             parameters)

# genus 1 bimodules: 

def init_ID1(torus_A):
    gen_by_name=AttrDict({
                "x": Generator("x"),
                "y": Generator("y")
                })
    gen_by_name.x.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.y.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)

    arrows=Bunch_of_arrows([
        (              gen_by_name.x,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r12,gen_by_name.x),

        (              gen_by_name.y,(torus_A.gen_by_name.r23,),
                torus_A.gen_by_name.r23,gen_by_name.y),

        (              gen_by_name.x,(torus_A.gen_by_name.r1,),
                torus_A.gen_by_name.r1,gen_by_name.y),

        (              gen_by_name.y,(torus_A.gen_by_name.r2,),
                torus_A.gen_by_name.r2,gen_by_name.x),

        (              gen_by_name.x,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r123,gen_by_name.y),

        (              gen_by_name.x,(torus_A.gen_by_name.r3,),
                torus_A.gen_by_name.r3,gen_by_name.y)])

    return DA_bimodule(gen_by_name,arrows,torus_A,torus_A,name="ID1")
def init_ID2(torus_A):
    gen_by_name=AttrDict({
                "x": Generator("x"),
                "y": Generator("y"),
                "z1": Generator("z1"),
                "z2": Generator("z2")
                })
    gen_by_name.x.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.y.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.z1.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.z2.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)

    arrows=Bunch_of_arrows([
        (       gen_by_name.x,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r12,gen_by_name.x
        ),

        (     gen_by_name.x,(torus_A.gen_by_name.r3,),
                1,gen_by_name.z2
        ),

        (       gen_by_name.x,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r12,gen_by_name.z2
        ),

        (                                       
            gen_by_name.x,(torus_A.gen_by_name.r1,),
            torus_A.gen_by_name.r1,gen_by_name.y
        ),

        (        gen_by_name.y,(torus_A.gen_by_name.r2,),
                torus_A.gen_by_name.r2,gen_by_name.x
        ),
        

        (        gen_by_name.y,(torus_A.gen_by_name.r23,),
                torus_A.gen_by_name.r2,gen_by_name.z2
        ),

        (        gen_by_name.z1,(),
                1,gen_by_name.z2
        ),

        (        gen_by_name.z1,(),
                torus_A.gen_by_name.r3,gen_by_name.y
        )

        ])

    return DA_bimodule(gen_by_name,arrows,torus_A,torus_A,name="ID2")
def init_ID3(torus_A):
    gen_by_name=AttrDict({
                "x": Generator("x"),
                "y": Generator("y"),
                "z1": Generator("z1"),
                "z2": Generator("z2"),
                "w1": Generator("w1"),
                "w2": Generator("w2")
                })
    gen_by_name.x.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.y.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.z1.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.z2.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.w1.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.w2.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)

    arrows=Bunch_of_arrows([
        (       gen_by_name.x,(),
                torus_A.gen_by_name.r1,gen_by_name.w2
        ),
        (     gen_by_name.x,(torus_A.gen_by_name.r3,),
                1,gen_by_name.z2
        ),



        (        gen_by_name.y,(torus_A.gen_by_name.r2,),
                torus_A.gen_by_name.r2,gen_by_name.x
        ),
        (        gen_by_name.y,(torus_A.gen_by_name.r23,),
                torus_A.gen_by_name.r2,gen_by_name.z2
        ),

        (        gen_by_name.z1,(),
                1,gen_by_name.z2
        ),
        (        gen_by_name.z1,(),
                torus_A.gen_by_name.r3,gen_by_name.y
        ),



        (        gen_by_name.w1,(),
                1,gen_by_name.w2
        ),
        (        gen_by_name.w1,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r2,gen_by_name.x
        ),
        (        gen_by_name.w1,(torus_A.gen_by_name.r1,),
                1,gen_by_name.y
        ),
        (        gen_by_name.w1,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r2,gen_by_name.z2
        )

        ])

    return DA_bimodule(gen_by_name,arrows,torus_A,torus_A,name="ID3")
def init_M_RHD(torus_A):
    gen_by_name=AttrDict({
                "p": Generator("p"),
                "q": Generator("q"),
                "r": Generator("r")
                })
    gen_by_name.p.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.q.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.r.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)

    arrows=Bunch_of_arrows([
        # from p to q
        (                         gen_by_name.p,(torus_A.gen_by_name.r1,),
                torus_A.gen_by_name.r1,gen_by_name.q),
        (                           gen_by_name.p,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r123,gen_by_name.q),
        (              gen_by_name.p,(torus_A.gen_by_name.r3,torus_A.gen_by_name.r23),
                torus_A.gen_by_name.r3,gen_by_name.q),

        # from p to r 
        (              gen_by_name.p,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r123,gen_by_name.r),
        (              gen_by_name.p,(torus_A.gen_by_name.r3,torus_A.gen_by_name.r2),
                torus_A.gen_by_name.r3,gen_by_name.r),

        #from r to p
        (              gen_by_name.r,(),
                torus_A.gen_by_name.r2,gen_by_name.p),

        #from r to q
        (              gen_by_name.r,(torus_A.gen_by_name.r3,),
                1,gen_by_name.q),

        #from q to r
        (              gen_by_name.q,(torus_A.gen_by_name.r2,),
                torus_A.gen_by_name.r23,gen_by_name.r),

        #from q to q
        (              gen_by_name.q,(torus_A.gen_by_name.r23,),
                torus_A.gen_by_name.r23,gen_by_name.q)
    ])

    return DA_bimodule(gen_by_name,arrows,torus_A,torus_A,name="M_RHD")
def init_M_LHD(torus_A):
    gen_by_name=AttrDict({
                "p": Generator("p"),
                "q": Generator("q"),
                "r": Generator("r")
                })
    gen_by_name.p.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.q.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.r.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)

    arrows=Bunch_of_arrows([
        # from p
        (                         gen_by_name.p,(torus_A.gen_by_name.r1,),
                torus_A.gen_by_name.r1,gen_by_name.q),
        (                           gen_by_name.p,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r123,gen_by_name.q),
        (              gen_by_name.p,(torus_A.gen_by_name.r123,torus_A.gen_by_name.r2),
                torus_A.gen_by_name.r12,gen_by_name.p),
        (              gen_by_name.p,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r1,gen_by_name.r),
        (              gen_by_name.p,(),
                torus_A.gen_by_name.r3,gen_by_name.r),

        #from r
        (              gen_by_name.r,(torus_A.gen_by_name.r3,torus_A.gen_by_name.r2),
                torus_A.gen_by_name.r2,gen_by_name.p),
        (              gen_by_name.r,(torus_A.gen_by_name.r3,),
                torus_A.gen_by_name.r23,gen_by_name.q),

        #from q
        (              gen_by_name.q,(torus_A.gen_by_name.r2,),
                1,gen_by_name.r),
        (              gen_by_name.q,(torus_A.gen_by_name.r23,),
                torus_A.gen_by_name.r23,gen_by_name.q),
        (              gen_by_name.q,(torus_A.gen_by_name.r23,torus_A.gen_by_name.r2),
                torus_A.gen_by_name.r2,gen_by_name.p),
    ])

    return DA_bimodule(gen_by_name,arrows,torus_A,torus_A,name="M_LHD")
def init_L_RHD(torus_A):
    gen_by_name=AttrDict({
                "p": Generator("p"),
                "q": Generator("q"),
                "s": Generator("s")
                })
    gen_by_name.p.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.q.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.s.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)

    arrows=Bunch_of_arrows([
        # from p to q
        (                         gen_by_name.q,(torus_A.gen_by_name.r2,torus_A.gen_by_name.r123),
                torus_A.gen_by_name.r23,gen_by_name.q),
        (                           gen_by_name.q,(torus_A.gen_by_name.r2,torus_A.gen_by_name.r12),
                torus_A.gen_by_name.r2,gen_by_name.p),
        (              gen_by_name.q,(torus_A.gen_by_name.r2,torus_A.gen_by_name.r1),
                torus_A.gen_by_name.r2,gen_by_name.s),

        # from p to r 
        (              gen_by_name.s,(torus_A.gen_by_name.r2,),
                1,gen_by_name.p),
        (              gen_by_name.s,(torus_A.gen_by_name.r23,),
                torus_A.gen_by_name.r3,gen_by_name.q),
        (              gen_by_name.s,(),
                torus_A.gen_by_name.r1,gen_by_name.q),

        #from r to p
        (              gen_by_name.p,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r123,gen_by_name.q),

        #from r to q
        (              gen_by_name.p,(torus_A.gen_by_name.r3,),
                torus_A.gen_by_name.r3,gen_by_name.q),

        #from q to r
        (              gen_by_name.p,(torus_A.gen_by_name.r1,),
                torus_A.gen_by_name.r12,gen_by_name.s),

        #from q to q
        (              gen_by_name.p,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r12,gen_by_name.p)
    ])

    return DA_bimodule(gen_by_name,arrows,torus_A,torus_A,name="L_RHD")
def init_L_LHD(torus_A):
    gen_by_name=AttrDict({
                "p": Generator("p"),
                "q": Generator("q"),
                "s": Generator("s")
                })
    gen_by_name.p.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.q.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.s.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)

    arrows=Bunch_of_arrows([
        # from q
        (                         gen_by_name.q,(),
                torus_A.gen_by_name.r2,gen_by_name.s),
        #from s
        (              gen_by_name.s,(torus_A.gen_by_name.r2,),
                torus_A.gen_by_name.r12,gen_by_name.p),
        (              gen_by_name.s,(torus_A.gen_by_name.r23,),
                torus_A.gen_by_name.r123,gen_by_name.q),
        (              gen_by_name.s,(torus_A.gen_by_name.r2,torus_A.gen_by_name.r1),
                torus_A.gen_by_name.r1,gen_by_name.q),

        #from p
        (              gen_by_name.p,(torus_A.gen_by_name.r1,),
                1,gen_by_name.s),
        (              gen_by_name.p,(torus_A.gen_by_name.r123,),
                torus_A.gen_by_name.r123,gen_by_name.q),
        (              gen_by_name.p,(torus_A.gen_by_name.r3,),
                torus_A.gen_by_name.r3,gen_by_name.q),
        (              gen_by_name.p,(torus_A.gen_by_name.r12,torus_A.gen_by_name.r1),
                torus_A.gen_by_name.r1,gen_by_name.q),
        (              gen_by_name.p,(torus_A.gen_by_name.r12,),
                torus_A.gen_by_name.r12,gen_by_name.p)
    ])

    return DA_bimodule(gen_by_name,arrows,torus_A,torus_A,name="L_LHD")

def init_AA_ID1(torus_A):
    gen_by_name=AttrDict({
                "y": Generator("y"),
                "x": Generator("x"),
                "z1": Generator("z1"),
                "z2": Generator("z2"),
                "w1": Generator("w1"),
                "w2": Generator("w2"),
                
                })

    gen_by_name.x.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.y.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.z1.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.z2.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.w1.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.w2.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)

    aa_arrows=Bunch_of_arrows([
        (       (torus_A.gen_by_name.r2,),gen_by_name.y,(torus_A.gen_by_name.r2,)
                        ,gen_by_name.x),

        (       (),gen_by_name.x,(torus_A.gen_by_name.r3,)
                        ,gen_by_name.z2),
        (       (torus_A.gen_by_name.r1,),gen_by_name.x,()
                        ,gen_by_name.w2),

        (       (torus_A.gen_by_name.r23,),gen_by_name.z1,(torus_A.gen_by_name.r2,)
                        ,gen_by_name.x),
        (       (torus_A.gen_by_name.r2,),gen_by_name.w1,(torus_A.gen_by_name.r12,)
                        ,gen_by_name.x),

        (      (torus_A.gen_by_name.r2,),gen_by_name.y,(torus_A.gen_by_name.r23,)
                        ,gen_by_name.z2),
        (       (torus_A.gen_by_name.r12,),gen_by_name.y,(torus_A.gen_by_name.r2,)
                        ,gen_by_name.w2),

        (       (torus_A.gen_by_name.r3,),gen_by_name.z1,()
                        ,gen_by_name.y),
        (       (),gen_by_name.w1,(torus_A.gen_by_name.r1,)
                        ,gen_by_name.y),

        (       (),gen_by_name.z1,()
                        ,gen_by_name.z2),
        (       (torus_A.gen_by_name.r23,),gen_by_name.z1,(torus_A.gen_by_name.r23,)
                        ,gen_by_name.z2),
        (       (),gen_by_name.w1,()
                        ,gen_by_name.w2),
        (       (torus_A.gen_by_name.r12,),gen_by_name.w1,(torus_A.gen_by_name.r12,)
                        ,gen_by_name.w2),

        (       (torus_A.gen_by_name.r2,),gen_by_name.w1,(torus_A.gen_by_name.r123,)
                        ,gen_by_name.z2),
        (       (torus_A.gen_by_name.r2,),gen_by_name.w1,(torus_A.gen_by_name.r3,torus_A.gen_by_name.r2,torus_A.gen_by_name.r1)
                        ,gen_by_name.z2),


        (       (torus_A.gen_by_name.r123,),gen_by_name.z1,(torus_A.gen_by_name.r2,)
                        ,gen_by_name.w2),
                                    ])
    return AA_bimodule(gen_by_name,aa_arrows,torus_A,torus_A,name="AA_id",to_check=True)
def init_DD_ID1(torus_A):
    gen_by_name=AttrDict({
                "p": Generator("p"),
                "q": Generator("q"),
                })

    gen_by_name.p.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.q.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)

    dd_arrows=Bunch_of_arrows([
        (                   gen_by_name.q,
        torus_A.gen_by_name.r2,gen_by_name.p,torus_A.gen_by_name.r2),

        (                   gen_by_name.p,
        torus_A.gen_by_name.r3,gen_by_name.q,torus_A.gen_by_name.r3),

        (                   gen_by_name.p,
        torus_A.gen_by_name.r1,gen_by_name.q,torus_A.gen_by_name.r1),
        (                   gen_by_name.p,
        torus_A.gen_by_name.r123,gen_by_name.q,torus_A.gen_by_name.r123)

                                    ])
    return DD_bimodule(gen_by_name,dd_arrows,torus_A,torus_A,name="DD_id")

def init_DD_TC(torus_A):
    gen_by_name=AttrDict({
                "afi": Generator("afi"),
                "afj": Generator("afj"),
                "afk": Generator("afk"),
                "afl": Generator("afl"),
                "ahi": Generator("ahi"),
                "ahj": Generator("ahj"),
                "ahk": Generator("ahk"),
                "ahl": Generator("ahl"),
                "ang": Generator("ang"),
                "amg": Generator("amg"),
                "ebi": Generator("ebi"),
                "ebj": Generator("ebj"),
                "ebk": Generator("ebk"),
                "ebl": Generator("ebl"),
                "edi": Generator("edi"),
                "edj": Generator("edj"),
                "edk": Generator("edk"),
                "edl": Generator("edl"),
                "enc": Generator("enc"),
                "emc": Generator("emc")
                })


    gen_by_name.afi.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.afj.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.afk.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.afl.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.ahi.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.ahj.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.ahk.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.ahl.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i1)
    gen_by_name.ang.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.amg.add_idems(torus_A.idem_by_name.i1,torus_A.idem_by_name.i0)
    gen_by_name.ebi.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.ebj.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.ebk.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.ebl.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.edi.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.edj.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.edk.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i0)
    gen_by_name.edl.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.enc.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)
    gen_by_name.emc.add_idems(torus_A.idem_by_name.i0,torus_A.idem_by_name.i1)

    dd_arrows=Bunch_of_arrows([
        (                   gen_by_name.ebl,
        torus_A.gen_by_name.r3,gen_by_name.afl,1),
        (                   gen_by_name.afl,
        torus_A.gen_by_name.r2,gen_by_name.edl,1),
        (                   gen_by_name.edl,
        torus_A.gen_by_name.r1,gen_by_name.ahl,1),

        (                   gen_by_name.ebk,
        torus_A.gen_by_name.r3,gen_by_name.afk,1),
        (                   gen_by_name.afk,
        torus_A.gen_by_name.r2,gen_by_name.edk,1),
        (                   gen_by_name.edk,
        torus_A.gen_by_name.r1,gen_by_name.ahk,1),

        (                   gen_by_name.ebj,
        torus_A.gen_by_name.r3,gen_by_name.afj,1),
        (                   gen_by_name.afj,
        torus_A.gen_by_name.r2,gen_by_name.edj,1),
        (                   gen_by_name.edj,
        torus_A.gen_by_name.r1,gen_by_name.ahj,1),

        (                   gen_by_name.ebi,
        torus_A.gen_by_name.r3,gen_by_name.afi,1),
        (                   gen_by_name.afi,
        torus_A.gen_by_name.r2,gen_by_name.edi,1),
        (                   gen_by_name.edi,
        torus_A.gen_by_name.r1,gen_by_name.ahi,1),

        (                   gen_by_name.ebl,
        1,gen_by_name.ebk,torus_A.gen_by_name.r1),
        (                   gen_by_name.ebk,
        1,gen_by_name.ebj,torus_A.gen_by_name.r2),
        (                   gen_by_name.ebj,
        1,gen_by_name.ebi,torus_A.gen_by_name.r3),

        (                   gen_by_name.afl,
        1,gen_by_name.afk,torus_A.gen_by_name.r1),
        (                   gen_by_name.afk,
        1,gen_by_name.afj,torus_A.gen_by_name.r2),
        (                   gen_by_name.afj,
        1,gen_by_name.afi,torus_A.gen_by_name.r3),

        (                   gen_by_name.edl,
        1,gen_by_name.edk,torus_A.gen_by_name.r1),
        (                   gen_by_name.edk,
        1,gen_by_name.edj,torus_A.gen_by_name.r2),
        (                   gen_by_name.edj,
        1,gen_by_name.edi,torus_A.gen_by_name.r3),

        (                   gen_by_name.ahl,
        1,gen_by_name.ahk,torus_A.gen_by_name.r1),
        (                   gen_by_name.ahk,
        1,gen_by_name.ahj,torus_A.gen_by_name.r2),
        (                   gen_by_name.ahj,
        1,gen_by_name.ahi,torus_A.gen_by_name.r3),

        (                   gen_by_name.ebl,
        1,gen_by_name.ebi,torus_A.gen_by_name.r123),
        (                   gen_by_name.afl,
        1,gen_by_name.afi,torus_A.gen_by_name.r123),
        (                   gen_by_name.edl,
        1,gen_by_name.edi,torus_A.gen_by_name.r123),
        (                   gen_by_name.ahl,
        1,gen_by_name.ahi,torus_A.gen_by_name.r123),


        # all connected to enc
        (                   gen_by_name.enc,
        torus_A.gen_by_name.r1,gen_by_name.ang,torus_A.gen_by_name.r1),
        (                   gen_by_name.enc,
        torus_A.gen_by_name.r123,gen_by_name.ang,torus_A.gen_by_name.r123),
        (                   gen_by_name.enc,
        torus_A.gen_by_name.r3,gen_by_name.ang,torus_A.gen_by_name.r3),
        (                   gen_by_name.ang,
        torus_A.gen_by_name.r2,gen_by_name.enc,torus_A.gen_by_name.r2),

        (                   gen_by_name.ebl,
        1,gen_by_name.enc,1),
        (                   gen_by_name.edl,
        1,gen_by_name.enc,torus_A.gen_by_name.r23),
        (                   gen_by_name.ebj,
        torus_A.gen_by_name.r12,gen_by_name.enc,1),

        # all connected to ang (ecxcept enc)

        (                   gen_by_name.afl,
        1,gen_by_name.ang,torus_A.gen_by_name.r3),
        (                   gen_by_name.ahl,
        1,gen_by_name.ang,torus_A.gen_by_name.r123),
        (                   gen_by_name.ebk,
        torus_A.gen_by_name.r1,gen_by_name.ang,1),
        (                   gen_by_name.ebi,
        torus_A.gen_by_name.r123,gen_by_name.ang,1),

        # all connected to emc

        (                   gen_by_name.emc,
        torus_A.gen_by_name.r1,gen_by_name.amg,torus_A.gen_by_name.r1),
        (                   gen_by_name.emc,
        torus_A.gen_by_name.r123,gen_by_name.amg,torus_A.gen_by_name.r123),
        (                   gen_by_name.emc,
        torus_A.gen_by_name.r3,gen_by_name.amg,torus_A.gen_by_name.r3),
        (                   gen_by_name.amg,
        torus_A.gen_by_name.r2,gen_by_name.emc,torus_A.gen_by_name.r2),

        (                   gen_by_name.emc,
        torus_A.gen_by_name.r3,gen_by_name.ahj,1),
        (                   gen_by_name.emc,
        torus_A.gen_by_name.r123,gen_by_name.ahl,1),
        (                   gen_by_name.emc,
        1,gen_by_name.edi,torus_A.gen_by_name.r1),
        (                   gen_by_name.emc,
        1,gen_by_name.ebi,torus_A.gen_by_name.r123),

        # all connected to amg (ecxcept emc)

        (                   gen_by_name.amg,
        torus_A.gen_by_name.r23,gen_by_name.ahk,1),
        (                   gen_by_name.amg,
        1,gen_by_name.ahi,1),
        (                   gen_by_name.amg,
        1,gen_by_name.afi,torus_A.gen_by_name.r12),

                                    ])

    return DD_bimodule(gen_by_name,dd_arrows,torus_A,torus_A,name="DD_TC")


ID1=init_ID1(torus_A) 
ID2=init_ID2(torus_A)
ID3=init_ID3(torus_A)
M_RHD=init_M_RHD(torus_A)
M_LHD=init_M_LHD(torus_A)
L_RHD=init_L_RHD(torus_A)
L_LHD=init_L_LHD(torus_A)
AA_id=init_AA_ID1(torus_A)
DD_id=init_DD_ID1(torus_A)
DD_TC=init_DD_TC(torus_A)

# genus 2 bimodules:

# this function reverses elements of algebra, useful when deriving say g2_K_LHD from g2_M_RHD
def g2_algebra_involution(a,g2_A):
    g2_algebra_involution_dict={
                        1:1,
                        g2_A.gen_by_name.r1:g2_A.gen_by_name.r7,
                        g2_A.gen_by_name.r2:g2_A.gen_by_name.r6,
                        g2_A.gen_by_name.r3:g2_A.gen_by_name.r5,
                        g2_A.gen_by_name.r4:g2_A.gen_by_name.r4,
                        g2_A.gen_by_name.r5:g2_A.gen_by_name.r3,
                        g2_A.gen_by_name.r6:g2_A.gen_by_name.r2,
                        g2_A.gen_by_name.r7:g2_A.gen_by_name.r1,
                        g2_A.gen_by_name.r12:g2_A.gen_by_name.r67,
                        g2_A.gen_by_name.r23:g2_A.gen_by_name.r56,
                        g2_A.gen_by_name.r34:g2_A.gen_by_name.r45,
                        g2_A.gen_by_name.r45:g2_A.gen_by_name.r34,
                        g2_A.gen_by_name.r56:g2_A.gen_by_name.r23,
                        g2_A.gen_by_name.r67:g2_A.gen_by_name.r12,
                        g2_A.gen_by_name.r123:g2_A.gen_by_name.r567,
                        g2_A.gen_by_name.r234:g2_A.gen_by_name.r456,
                        g2_A.gen_by_name.r345:g2_A.gen_by_name.r345,
                        g2_A.gen_by_name.r456:g2_A.gen_by_name.r234,
                        g2_A.gen_by_name.r567:g2_A.gen_by_name.r123,
                        g2_A.gen_by_name.r1234:g2_A.gen_by_name.r4567,
                        g2_A.gen_by_name.r2345:g2_A.gen_by_name.r3456,
                        g2_A.gen_by_name.r3456:g2_A.gen_by_name.r2345,
                        g2_A.gen_by_name.r4567:g2_A.gen_by_name.r1234,
                        g2_A.gen_by_name.r12345:g2_A.gen_by_name.r34567,
                        g2_A.gen_by_name.r23456:g2_A.gen_by_name.r23456,
                        g2_A.gen_by_name.r34567:g2_A.gen_by_name.r12345,
                        g2_A.gen_by_name.r123456:g2_A.gen_by_name.r234567,
                        g2_A.gen_by_name.r234567:g2_A.gen_by_name.r123456,
                        g2_A.gen_by_name.r1234567:g2_A.gen_by_name.r1234567,
                        }
    try: return g2_algebra_involution_dict[a]
    except: print (a)

def init_g2_ID(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)

    #this is for convenience
    idem_to_gens={g2_A.idem_by_name.i0:gen_by_name.x0,
                    g2_A.idem_by_name.i1:gen_by_name.x1,
                    g2_A.idem_by_name.i2:gen_by_name.x2,
                    g2_A.idem_by_name.i3:gen_by_name.x3}

    #initializing all actions
    arrows=Bunch_of_arrows()
    for algebra_element in (set(g2_A.genset) - set(g2_A.idemset)):
        gen_in=idem_to_gens[algebra_element.idem.left]
        gen_out=idem_to_gens[algebra_element.idem.right]
        
        arrows[(gen_in,(algebra_element,),
                algebra_element,gen_out)]+=1

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_ID")
def init_g2_ID_bounded(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "z1": Generator("z1"),
                "z2": Generator("z2"),
                "w1": Generator("w1"),
                "w2": Generator("w2"),
                "t1": Generator("t1"),
                "t2": Generator("t2"),
                "c1": Generator("c1"),
                "c2": Generator("c2")

                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.z1.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i1)
    gen_by_name.z2.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i1)
    gen_by_name.w1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i0)
    gen_by_name.w2.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i0)
    gen_by_name.t1.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i3)
    gen_by_name.t2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i3)
    gen_by_name.c1.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i2)
    gen_by_name.c2.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i2)

    arrows=Bunch_of_arrows([
        # higher torus
        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),
        (              gen_by_name.t1,(),
                1,gen_by_name.t2),
        (              gen_by_name.t1,(),
                g2_A.gen_by_name.r7,gen_by_name.x3),
        (               gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r6,gen_by_name.t2),

        (               gen_by_name.x2,(),
                g2_A.gen_by_name.r5,gen_by_name.c2),

        (               gen_by_name.x2,(g2_A.gen_by_name.r7,),
                1,gen_by_name.t2),

        (               gen_by_name.c1,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (               gen_by_name.c1,(),
                    1,gen_by_name.c2),
        (               gen_by_name.c1,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r6,gen_by_name.t2),
        (               gen_by_name.c1,(g2_A.gen_by_name.r5,),
                1,gen_by_name.x3),
        
        # higher torus
        # (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
        #         g2_A.gen_by_name.r6,gen_by_name.x2),
        # (               gen_by_name.t1,(),
        #         g2_A.gen_by_name.r7,gen_by_name.x3),
        # (               gen_by_name.t1,(),
        #         1,gen_by_name.t2),
        # (               gen_by_name.x3,(g2_A.gen_by_name.r67,),
        #         g2_A.gen_by_name.r6,gen_by_name.t2),
        # (               gen_by_name.x2,(g2_A.gen_by_name.r567,),
        #         g2_A.gen_by_name.r56,gen_by_name.t2),
        # (               gen_by_name.x2,(g2_A.gen_by_name.r7,),
        #         1,gen_by_name.t2),

        #lower torus
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,),
                g2_A.gen_by_name.r2,gen_by_name.x0),
        (              gen_by_name.z1,(),
                1,gen_by_name.z2),
        (              gen_by_name.z1,(),
                g2_A.gen_by_name.r3,gen_by_name.x1),
        (               gen_by_name.x1,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r2,gen_by_name.z2),

        (               gen_by_name.x0,(),
                g2_A.gen_by_name.r1,gen_by_name.w2),
        (               gen_by_name.x0,(g2_A.gen_by_name.r3,),
                1,gen_by_name.z2),
        (               gen_by_name.w1,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r2,gen_by_name.x0),
        (               gen_by_name.w1,(),
                    1,gen_by_name.w2),
        (               gen_by_name.w1,(g2_A.gen_by_name.r1,),
                1,gen_by_name.x1),
        (               gen_by_name.w1,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r2,gen_by_name.z2),


        # from x0
        (              gen_by_name.x0,(g2_A.gen_by_name.r34,g2_A.gen_by_name.r7),
                g2_A.gen_by_name.r34,gen_by_name.t1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r3456,gen_by_name.t2),

        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r345,gen_by_name.c1),

            (              gen_by_name.x0,(g2_A.gen_by_name.r345,),
                    g2_A.gen_by_name.r345,gen_by_name.x3),

        #from x1

            (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                    g2_A.gen_by_name.r45,gen_by_name.x3),
            (              gen_by_name.x1,(g2_A.gen_by_name.r2345,),
                    g2_A.gen_by_name.r2345,gen_by_name.x3),
        
        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r2345,gen_by_name.c1),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,g2_A.gen_by_name.r7),
                g2_A.gen_by_name.r4,gen_by_name.t1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234,g2_A.gen_by_name.r7),
                g2_A.gen_by_name.r234,gen_by_name.t1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r456,gen_by_name.t2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r23456,gen_by_name.t2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r234,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r45,gen_by_name.c1),

        # added after perturbation z1 and z2, w1 and w2, t1 and t2

        (              gen_by_name.z2,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.z2,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),
            (              gen_by_name.z2,(g2_A.gen_by_name.r45,),
                    g2_A.gen_by_name.r345,gen_by_name.x3),
        

        (               gen_by_name.w1,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r234,gen_by_name.x2),
        (               gen_by_name.w1,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),
            (               gen_by_name.w1,(g2_A.gen_by_name.r12345,),
                    g2_A.gen_by_name.r2345,gen_by_name.x3),

        (               gen_by_name.w1,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r2345,gen_by_name.c1),    
        (               gen_by_name.w1,(g2_A.gen_by_name.r1234,g2_A.gen_by_name.r7),
                g2_A.gen_by_name.r234,gen_by_name.t1),
        (               gen_by_name.w1,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r23456,gen_by_name.t2),


        (               gen_by_name.w2,(g2_A.gen_by_name.r1,g2_A.gen_by_name.r4),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (               gen_by_name.w2,(g2_A.gen_by_name.r1,g2_A.gen_by_name.r4),
                g2_A.gen_by_name.r45,gen_by_name.c1),
        (               gen_by_name.w2,(g2_A.gen_by_name.r1,g2_A.gen_by_name.r456),
                g2_A.gen_by_name.r456,gen_by_name.x2),

        (               gen_by_name.w2,(g2_A.gen_by_name.r1,g2_A.gen_by_name.r4,g2_A.gen_by_name.r7),
                g2_A.gen_by_name.r4,gen_by_name.t1),
        (               gen_by_name.w2,(g2_A.gen_by_name.r1,g2_A.gen_by_name.r4567),
                g2_A.gen_by_name.r456,gen_by_name.t2),


            (               gen_by_name.w2,(g2_A.gen_by_name.r1,g2_A.gen_by_name.r45),
                    g2_A.gen_by_name.r45,gen_by_name.x3),


        (              gen_by_name.z2,(g2_A.gen_by_name.r4,g2_A.gen_by_name.r7),
                g2_A.gen_by_name.r34,gen_by_name.t1),
        (              gen_by_name.z2,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r345,gen_by_name.c1),
        (              gen_by_name.z2,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r3456,gen_by_name.t2),


        



        ])
    

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_ID_bounded")

def init_g2_M_RHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "r":Generator("r")
                })
    
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.r.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i0)

    arrows=Bunch_of_arrows([

        #short near-chords

        (                         gen_by_name.x0,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r1,gen_by_name.x1),

        (              gen_by_name.r,(),
                g2_A.gen_by_name.r2,gen_by_name.x0),

        (              gen_by_name.r,(g2_A.gen_by_name.r3,),
                1,gen_by_name.x1),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        #two extra

        (              gen_by_name.x1,(g2_A.gen_by_name.r2,),
                g2_A.gen_by_name.r23,gen_by_name.r),

        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r3,gen_by_name.r),
        

        # other

        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r23),
                g2_A.gen_by_name.r3,gen_by_name.x1),

        (              gen_by_name.r,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r4,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x1,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r23,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),

        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r123,gen_by_name.r),
        

        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r2345),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r234567),
                g2_A.gen_by_name.r34567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r234),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r23456),
                g2_A.gen_by_name.r3456,gen_by_name.x2),
        
        (              gen_by_name.r,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),
        (              gen_by_name.r,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.r,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),


        
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r234,gen_by_name.x2)



        ])

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_M_RHD")
# next bimodule we derive from the previous one by horizontal reflection
def init_g2_K_LHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "r":Generator("r")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.r.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i3)

    arrows_old=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        # from x0
        (                         gen_by_name.x0,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r23),
                g2_A.gen_by_name.r3,gen_by_name.x1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r123,gen_by_name.r),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r3,gen_by_name.r),

        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r2345),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r234567),
                g2_A.gen_by_name.r34567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r234),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r23456),
                g2_A.gen_by_name.r3456,gen_by_name.x2),



        #from r
        (              gen_by_name.r,(),
                g2_A.gen_by_name.r2,gen_by_name.x0),
        (              gen_by_name.r,(g2_A.gen_by_name.r3,),
                1,gen_by_name.x1),
        (              gen_by_name.r,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.r,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),
        (              gen_by_name.r,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.r,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),

        #from x1
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,),
                g2_A.gen_by_name.r23,gen_by_name.r),
        (              gen_by_name.x1,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r23,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r234,gen_by_name.x2)



        ])
    

    arrows=Bunch_of_arrows()
    for old_arrow in arrows_old:
        in_gen=old_arrow[3]
        out_gen=old_arrow[0]

        in_alg_tuple=[]
        for a in reversed(old_arrow[1]):
            in_alg_tuple.append(g2_algebra_involution(a,g2_A))
        in_alg_tuple=tuple(in_alg_tuple)
        out_alg=g2_algebra_involution(old_arrow[2],g2_A)
        
        new_arrow=(              in_gen,in_alg_tuple,
                out_alg,out_gen)

        arrows[new_arrow]+=1

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_K_LHD")

def init_g2_M_LHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "r":Generator("r")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.r.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i0)

    arrows=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        # from x0
        (                         gen_by_name.x0,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (              gen_by_name.x0,(),
                g2_A.gen_by_name.r3,gen_by_name.r),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r1,gen_by_name.r),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),


        #from x1
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,),
                1,gen_by_name.r),
        (              gen_by_name.x1,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r23,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r23,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r2,gen_by_name.x0),

        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r234,gen_by_name.x2),

        #from r
        (              gen_by_name.r,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r2,gen_by_name.x0),
        (              gen_by_name.r,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r23,gen_by_name.x1),

        (              gen_by_name.r,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r234,gen_by_name.x2),
        (              gen_by_name.r,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        (              gen_by_name.r,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.r,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3)



        ])

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_M_LHD")
# next bimodule we derive from the previous one by horizontal reflection
def init_g2_K_RHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "r":Generator("r")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.r.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i3)

    arrows_old=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        # from x0
        (                         gen_by_name.x0,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (              gen_by_name.x0,(),
                g2_A.gen_by_name.r3,gen_by_name.r),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r1,gen_by_name.r),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),


        #from x1
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,),
                1,gen_by_name.r),
        (              gen_by_name.x1,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r23,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r23,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r2,gen_by_name.x0),

        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r234,gen_by_name.x2),

        #from r
        (              gen_by_name.r,(g2_A.gen_by_name.r3,g2_A.gen_by_name.r2),
                g2_A.gen_by_name.r2,gen_by_name.x0),
        (              gen_by_name.r,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r23,gen_by_name.x1),

        (              gen_by_name.r,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r234,gen_by_name.x2),
        (              gen_by_name.r,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        (              gen_by_name.r,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.r,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3)



        ])

    arrows=Bunch_of_arrows()
    for old_arrow in arrows_old:
        in_gen=old_arrow[3]
        out_gen=old_arrow[0]

        in_alg_tuple=[]
        for a in reversed(old_arrow[1]):
            in_alg_tuple.append(g2_algebra_involution(a,g2_A))
        in_alg_tuple=tuple(in_alg_tuple)
        out_alg=g2_algebra_involution(old_arrow[2],g2_A)
        
        new_arrow=(              in_gen,in_alg_tuple,
                out_alg,out_gen)

        arrows[new_arrow]+=1


    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_K_RHD")

def init_g2_L_RHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "s":Generator("s")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.s.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i1)

    arrows=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        # from x0
        (                         gen_by_name.x0,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r3,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r12,gen_by_name.s),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r34567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),


        #from x1
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1),
                g2_A.gen_by_name.r2,gen_by_name.s),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r123),
                g2_A.gen_by_name.r23,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r12),
                g2_A.gen_by_name.r2,gen_by_name.x0),

        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r12345),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1234567),
                g2_A.gen_by_name.r234567,gen_by_name.x3),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1234),
                g2_A.gen_by_name.r234,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r123456),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        #from s
        (              gen_by_name.s,(g2_A.gen_by_name.r2,),
                1,gen_by_name.x0),
        (              gen_by_name.s,(),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (              gen_by_name.s,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r3,gen_by_name.x1),

        (              gen_by_name.s,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.s,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),

        (              gen_by_name.s,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.s,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r34567,gen_by_name.x3)



        ])

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_L_RHD")
# next bimodule we derive from the previous one by horizontal reflection
def init_g2_N_LHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "s":Generator("s")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.s.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i2)

    arrows_old=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        # from x0
        (                         gen_by_name.x0,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r3,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r12,gen_by_name.s),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r34567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),


        #from x1
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1),
                g2_A.gen_by_name.r2,gen_by_name.s),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r123),
                g2_A.gen_by_name.r23,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r12),
                g2_A.gen_by_name.r2,gen_by_name.x0),

        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r12345),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1234567),
                g2_A.gen_by_name.r234567,gen_by_name.x3),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1234),
                g2_A.gen_by_name.r234,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r123456),
                g2_A.gen_by_name.r23456,gen_by_name.x2),

        #from s
        (              gen_by_name.s,(g2_A.gen_by_name.r2,),
                1,gen_by_name.x0),
        (              gen_by_name.s,(),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (              gen_by_name.s,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r3,gen_by_name.x1),

        (              gen_by_name.s,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.s,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),

        (              gen_by_name.s,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.s,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r34567,gen_by_name.x3)



        ])

    arrows=Bunch_of_arrows()
    for old_arrow in arrows_old:
        in_gen=old_arrow[3]
        out_gen=old_arrow[0]

        in_alg_tuple=[]
        for a in reversed(old_arrow[1]):
            in_alg_tuple.append(g2_algebra_involution(a,g2_A))
        in_alg_tuple=tuple(in_alg_tuple)
        out_alg=g2_algebra_involution(old_arrow[2],g2_A)
        
        new_arrow=(              in_gen,in_alg_tuple,
                out_alg,out_gen)

        arrows[new_arrow]+=1

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_N_LHD")

def init_g2_L_LHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "s":Generator("s")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.s.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i1)

    arrows=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        # from x0
        (                         gen_by_name.x0,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r3,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r12,g2_A.gen_by_name.r1),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1,),
                1,gen_by_name.s),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r34567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),


        #from x1
        (              gen_by_name.x1,(),
                g2_A.gen_by_name.r2,gen_by_name.s),

        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),

        #from s
        (              gen_by_name.s,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (              gen_by_name.s,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (              gen_by_name.s,(g2_A.gen_by_name.r2,),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.s,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.s,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),

        (              gen_by_name.s,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.s,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3)



        ])

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_L_LHD")
# next bimodule we derive from the previous one by horizontal reflection
def init_g2_N_RHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "s":Generator("s")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.s.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i2)

    arrows_old=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),

        (              gen_by_name.x3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r67,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),

        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),

        (              gen_by_name.x2,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r567,gen_by_name.x3),

        (              gen_by_name.x2,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x3),

        # from x0
        (                         gen_by_name.x0,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r3,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (                           gen_by_name.x0,(g2_A.gen_by_name.r12,g2_A.gen_by_name.r1),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1,),
                1,gen_by_name.s),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r34567,gen_by_name.x3),

        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),


        #from x1
        (              gen_by_name.x1,(),
                g2_A.gen_by_name.r2,gen_by_name.s),

        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r4567,gen_by_name.x3),

        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),

        #from s
        (              gen_by_name.s,(g2_A.gen_by_name.r2,g2_A.gen_by_name.r1),
                g2_A.gen_by_name.r1,gen_by_name.x1),
        (              gen_by_name.s,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r123,gen_by_name.x1),
        (              gen_by_name.s,(g2_A.gen_by_name.r2,),
                g2_A.gen_by_name.r12,gen_by_name.x0),

        (              gen_by_name.s,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r1234,gen_by_name.x2),
        (              gen_by_name.s,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r123456,gen_by_name.x2),

        (              gen_by_name.s,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r12345,gen_by_name.x3),
        (              gen_by_name.s,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3)



        ])

    arrows=Bunch_of_arrows()
    for old_arrow in arrows_old:
        in_gen=old_arrow[3]
        out_gen=old_arrow[0]

        in_alg_tuple=[]
        for a in reversed(old_arrow[1]):
            in_alg_tuple.append(g2_algebra_involution(a,g2_A))
        in_alg_tuple=tuple(in_alg_tuple)
        out_alg=g2_algebra_involution(old_arrow[2],g2_A)
        
        new_arrow=(              in_gen,in_alg_tuple,
                out_alg,out_gen)

        arrows[new_arrow]+=1

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_N_RHD")

def init_g2_T_LHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "t1":Generator("t1"),
                "t2":Generator("t2"),
                "t3":Generator("t3"),
                "t4":Generator("t4"),
                "t5":Generator("t5"),
                "t6":Generator("t6"),
                "t7":Generator("t7"),
                "t8":Generator("t8"),
                "t9":Generator("t9"),
                "t10":Generator("t10"),
                "t11":Generator("t11"),
                "t12":Generator("t12")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)

    gen_by_name.t1.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.t2.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i1)
    gen_by_name.t3.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i0)
    gen_by_name.t4.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i0)

    gen_by_name.t5.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.t6.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i1)
    gen_by_name.t7.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i2)
    gen_by_name.t8.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)

    gen_by_name.t9.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i3)
    gen_by_name.t10.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i3)
    gen_by_name.t11.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i2)
    gen_by_name.t12.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)

    arrows=Bunch_of_arrows([
        (              gen_by_name.x2,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x3),
        (              gen_by_name.x3,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x2),
        (              gen_by_name.x2,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r3,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r2,),
                g2_A.gen_by_name.r2,gen_by_name.x0),
        
        (              gen_by_name.x1,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r23,gen_by_name.x1),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r34,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r345,gen_by_name.x3),
        
        (              gen_by_name.x1,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r2345,gen_by_name.x3),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r234,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r23456,gen_by_name.x2),
        (              gen_by_name.x0,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r3456,gen_by_name.x2),

        # between odd and even t's
        (              gen_by_name.t1,(),
                g2_A.gen_by_name.r4,gen_by_name.t2),
        (              gen_by_name.t3,(),
                g2_A.gen_by_name.r4,gen_by_name.t4),
        (              gen_by_name.t5,(),
                g2_A.gen_by_name.r4,gen_by_name.t6),
        (              gen_by_name.t7,(),
                g2_A.gen_by_name.r4,gen_by_name.t8),
        (              gen_by_name.t9,(),
                g2_A.gen_by_name.r4,gen_by_name.t10),
        (              gen_by_name.t11,(),
                g2_A.gen_by_name.r4,gen_by_name.t12),

        # between odd t's
        (              gen_by_name.t1,(g2_A.gen_by_name.r2,),
                1,gen_by_name.t3),
        (              gen_by_name.t3,(g2_A.gen_by_name.r3,),
                1,gen_by_name.t5),
        (              gen_by_name.t5,(g2_A.gen_by_name.r4,),
                1,gen_by_name.t7),
        (              gen_by_name.t7,(g2_A.gen_by_name.r5,),
                1,gen_by_name.t9),
        (              gen_by_name.t9,(g2_A.gen_by_name.r6,),
                1,gen_by_name.t11),
        
        (              gen_by_name.t1,(g2_A.gen_by_name.r23,),
                1,gen_by_name.t5),
        (              gen_by_name.t3,(g2_A.gen_by_name.r34,),
                1,gen_by_name.t7),
        (              gen_by_name.t5,(g2_A.gen_by_name.r45,),
                1,gen_by_name.t9),
        (              gen_by_name.t7,(g2_A.gen_by_name.r56,),
                1,gen_by_name.t11),
        (              gen_by_name.t1,(g2_A.gen_by_name.r234,),
                1,gen_by_name.t7),
        
        (              gen_by_name.t3,(g2_A.gen_by_name.r345,),
                1,gen_by_name.t9),
        (              gen_by_name.t5,(g2_A.gen_by_name.r456,),
                1,gen_by_name.t11),
        (              gen_by_name.t1,(g2_A.gen_by_name.r2345,),
                1,gen_by_name.t9),
        (              gen_by_name.t3,(g2_A.gen_by_name.r3456,),
                1,gen_by_name.t11),
        (              gen_by_name.t1,(g2_A.gen_by_name.r23456,),
                1,gen_by_name.t11),

        # between even t's
        (              gen_by_name.t2,(g2_A.gen_by_name.r2,),
                1,gen_by_name.t4),
        (              gen_by_name.t4,(g2_A.gen_by_name.r3,),
                1,gen_by_name.t6),
        (              gen_by_name.t6,(g2_A.gen_by_name.r4,),
                1,gen_by_name.t8),
        (              gen_by_name.t8,(g2_A.gen_by_name.r5,),
                1,gen_by_name.t10),
        (              gen_by_name.t10,(g2_A.gen_by_name.r6,),
                1,gen_by_name.t12),
        
        (              gen_by_name.t2,(g2_A.gen_by_name.r23,),
                1,gen_by_name.t6),
        (              gen_by_name.t4,(g2_A.gen_by_name.r34,),
                1,gen_by_name.t8),
        (              gen_by_name.t6,(g2_A.gen_by_name.r45,),
                1,gen_by_name.t10),
        (              gen_by_name.t8,(g2_A.gen_by_name.r56,),
                1,gen_by_name.t12),
        (              gen_by_name.t2,(g2_A.gen_by_name.r234,),
                1,gen_by_name.t8),
        
        (              gen_by_name.t4,(g2_A.gen_by_name.r345,),
                1,gen_by_name.t10),
        (              gen_by_name.t6,(g2_A.gen_by_name.r456,),
                1,gen_by_name.t12),
        (              gen_by_name.t2,(g2_A.gen_by_name.r2345,),
                1,gen_by_name.t10),
        (              gen_by_name.t4,(g2_A.gen_by_name.r3456,),
                1,gen_by_name.t12),
        (              gen_by_name.t2,(g2_A.gen_by_name.r23456,),
                1,gen_by_name.t12),



        #actual non-trivial stuff

        # ??
        (              gen_by_name.x0,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r1,gen_by_name.t1),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r1,gen_by_name.t3),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r1,gen_by_name.t5),
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r1,gen_by_name.t7),
        (              gen_by_name.x0,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r1,gen_by_name.t9),
        (              gen_by_name.x0,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r1,gen_by_name.t11),
        # ??

        # ??
        (              gen_by_name.t11,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        (              gen_by_name.t9,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        (              gen_by_name.t7,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        (              gen_by_name.t5,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        (              gen_by_name.t3,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        (              gen_by_name.t1,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r234567,gen_by_name.x3),
        # ??

        # ??
        (              gen_by_name.x0,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x3),
        (              gen_by_name.t2,(g2_A.gen_by_name.r4,g2_A.gen_by_name.r7),
                g2_A.gen_by_name.r7,gen_by_name.x3),
        # ??

        (              gen_by_name.x2,(),
                1,gen_by_name.t12),
        (              gen_by_name.x3,(),
                g2_A.gen_by_name.r6,gen_by_name.t10),
        (              gen_by_name.x1,(g2_A.gen_by_name.r4,),
                1,gen_by_name.t11),
        (              gen_by_name.x1,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r23,gen_by_name.t11),

        (              gen_by_name.x0,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r3,gen_by_name.t11),
        (              gen_by_name.x2,(),
                g2_A.gen_by_name.r56,gen_by_name.t8),
        (              gen_by_name.x1,(),
                g2_A.gen_by_name.r456,gen_by_name.t6),
        (              gen_by_name.x1,(),
                g2_A.gen_by_name.r23456,gen_by_name.t2),
        (              gen_by_name.x0,(),
                g2_A.gen_by_name.r3456,gen_by_name.t4),


        ])

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_T_LHD")
# next bimodule we derive from the previous one by horizontal reflection
def init_g2_T_RHD(g2_A):
    gen_by_name=AttrDict({
                "x0": Generator("x0"),
                "x1": Generator("x1"),
                "x2": Generator("x2"),
                "x3": Generator("x3"),
                "t1":Generator("t1"),
                "t2":Generator("t2"),
                "t3":Generator("t3"),
                "t4":Generator("t4"),
                "t5":Generator("t5"),
                "t6":Generator("t6"),
                "t7":Generator("t7"),
                "t8":Generator("t8"),
                "t9":Generator("t9"),
                "t10":Generator("t10"),
                "t11":Generator("t11"),
                "t12":Generator("t12")
                })
    gen_by_name.x0.add_idems(g2_A.idem_by_name.i3,g2_A.idem_by_name.i3)
    gen_by_name.x1.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.x2.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    gen_by_name.x3.add_idems(g2_A.idem_by_name.i0,g2_A.idem_by_name.i0)

    gen_by_name.t1.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.t2.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i2)
    gen_by_name.t3.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i3)
    gen_by_name.t4.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i3)

    gen_by_name.t5.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i2)
    gen_by_name.t6.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i2)
    gen_by_name.t7.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i1)
    gen_by_name.t8.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)

    gen_by_name.t9.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i0)
    gen_by_name.t10.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i0)
    gen_by_name.t11.add_idems(g2_A.idem_by_name.i2,g2_A.idem_by_name.i1)
    gen_by_name.t12.add_idems(g2_A.idem_by_name.i1,g2_A.idem_by_name.i1)
    
    arrows=Bunch_of_arrows([
        (              gen_by_name.x3,(g2_A.gen_by_name.r3,),
                g2_A.gen_by_name.r3,gen_by_name.x2),
        (              gen_by_name.x2,(g2_A.gen_by_name.r2,),
                g2_A.gen_by_name.r2,gen_by_name.x3),
        (              gen_by_name.x2,(g2_A.gen_by_name.r23,),
                g2_A.gen_by_name.r23,gen_by_name.x2),
        (              gen_by_name.x1,(g2_A.gen_by_name.r5,),
                g2_A.gen_by_name.r5,gen_by_name.x0),
        (              gen_by_name.x0,(g2_A.gen_by_name.r6,),
                g2_A.gen_by_name.r6,gen_by_name.x1),
        
        (              gen_by_name.x1,(g2_A.gen_by_name.r56,),
                g2_A.gen_by_name.r56,gen_by_name.x1),
        (              gen_by_name.x2,(g2_A.gen_by_name.r4,),
                g2_A.gen_by_name.r4,gen_by_name.x1),
        (              gen_by_name.x2,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r45,gen_by_name.x0),
        (              gen_by_name.x3,(g2_A.gen_by_name.r34,),
                g2_A.gen_by_name.r34,gen_by_name.x1),
        (              gen_by_name.x3,(g2_A.gen_by_name.r345,),
                g2_A.gen_by_name.r345,gen_by_name.x0),
        
        (              gen_by_name.x3,(g2_A.gen_by_name.r3456,),
                g2_A.gen_by_name.r3456,gen_by_name.x1),
        (              gen_by_name.x2,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r456,gen_by_name.x1),
        (              gen_by_name.x2,(g2_A.gen_by_name.r234,),
                g2_A.gen_by_name.r234,gen_by_name.x1),
        (              gen_by_name.x2,(g2_A.gen_by_name.r23456,),
                g2_A.gen_by_name.r23456,gen_by_name.x1),
        (              gen_by_name.x2,(g2_A.gen_by_name.r2345,),
                g2_A.gen_by_name.r2345,gen_by_name.x0),

        # between odd and even t's
        (              gen_by_name.t2,(),
                g2_A.gen_by_name.r4,gen_by_name.t1),
        (              gen_by_name.t4,(),
                g2_A.gen_by_name.r4,gen_by_name.t3),
        (              gen_by_name.t6,(),
                g2_A.gen_by_name.r4,gen_by_name.t5),
        (              gen_by_name.t8,(),
                g2_A.gen_by_name.r4,gen_by_name.t7),
        (              gen_by_name.t10,(),
                g2_A.gen_by_name.r4,gen_by_name.t9),
        (              gen_by_name.t12,(),
                g2_A.gen_by_name.r4,gen_by_name.t11),

        # between odd t's
        (              gen_by_name.t3,(g2_A.gen_by_name.r6,),
                1,gen_by_name.t1),
        (              gen_by_name.t5,(g2_A.gen_by_name.r5,),
                1,gen_by_name.t3),
        (              gen_by_name.t7,(g2_A.gen_by_name.r4,),
                1,gen_by_name.t5),
        (              gen_by_name.t9,(g2_A.gen_by_name.r3,),
                1,gen_by_name.t7),
        (              gen_by_name.t11,(g2_A.gen_by_name.r2,),
                1,gen_by_name.t9),
        
        (              gen_by_name.t5,(g2_A.gen_by_name.r56,),
                1,gen_by_name.t1),
        (              gen_by_name.t7,(g2_A.gen_by_name.r45,),
                1,gen_by_name.t3),
        (              gen_by_name.t9,(g2_A.gen_by_name.r34,),
                1,gen_by_name.t5),
        (              gen_by_name.t11,(g2_A.gen_by_name.r23,),
                1,gen_by_name.t7),
        (              gen_by_name.t7,(g2_A.gen_by_name.r456,),
                1,gen_by_name.t1),
        
        (              gen_by_name.t9,(g2_A.gen_by_name.r345,),
                1,gen_by_name.t3),
        (              gen_by_name.t11,(g2_A.gen_by_name.r234,),
                1,gen_by_name.t5),
        (              gen_by_name.t9,(g2_A.gen_by_name.r3456,),
                1,gen_by_name.t1),
        (              gen_by_name.t11,(g2_A.gen_by_name.r2345,),
                1,gen_by_name.t3),
        (              gen_by_name.t11,(g2_A.gen_by_name.r23456,),
                1,gen_by_name.t1),

        # between even t's
        (              gen_by_name.t4,(g2_A.gen_by_name.r6,),
                1,gen_by_name.t2),
        (              gen_by_name.t6,(g2_A.gen_by_name.r5,),
                1,gen_by_name.t4),
        (              gen_by_name.t8,(g2_A.gen_by_name.r4,),
                1,gen_by_name.t6),
        (              gen_by_name.t10,(g2_A.gen_by_name.r3,),
                1,gen_by_name.t8),
        (              gen_by_name.t12,(g2_A.gen_by_name.r2,),
                1,gen_by_name.t10),
        
        (              gen_by_name.t6,(g2_A.gen_by_name.r56,),
                1,gen_by_name.t2),
        (              gen_by_name.t8,(g2_A.gen_by_name.r45,),
                1,gen_by_name.t4),
        (              gen_by_name.t10,(g2_A.gen_by_name.r34,),
                1,gen_by_name.t6),
        (              gen_by_name.t12,(g2_A.gen_by_name.r23,),
                1,gen_by_name.t8),
        (              gen_by_name.t8,(g2_A.gen_by_name.r456,),
                1,gen_by_name.t2),
        
        (              gen_by_name.t10,(g2_A.gen_by_name.r345,),
                1,gen_by_name.t4),
        (              gen_by_name.t12,(g2_A.gen_by_name.r234,),
                1,gen_by_name.t6),
        (              gen_by_name.t10,(g2_A.gen_by_name.r3456,),
                1,gen_by_name.t2),
        (              gen_by_name.t12,(g2_A.gen_by_name.r2345,),
                1,gen_by_name.t4),
        (              gen_by_name.t12,(g2_A.gen_by_name.r23456,),
                1,gen_by_name.t2),



        #actual non-trivial stuff
        # ??
        (              gen_by_name.t1,(g2_A.gen_by_name.r7,),
                g2_A.gen_by_name.r7,gen_by_name.x0),
        (              gen_by_name.t3,(g2_A.gen_by_name.r67,),
                g2_A.gen_by_name.r7,gen_by_name.x0),
        (              gen_by_name.t5,(g2_A.gen_by_name.r567,),
                g2_A.gen_by_name.r7,gen_by_name.x0),
        (              gen_by_name.t7,(g2_A.gen_by_name.r4567,),
                g2_A.gen_by_name.r7,gen_by_name.x0),
        (              gen_by_name.t9,(g2_A.gen_by_name.r34567,),
                g2_A.gen_by_name.r7,gen_by_name.x0),
        (              gen_by_name.t11,(g2_A.gen_by_name.r234567,),
                g2_A.gen_by_name.r7,gen_by_name.x0),
        # ??

        # ??
        (              gen_by_name.x3,(g2_A.gen_by_name.r1,),
                g2_A.gen_by_name.r123456,gen_by_name.t11),
        (              gen_by_name.x3,(g2_A.gen_by_name.r12,),
                g2_A.gen_by_name.r123456,gen_by_name.t9),
        (              gen_by_name.x3,(g2_A.gen_by_name.r123,),
                g2_A.gen_by_name.r123456,gen_by_name.t7),
        (              gen_by_name.x3,(g2_A.gen_by_name.r1234,),
                g2_A.gen_by_name.r123456,gen_by_name.t5),
        (              gen_by_name.x3,(g2_A.gen_by_name.r12345,),
                g2_A.gen_by_name.r123456,gen_by_name.t3),
        (              gen_by_name.x3,(g2_A.gen_by_name.r123456,),
                g2_A.gen_by_name.r123456,gen_by_name.t1),
        # ??

        # ??    
        (              gen_by_name.x3,(g2_A.gen_by_name.r1234567,),
                g2_A.gen_by_name.r1234567,gen_by_name.x0),
        (              gen_by_name.x3,(g2_A.gen_by_name.r1,g2_A.gen_by_name.r4),
                g2_A.gen_by_name.r1,gen_by_name.t2),
        # ??



        (              gen_by_name.t12,(),
                1,gen_by_name.x2),
        (              gen_by_name.t10,(),
                g2_A.gen_by_name.r2,gen_by_name.x3),
        (              gen_by_name.t11,(g2_A.gen_by_name.r4,),
                1,gen_by_name.x1),
        (              gen_by_name.t11,(g2_A.gen_by_name.r456,),
                g2_A.gen_by_name.r56,gen_by_name.x1),

        (              gen_by_name.t11,(g2_A.gen_by_name.r45,),
                g2_A.gen_by_name.r5,gen_by_name.x0),
        (              gen_by_name.t8,(),
                g2_A.gen_by_name.r23,gen_by_name.x2),
        (              gen_by_name.t6,(),
                g2_A.gen_by_name.r234,gen_by_name.x1),
        (              gen_by_name.t2,(),
                g2_A.gen_by_name.r23456,gen_by_name.x1),
        (              gen_by_name.t4,(),
                g2_A.gen_by_name.r2345,gen_by_name.x0),



        ])

    return DA_bimodule(gen_by_name,arrows,g2_A,g2_A,name="g2_T_RHD")


g2_ID=init_g2_ID(g2_A)
g2_ID_bounded=init_g2_ID_bounded(g2_A)
g2_M_RHD=init_g2_M_RHD(g2_A)
g2_M_LHD=init_g2_M_LHD(g2_A)
g2_L_RHD=init_g2_L_RHD(g2_A)
g2_L_LHD=init_g2_L_LHD(g2_A)
g2_K_LHD=init_g2_K_LHD(g2_A)
g2_K_RHD=init_g2_K_RHD(g2_A)
g2_N_LHD=init_g2_N_LHD(g2_A)
g2_N_RHD=init_g2_N_RHD(g2_A)
g2_T_RHD=init_g2_T_RHD(g2_A)
g2_T_LHD=init_g2_T_LHD(g2_A)



