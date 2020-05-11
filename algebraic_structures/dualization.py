import sys
from .basics import AttrDict, Generator, Bunch_of_arrows, debug
sys.path.append('../')
from .right_a_module import Right_A_module
from .right_d_module import Right_D_module
from .left_d_module import Left_D_module
from .left_a_module import Left_A_module
from .dd_bimodule import DD_bimodule
from .aa_bimodule import AA_bimodule

def dualization_right_to_left_A(Right_A):
    gen_by_name=AttrDict({})
    for gen in Right_A.genset:
        gen_by_name[gen.name+'*']=Generator(gen.name+'*')
        gen_by_name[gen.name+'*'].add_idems(gen.idem.right,0)

    left_a_arrows=Bunch_of_arrows([])
    for arrow in Right_A.right_a_arrows:
        left_a_arrows[(arrow[1],gen_by_name[arrow[2].name+'*'],gen_by_name[arrow[0].name+'*'])]+=1

    return Left_A_module(gen_by_name,left_a_arrows,left_algebra=Right_A.right_algebra,name=Right_A.name + '_dual')

def dualization_left_to_right_D(Left_D):
    gen_by_name=AttrDict({})
    for gen in Left_D.genset:
        gen_by_name[gen.name+'*']=Generator(gen.name+'*')
        gen_by_name[gen.name+'*'].add_idems(0,gen.idem.left)

    right_d_arrows=Bunch_of_arrows([])
    for arrow in Left_D.left_d_arrows:
        right_d_arrows[(gen_by_name[arrow[2].name+'*'],gen_by_name[arrow[0].name+'*'],arrow[1])]+=1

    return Right_D_module(gen_by_name,right_d_arrows,right_algebra=Left_D.left_algebra,name=Left_D.name + '_dual')


def dualization_of_DD(DD):
    gen_by_name=AttrDict({})
    for gen in DD.genset:
        gen_by_name[gen.name+'*']=Generator(gen.name+'*')
        gen_by_name[gen.name+'*'].add_idems(gen.idem.right,gen.idem.left)

    new_arrows=Bunch_of_arrows([])
    for old_arrow in DD.dd_arrows:
        new_arrows[(gen_by_name[old_arrow[2].name+'*'],
            old_arrow[3],gen_by_name[old_arrow[0].name+'*'],old_arrow[1])]+=1

    return DD_bimodule(gen_by_name,new_arrows,left_algebra=DD.right_algebra,right_algebra=DD.left_algebra,name=DD.name + '_dual')

def dualization_of_AA(AA,to_check=True):
    gen_by_name=AttrDict({})
    for gen in AA.genset:
        gen_by_name[gen.name+'*']=Generator(gen.name+'*')
        gen_by_name[gen.name+'*'].add_idems(gen.idem.right,gen.idem.left)

    new_arrows=Bunch_of_arrows([])
    for old_arrow in AA.aa_arrows:
        new_arrows[(old_arrow[2],gen_by_name[old_arrow[3].name+'*'],old_arrow[0]
            ,gen_by_name[old_arrow[1].name+'*'],)]+=1

    return AA_bimodule(gen_by_name,new_arrows,left_algebra=AA.right_algebra,right_algebra=AA.left_algebra,name=AA.name + '_dual',to_check=to_check)




