# -*- coding: utf-8 -*- 
import sys
# sys.path.append('../')

from algebraic_structures.algebra import torus_A
from algebraic_structures.basics import  Bunch_of_arrows,in_red
from algebraic_structures.da_bimodule import  (
    da_randomly_cancel_until_possible, are_equal_smart_da,
    da_check_df_is_0, composition)
from algebraic_structures.tensoring import (
    da_da_box_tensor_many_efficient_cancelations, 
    da_da_box_tensor_many_no_cancelations,
    dd_aa_box_tensor_product,
    )
from algebraic_structures.hochschild_homology import is_bounded, CH, dimHH
from algebraic_structures.visual import draw_DA_bimodule, draw_chain_complex


from input_DA_bimodules import (
    ID1,ID2,ID3,M_RHD,M_LHD,L_RHD,L_LHD,AA_id,DD_id,DD_TC
    )


DA_id_closed=dd_aa_box_tensor_product(DD_TC,AA_id)
DA_id_closed=da_randomly_cancel_until_possible(DA_id_closed)

##### computing HH=HF(Y_phi)
X=da_da_box_tensor_many_efficient_cancelations(DA_id_closed)
X=da_da_box_tensor_many_efficient_cancelations(X, X)
X=da_da_box_tensor_many_no_cancelations(ID3,X,ID3)
X.show_short()
print "\ndim(HH)=" + str(dimHH(X))
# HC=CH(X)
# HC.show()