# -*- coding: utf-8 -*- 
from chain_complex import ChainComplex
from basics import AttrDict, debug
from basics import check_idempotents_match_left_left, check_idempotents_match_right_left,check_idempotents_match_right_right
from basics import Bunch_of_arrows
from random import shuffle
from algebra import Generator, dg_algebra
from aa_bimodule import AA_bimodule
from itertools import product
from algebraic_structures.left_d_module import (
    Left_D_module, 
    left_d_randomly_cancel_until_possible,
    morphism_space_for_left_D_structures,
    left_d_in_mod_gen,
    left_d_out_alg_gen,
    left_d_out_mod_gen
    )

def print_algebra_action(arrow):
    s=''
    for el in arrow[:-2]:
        s+= el.name + ' ⊗ '
    s+= arrow[-2].name
    s+= ' ⟶   '+arrow[-1].name
    print(s)

class simpler_A_inf_Algebra(object):
    def __init__(self,
                gen_by_name,
                name,
                a_inf_actions=Bunch_of_arrows([])):
        self.name=name
        self.gen_by_name=gen_by_name #dictionary of generators, where keys are their names
        self.genset=self.gen_by_name.values() #list of generators
        self.a_inf_actions=a_inf_actions
        self.check()

    def check(self):
        d_squared=Bunch_of_arrows([])
        for second_action in self.a_inf_actions:
            list_of_sets_of_first_actions=[]
            for el in second_action[:-1]:
                list_of_sets_of_first_actions.append(Bunch_of_arrows( [act 
                                                                for act 
                                                                in self.a_inf_actions 
                                                                if act[-1]==el]  )  )
            
            for (i, set_of_first_action) in enumerate(list_of_sets_of_first_actions):
                for act in set_of_first_action:\
                    d_squared[second_action[:i]+act[:-1]+second_action[i+1:] ] +=1
            
        d_squared.delete_arrows_with_even_coeff()
        if d_squared:
            debug("d_squared is not 0! the terms that are not canceled:")
            m=min([len(el) for el in d_squared.elements()])
            debug("the smallest arrow in d_squared is " + str(m))
            d_squared.show()

            raise NameError("Algebra " + self.name + " doesn't satisfy d_squared=0.")


    def show(self, restrict=[], only_diff=False):
        if restrict: gens=[gen for gen in self.genset if gen in restrict]
        else: gens=self.genset
        print('========\n'+ self.name +'.')
        if restrict: debug("RESTRICTED TO " + str(restrict) +'\n')
        print('{} generators in this algebra:'.format(len(gens)))
        for gen in sorted(gens, key=lambda gen: str(gen)) :
            print gen
        if only_diff:
            print ('A∞ actions: ...')
            print ('Differentials:')
            diffs=[act 
                for act 
                in self.a_inf_actions 
                if len(act)==2]
            for dif in diffs:
                k=1
                if restrict:
                    for a in dif:
                        if not (a in restrict): k=0
                if k: print_algebra_action(dif)
        else: 
            print ('A∞ actions:')
            for act in sorted(self.a_inf_actions, key=lambda action: (-len(action),str(action))):
                k=1
                if restrict:
                    for a in act:
                        if not (a in restrict): k=0
                if k: print_algebra_action(act)

# Fuk is simpler_A_inf_Algebra
def base_change_AA_from_fuk(Fuk,generators,left_dg_algebra,right_dg_algebra):
    def check_validity(x,action):
        if action[:-1].count(x)==1 and (action[-1] in generators):
            index_of_x=action[:-1].index(x)
            for el in action[: index_of_x]:
                if not el.name in left_dg_algebra.gen_by_name.keys(): return False
            for el in action[index_of_x+1: -1]:
                if not el.name in right_dg_algebra.gen_by_name.keys(): return False
            return True
        else: return False
    gen_by_name=AttrDict({})
    for x in generators:
        gen_by_name[''+x.name+'']=Generator(''+x.name+'')
        # gen_by_name[''+x.name+''].add_idems(x.idem.left, x.idem.right)        

        # adding idempotents
        for action in [act for act in Fuk.a_inf_actions if (len(act)==3  and (act[0].name in left_dg_algebra.idem_by_name.keys() ) and (x == act[1]) and (x == act[-1]))]:
            left_idem=left_dg_algebra.idem_by_name[action[0].name]
        for action in [act for act in Fuk.a_inf_actions if (len(act)==3 and (act[1].name in right_dg_algebra.idem_by_name.keys() ) and (x == act[0]) and (x == act[-1]))]:
            right_idem=right_dg_algebra.idem_by_name[action[1].name]
        gen_by_name[''+x.name+''].add_idems(left_idem, right_idem)

    arrows=Bunch_of_arrows([])
    for y in generators:
        for action in Fuk.a_inf_actions:
            if (len(action)==3  and (action[0].name in left_dg_algebra.idem_by_name.keys() ) and (y == action[1]) and (y == action[-1])):
                continue
            if  (len(action)==3 and (action[1].name in right_dg_algebra.idem_by_name.keys() ) and (y == action[0]) and (y == action[-1])):
                continue
            if check_validity(y,action):
                index_of_y=action[:-1].index(y)
                tuple_from_left= tuple([ (lambda z: left_dg_algebra.gen_by_name[z.name])(z) for z in action[:index_of_y]])
                tuple_from_right= tuple([ (lambda z: right_dg_algebra.gen_by_name[z.name])(z) for z in action[index_of_y+1:-1]])
                arrows[( tuple_from_left,gen_by_name[''+y.name+''],tuple_from_right,
                    gen_by_name[''+action[-1].name+''])]+=1

    return AA_bimodule(gen_by_name,
                    arrows,
                    left_dg_algebra,
                    right_dg_algebra,
            name= 'AA_from_Fuk',
            to_check=True)


# def apply(linear_map, #it is just a bunch of arrows (x_1,x_2)
#             element):
#     for arrow in linear_map:
#         if arrow==element:

def act_by_G2(B,G2,max_a_infty_action=15): # check Seidel's paper on quartic
    # form generators as a new dictionary
    B_new_gen_by_name=AttrDict(dict(B.gen_by_name))

    #form old differentials
    new_B_inf_actions=Bunch_of_arrows()
    
    # max_a_infty_action is the maximum a inf action we are computing

    new_mu=[Bunch_of_arrows([]) for i in range(max_a_infty_action)]
    new_mu[1]=Bunch_of_arrows([translate_arrow(B_new_gen_by_name,ar) for ar in B.a_inf_actions if (len(ar)==2) ]) 
    new_mu[2]=Bunch_of_arrows([translate_arrow(B_new_gen_by_name,ar) for ar in B.a_inf_actions if (len(ar)==3) ]) 
    new_mu[3]=Bunch_of_arrows([translate_arrow(B_new_gen_by_name,ar) for ar in B.a_inf_actions if (len(ar)==4) ]) 
    for gauge in G2:
        a1=B_new_gen_by_name[gauge[0].name]
        a2=B_new_gen_by_name[gauge[1].name]
        a3=B_new_gen_by_name[gauge[2].name]
        for action in [act for act in new_mu[2] if act[2]==a1]:
            new_mu[3][ action[:-1]+(a2,a3) ]+=1
        for action in [act for act in new_mu[2] if act[2]==a2]:
            new_mu[3][ (a1,)+ action[:-1]+(a3,) ]+=1
    for action in new_mu[2]:
        for gauge in G2:
            a1=B_new_gen_by_name[gauge[0].name]
            a2=B_new_gen_by_name[gauge[1].name]
            a3=B_new_gen_by_name[gauge[2].name]

            if a3==action[0]:
                new_mu[3][ (a1,a2,action[1],action[2])]+=1
            if a3==action[1]:
                new_mu[3][ (action[0],a1,a2,action[2])]+=1

    new_B_inf_actions=Bunch_of_arrows(new_B_inf_actions+new_mu[1]+new_mu[2]+new_mu[3])
    new_B_inf_actions.delete_arrows_with_even_coeff()
    return simpler_A_inf_Algebra(B_new_gen_by_name, B.name +'_gauged',new_B_inf_actions)


def algebra_from_mor_spaces_of_left_D_structures(*D_structures):
    A=D_structures[0].left_algebra

    # want name
    algebra_name='Morphism algbebra on D-structures'
    for D in D_structures:
        algebra_name=algebra_name+' '+D.name

    
    a_inf_actions=Bunch_of_arrows()
    a_inf_actions.delete_arrows_with_even_coeff()
    # want gen_by_name 
    algebra_gen_by_name=AttrDict({})
    for D1 in D_structures:
        for D2 in D_structures:
            mor_space=morphism_space_for_left_D_structures(D1,D2)
            algebra_gen_by_name.update(mor_space.gen_by_name)
            a_inf_actions=Bunch_of_arrows(a_inf_actions+mor_space.arrows)

    # want higher a_inf_actions (multiplication only, since A in the beginning was dg-algebra)
    for mor1 in algebra_gen_by_name.values():
        for mor2 in algebra_gen_by_name.values():
            f1=mor1.aux_info
            f2=mor2.aux_info
            if not left_d_out_mod_gen(f1)==left_d_in_mod_gen(f2): continue
            a_new=A.multiply(left_d_out_alg_gen(f1),left_d_out_alg_gen(f2))
            if a_new:
                a_inf_actions[(mor1,mor2,algebra_gen_by_name[str((left_d_in_mod_gen(f1),a_new,left_d_out_mod_gen(f2)))])]+=1

    a_inf_actions.delete_arrows_with_even_coeff()

    return simpler_A_inf_Algebra(algebra_gen_by_name,algebra_name,a_inf_actions)

def arguments_replace_in_action(action, old_element, new_elements):
    indices=[i for (i, e) in enumerate(action[:-1]) if e==old_element]
    if not [indices]: return Bunch_of_arrows([action])

    new_actions=[]
    for sub_elements_tuple in product(new_elements,repeat=len(indices)):
        sub_elements_list=list(sub_elements_tuple)[::]
        new_action=tuple([sub_elements_list.pop() if (e==old_element and sub_elements_list) else e for e in action])
        new_actions.append(new_action)

    return Bunch_of_arrows(new_actions)

def change_of_basis(A,a,x): # a ---> a+x, x ---> x
    A.a_inf_actions.delete_arrows_with_even_coeff()
    new_A_inf_actions=Bunch_of_arrows()
    
    # changing basis in domain
    for action in A.a_inf_actions:
        possible_actions=Bunch_of_arrows ( arguments_replace_in_action(action,  x, [x,a]) )
        new_A_inf_actions=Bunch_of_arrows( new_A_inf_actions+possible_actions)
        new_A_inf_actions.delete_arrows_with_even_coeff()
    # print new_A_inf_actions

    # changing basis in target
    arrows_to_tweak=[action for action in new_A_inf_actions if action[-1]==a]
    for arrow in arrows_to_tweak:
        new_A_inf_actions[arrow[:-1]+(x,)]+=1
    new_A_inf_actions.delete_arrows_with_even_coeff()

    return simpler_A_inf_Algebra(A.gen_by_name, A.name +'_cb',new_A_inf_actions)

def translate_arrow(new_gens_by_name,arrow):
    new_arrow=()
    for el in arrow:
        new_arrow=new_arrow+(new_gens_by_name[el.name],)
    return new_arrow

def perturb_algebra(A_old,pure_action,max_a_infty_action=25):
    x=pure_action[0]
    y=pure_action[1]

    # change of basis
    A_old.a_inf_actions.delete_arrows_with_even_coeff()
    for action in [action for action in A_old.a_inf_actions if (len(action)==2 and action[1]==y and action[0]!=x)]:
        A_old=change_of_basis(A_old, action[0], x)

    for action in [action for action in A_old.a_inf_actions if (len(action)==2 and action[0]==x and action[1]!=y)]:
        A_old=change_of_basis(A_old, y, action[1])
    A_cb=A_old

    # form generators as a new dictionary
    A_new_gen_by_name=AttrDict(dict(A_cb.gen_by_name))

    del A_new_gen_by_name[x.name]
    del A_new_gen_by_name[y.name]


    #form old differentials
    new_A_inf_actions=Bunch_of_arrows()
    
    # max_a_infty_action this is the maximum a inf action we are computing

    I=[Bunch_of_arrows([]) for i in range(max_a_infty_action)]
    I[1]=Bunch_of_arrows([ (gen, A_cb.gen_by_name[gen.name]) for gen in A_new_gen_by_name.values() ])

    new_mu=[Bunch_of_arrows([]) for i in range(max_a_infty_action)]
    new_mu[1]=Bunch_of_arrows([translate_arrow(A_new_gen_by_name,ar) for ar in A_cb.a_inf_actions if (len(ar)==2 and ar!=(x,y) ) ]) 
    new_A_inf_actions=Bunch_of_arrows(new_A_inf_actions+new_mu[1])

    all_I=I[1]

    # max_length_action=max_a_infty_action([len(action) for action in A_cb.a_inf_actions ])

    for k in range(max_a_infty_action)[2:] :
        # Тут интересно, надо композицию наборов стрелок 
        # искать в обратном порядке, тогда будет эффективно!

        # Вычисляем I[max_a_infty_action]
        # пример: I[3](a1,a2,a3)=H ( old_mu( I[2](a1,a2), I[1](a3) ) )
        old_mu_incoming_possibilities=[action[:-1] for action in A_cb.a_inf_actions if (action[-1]==y and len(action)<=k+1 )]
        for old_mu_incoming_action in old_mu_incoming_possibilities:
            list_of_lists_of_i_possibilities=[]
            for target_for_i in old_mu_incoming_action:
                list_of_lists_of_i_possibilities.append( [i_k for i_k in all_I if (i_k[-1]==target_for_i)] )

            resulting_i_actions=[p 
                                for p 
                                in product(*list_of_lists_of_i_possibilities) 
                                if sum(len(ar) for ar in p)==k+len(old_mu_incoming_action) ]
            for set_of_i_actions in resulting_i_actions:
                i_k_action=()
                for ar in set_of_i_actions:
                    i_k_action=i_k_action + ar[:-1]
                i_k_action=i_k_action + (x,)
                I[k][i_k_action]+=1

        # Вычисляем new_mu[max_a_infty_action]
        # пример: new_mu[3](a1,a2,a3)=P ( old_mu( I[2](a1,a2), I[1](a3) ) )
        old_mu_incoming_possibilities=[action for action in A_cb.a_inf_actions if (action[-1]!=y and action[-1]!=x and len(action)<=k+1 )]
        for old_mu_incoming_action in old_mu_incoming_possibilities:
            list_of_lists_of_i_possibilities=[]
            for target_for_i in old_mu_incoming_action[:-1]:
                list_of_lists_of_i_possibilities.append( [i_k for i_k in all_I if (i_k[-1]==target_for_i)] )

            resulting_i_actions=[p 
                                for p 
                                in product(*list_of_lists_of_i_possibilities) 
                                if sum(len(ar) for ar in p)==k+len(old_mu_incoming_action)-1 ]

            for set_of_i_actions in resulting_i_actions:
                new_mu_k_action=()
                for ar in set_of_i_actions:
                    new_mu_k_action=new_mu_k_action + ar[:-1]
                new_mu_k_action=new_mu_k_action + (A_new_gen_by_name[old_mu_incoming_action[-1].name],)
                new_mu[k][new_mu_k_action]+=1


        all_I=Bunch_of_arrows(all_I+I[k])
        new_A_inf_actions=Bunch_of_arrows(new_A_inf_actions+new_mu[k])
    all_I.delete_arrows_with_even_coeff()
    new_A_inf_actions.delete_arrows_with_even_coeff()
    debug("maximum length of I is "+str(max([len(a)-1 for a in all_I])))
    debug("maximum length of mu is "+str(max([len(a)-1 for a in new_A_inf_actions])))
    return simpler_A_inf_Algebra(A_new_gen_by_name, A_cb.name +'_pert',new_A_inf_actions)

def rename_generators(A,list_of_tuples_to_rename):
    new_generators=AttrDict({})
    # for t in list_of_tuples_to_rename:
    #     new_generators[t[1]]=Generator(name=t[1], aux_info=A.gen_by_name[t[0]])
    
    for gen in A.gen_by_name.values():
        if gen.name in [t[0] for t in list_of_tuples_to_rename]:
            t=next(t for t in list_of_tuples_to_rename if t[0]==gen.name)
            new_generators[t[1]]=Generator(name=t[1], aux_info=A.gen_by_name[t[0]])
        else: new_generators[gen.name]=Generator(name=gen.name, aux_info=gen)

    new_A_inf_actions=Bunch_of_arrows({})
    for action in A.a_inf_actions:
        new_action=()
        for el in action:
            new_action=new_action+( next(gen for gen in new_generators.values() if gen.aux_info==el),)
        new_A_inf_actions[new_action]+=1

    return simpler_A_inf_Algebra(new_generators, A.name +'_renamed',new_A_inf_actions)


def u_i_rename_generators(A):
    new_generators=AttrDict({})
    i=0
    for gen in A.gen_by_name.values():
        new_generators["u_"+str(i)]=Generator(name="u_"+str(i), aux_info=gen)
        i+=1


    new_A_inf_actions=Bunch_of_arrows({})
    for action in A.a_inf_actions:
        new_action=()
        for el in action:
            new_action=new_action+( next(gen for gen in new_generators.values() if gen.aux_info==el),)
        new_A_inf_actions[new_action]+=1

    return simpler_A_inf_Algebra(new_generators, A.name +'_renamed',new_A_inf_actions)

def opposite_algebra(A_inf):
    gen_by_name=AttrDict({})
    for gen in A_inf.genset:
        gen_by_name[gen.name+'*']=Generator(gen.name+'*')

    new_A_inf_actions=Bunch_of_arrows([])

    for action in A_inf.a_inf_actions:
        new_action=()
        for element in action:
            new_action=new_action+(gen_by_name[element.name+'*'],)
        new_A_inf_actions[ new_action[:-1][::-1]+(new_action[-1],)]+=1

    return simpler_A_inf_Algebra(gen_by_name,'opposite of '+A_inf.name,new_A_inf_actions)