# -*- coding: utf-8 -*- 
from collections import Counter

def simpler_names_for_generators(module):
    for ind,gen in enumerate(sorted(module.gen_by_name.values(), key=str)):
        module.gen_by_name.pop(gen.name, None)
        gen.name='x'+str(ind)
        module.gen_by_name['x'+str(ind)]=gen

def rename_gen_in_module(module,name1,name2):
    gen=module.gen_by_name[name1]
    module.gen_by_name.pop(name1, None)
    gen.name=name2
    module.gen_by_name[name2]=gen

# turns printing in red:
def debug(whatever):
    print('\033[91m' + str(whatever) + '\033[0m')

# prints in red:
def in_red(whatever):
    print('\033[91m' + str(whatever) + '\033[0m')

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class Bunch_of_arrows(Counter):
    def show(self):
        for arrow_as_a_tuple in sorted(self, key=str):
            print(str(arrow_as_a_tuple) + '     *' + str(self[arrow_as_a_tuple]))

    def delete_arrows_with_even_coeff(self): 
        arrows_to_delete=[]
        for arrow in self:
            if self[arrow] % 2 ==0:
                arrows_to_delete.append(arrow)
            else:
                self[arrow]=1

        for ar in arrows_to_delete:
            del self[ar]

class Generator(object):
    def __init__(self, name,tex_name="WRONG:TEX NAME NOT SPECIFIED",reverse_tex_name="WRONG:REVERSE TEX NAME NOT SPECIFIED",aux_info=None):
        self.name = name
        self.tex_name=tex_name
        self.factorizations=[]
        self.tex_name=tex_name
        self.reverse_tex_name=reverse_tex_name
        self.aux_info=aux_info

    def add_idems(self,idem1,idem2):
        self.idem=AttrDict({"left":idem1, "right":idem2})

    def add_factorizations(self,*factorizations): #for algebra generators only!
        self.factorizations=getattr(self,'factorizations', [])
        for factorization in factorizations:
            self.factorizations.append(factorization)

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

def check_idempotents_match_right_left(l,m):
    if l==1: return True
    if m==1: return True
    if not (l.idem.right==m.idem.left): 
        print("These elements right_left idempotents don't match! " + str((l,m)))
    return (l.idem.right==m.idem.left)


def check_idempotents_match_right_right(l,m):
    if l==1: return True
    if m==1: return True
    if not (l.idem.right==m.idem.right): 
        print("These elements right_right idempotents don't match! " + str((l,m)))
    return (l.idem.right==m.idem.right)


def check_idempotents_match_left_left(l,m):
    if l==1: return True
    if m==1: return True
    if not (l.idem.left==m.idem.left): 
        print("These elements left_left idempotents don't match! " + str((l,m)))
    return (l.idem.left==m.idem.left)



def is_cyclic(g):
    """Return True if the directed graph g has a cycle.
    g must be represented as a dictionary mapping vertices to
    iterables of neighbouring vertices. For example:

    >>> cyclic({1: (2,), 2: (3,), 3: (1,)})
    True
    >>> cyclic({1: (2,), 2: (3,), 3: (4,)})
    False
    """
    path = set()

    def visit(vertex):
        path.add(vertex)
        for neighbour in g.get(vertex, ()):
            if neighbour in path or visit(neighbour):
                return True
        path.remove(vertex)

        return False

    return any(visit(v) for v in g)

