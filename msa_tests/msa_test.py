import sys
import pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd() / "classes"))


from super_sequence import SuperSequence
from block import Block
from sequence import Sequence
from msa import Msa

super_seq = SuperSequence(50,3)

root_seq = Sequence(super_seq, True, 0)
root_seq.init_root_seq()

child_seq = Sequence(super_seq, True, 1)
blocks = [Block(0,10,5),Block(15,35,4)]
child_seq.generate_sequence(blocks, root_seq)

grandchild_seq = Sequence(super_seq, True, 2)
blocks = [Block(0,5,2),Block(5,15,0),Block(30,24,1)]
grandchild_seq.generate_sequence(blocks, child_seq)

msa = Msa(super_seq)

msa.compute_msa([root_seq, child_seq, grandchild_seq])

print(msa.generate_msa_string_without_subs())