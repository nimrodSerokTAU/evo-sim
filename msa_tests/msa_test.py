import sys
import pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd() / "classes"))


from super_sequence import SuperSequence
from block import Block
from sequence import Sequence
from msa import Msa

super_seq = SuperSequence(20,3)

root_seq = Sequence(super_seq, True, 0)
root_seq.init_root_seq()

child_seq = Sequence(super_seq, True, 1)
blocks = [Block(-1,0,1),Block(0,11,5),Block(11,10,4)]
child_seq.generate_sequence(blocks, root_seq)

grandchild_seq = Sequence(super_seq, True, 2)
blocks = [Block(0,5,2),Block(5,15,0),Block(23,8,1)]
grandchild_seq.generate_sequence(blocks, child_seq)


msa = Msa(super_seq)
msa.compute_msa([root_seq, child_seq, grandchild_seq])
print(msa)