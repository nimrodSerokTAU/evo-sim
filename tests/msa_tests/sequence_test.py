
from indelsim.classes.super_sequence import SuperSequence
from indelsim.classes.block import Block
from indelsim.classes.sequence import Sequence


super_seq = SuperSequence(50,3)

root_seq = Sequence(super_seq, True, 0)
root_seq.init_root_seq()
print(root_seq)

child_seq = Sequence(super_seq, False, 1)
blocks = [Block(0,10,5),Block(15,35,4)]
child_seq.generate_sequence(blocks, root_seq)

print(child_seq)
grandchild_seq = Sequence(super_seq, True, 1)
blocks = [Block(0,5,2),Block(5,15,0),Block(30,24,1)]
grandchild_seq.generate_sequence(blocks, child_seq)
print(grandchild_seq)
print(super_seq)