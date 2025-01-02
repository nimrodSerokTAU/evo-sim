import sys
import pathlib
print(pathlib.Path.cwd())
print(sys.path)
sys.path.append(str(pathlib.Path.cwd() / "classes"))

from super_sequence import SuperSequence


seq = SuperSequence(10, 5)
seq.insert_item_at_position(seq[9], 11, True)

seq.insert_item_at_position(seq[3], 12, True)

seq.insert_item_at_position(seq[0], 13, True)



seq.print_seq()
print(str(seq))

assert str(seq) == "13·1·2·3·4·12·5·6·7·8·9·10·11"