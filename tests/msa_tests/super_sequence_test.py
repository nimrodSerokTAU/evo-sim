
from indelsim.classes.super_sequence import SuperSequence


seq = SuperSequence(10, 5)
seq.insert_item_at_position(seq[10], 11, True)

seq.insert_item_at_position(seq[3], 12, True)

seq.insert_item_at_position(seq[0], 13, True)



assert str(seq) == "0·13·1·2·3·12·4·5·6·7·8·9·10·11"