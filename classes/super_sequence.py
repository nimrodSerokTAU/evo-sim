
from llist import sllist,sllistnode #https://ajakubek.github.io/python-llist/index.html#llist.sllist


class SuperSequence:
    _sequence: sllist
    _msa_seq_length: int
    _leaf_num: int
    _original_seq_size: int
    _inserted_sequence_counter: int

    def __init__(self, root_sequence_size, num_sequences):
        self._original_seq_size = root_sequence_size
        self._msa_seq_length = 0
        self._leaf_num = 0
        self._num_sequences = num_sequences
        self._sequence = sllist()
        
        for i in range(1, root_sequence_size + 1):
            column = {'position': i, 'is_column': False}
            self._sequence.append(column)
            
        self._inserted_sequence_counter = root_sequence_size + 1

    def reference_position(self, position_ref: sllistnode):
        if (not position_ref().isColumn):
            position_ref().isColumn = True
            self._msa_seq_length += 1
        

    def set_absolute_positions(self):
        i = 0
        for column in self._sequence:
            if not column.get('is_column', False): 
                continue
            column['absolute_position'] = i
            i += 1

    def insert_item_at_position(self, ref_position, item, is_to_save):
        new_column = {'item': item, 'is_column': False}

        if is_to_save:
            new_column['is_column'] = True
            self._msa_seq_length += 1

        return self._sequence.insertafter(new_column, ref_position)



    def get_num_inserted_positions(self):
        return self._inserted_sequence_counter

    def increment_num_inserted_positions(self):
        self._inserted_sequence_counter += 1

    def get_original_seq_size(self):
        return self._original_seq_size
    
    def get_msa_length(self):
        return self._msa_seq_length
    
    def print_seq(self):
        print(self._sequence)

    def __getitem__(self, item):
        return self._sequence.nodeat(item)
    
    def get_iterator(self):
        return self._sequence.iternodes()


# seq = SuperSequence(10, 5)
# seq.insert_item_at_position(seq[3], 30, True)
# seq.print_seq()