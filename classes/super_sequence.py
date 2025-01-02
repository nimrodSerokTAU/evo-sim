
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
        
        for i in range(0, root_sequence_size + 1):
            column = {'position': i, 'is_column': False}
            self._sequence.append(column)
            
        self._inserted_sequence_counter = root_sequence_size + 1

    def reference_position(self, position_ref: sllistnode):
        if (not position_ref()['is_column']):
            position_ref()['is_column'] = True
            self._msa_seq_length += 1
        

    def set_absolute_positions(self):
        i = 0
        for column in self._sequence:
            if not column.get('is_column', False): 
                continue
            column['absolute_position'] = i
            i += 1

    def insert_item_at_position(self, ref_position, position, is_to_save):
        new_column = {'position': position, 'is_column': False}

        if is_to_save:
            new_column['is_column'] = True
            self._msa_seq_length += 1
        self._sequence.insertafter(new_column, ref_position)
        return ref_position.next
    
    def insert_item_before_position(self, ref_position, position, is_to_save):
        new_column = {'position': position, 'is_column': False}

        if is_to_save:
            new_column['is_column'] = True
            self._msa_seq_length += 1
        self._sequence.insertbefore(new_column, ref_position)




    def get_num_inserted_positions(self):
        return self._inserted_sequence_counter

    def increment_num_inserted_positions(self):
        self._inserted_sequence_counter += 1
        return self._inserted_sequence_counter
    
    def increment_leaf_num(self):
        self._leaf_num += 1

    def get_original_seq_size(self):
        return self._original_seq_size
    
    def get_msa_length(self) -> int:
        return self._msa_seq_length
    
    def get_number_of_sequences(self) -> int:
        return self._num_sequences
    
    def print_seq(self):
        print(self._sequence)

    def __getitem__(self, item):
        return self._sequence.nodeat(item)
    
    def get_iterator(self):
        return self._sequence.iternodes()
    
    def __repr__(self):
        zip()
        # positions = [str(node["position"]) for node in self._sequence]
        # absolute_positions = [str(node["absolute_position"]) for node in self._sequence]
        # super_seq_str = "\n".join(map(str,zip(positions, absolute_positions)))

        super_seq_str = "·".join([str(node["position"]) for node in self._sequence])

        # super_seq_str += "\n"
        # super_seq_str += "~".join([str(node["absolute_position"]) for node in self._sequence])

        return super_seq_str

