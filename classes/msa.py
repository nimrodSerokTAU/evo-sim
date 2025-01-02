from super_sequence import SuperSequence
from sequence import Sequence


class Msa:
    _aligned_sequences: dict[int, str]
    _substitutions: list[str]
    _msa_length: int
    _number_of_sequences: int
    _sequences_to_save: list[int]

    def __init__(self, super_seq: SuperSequence):
        self._number_of_sequences = super_seq.get_number_of_sequences()
        self._msa_length = super_seq.get_msa_length()
        super_seq.set_absolute_positions()
        self._aligned_sequences = {}
        
    def compute_msa(self, sequences: list[Sequence]):
        """
        Compute Multiple Sequence Alignment (MSA) based on sequences.
        
        Args:
            sequences: List of Sequence objects
        """
        for seq in sequences:            
            sequence_node_id = seq.get_sequence_node_id()
            
            # Initialize the aligned sequence for this ID if not exists
            if sequence_node_id not in self._aligned_sequences:
                self._aligned_sequences[sequence_node_id] = []
            
            # Handle empty sequence (only gaps)
            if len(seq) == 0:
                self._aligned_sequences[sequence_node_id].append(-self._msa_length)
                continue
            
            previous_absolute_position = -1
            seq_str = ""
            for site in seq[1:]:
                current_absolute_position = site()['absolute_position']
                position_difference = current_absolute_position-previous_absolute_position
                    
                if (position_difference) > 1:
                    seq_str += "-" * ((position_difference) - 1)
                    seq_str += "X"
                else:
                    seq_str += "X"

                previous_absolute_position = current_absolute_position
            # print()
            # print(self._msa_length, len(seq_str))
            if (self._msa_length - len(seq_str))  >  0:
                seq_str += "-" * (self._msa_length - len(seq_str))


            self._aligned_sequences[sequence_node_id] = seq_str
    
    def get_aligned_sequence(self):
        """Return the aligned sequence dictionary"""
        return self._aligned_sequence
    
    def generate_msa_string_without_subs(self):
        """
        Generate MSA string representation using 'A' for matches and '-' for gaps.
        Returns a string with each sequence on a new line.
        """
        msa_lines = []
        for id in [0,1,2]:#self._sequences_to_save:
            line = ''
            for str_size in self._aligned_sequence[id]:
                # Convert negative numbers to positive for string multiplication
                size = abs(str_size)
                # Use '-' for negative numbers (gaps) and 'A' for positive (matches)
                char = '-' if str_size < 0 else 'A'
                line += char * size
            msa_lines.append(line)
        
        return '\n'.join(msa_lines)

    def __repr__(self):
        msa_str = ""
        for key,val in self._aligned_sequences.items():
            msa_str += f">{key}\n{val}\n"
        return msa_str
