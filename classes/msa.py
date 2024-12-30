from super_sequence import SuperSequence
from sequence import Sequence


class Msa:
    _aligned_sequences: list[list[bool]]
    _substitutions: list[str]
    _msa_length: int
    _number_of_sequences: int
    _sequences_to_save: list[int]

    def __init__(self, super_seq: SuperSequence):
        self._number_of_sequences = super_seq.get_number_of_sequences()
        self._msa_length = super_seq.get_msa_length()
        super_seq.set_absolute_positions()
        self._aligned_sequence = {}  # Using dict instead of vector for indexed access
        
    def compute_msa(self, sequences: list[Sequence]):
        """
        Compute Multiple Sequence Alignment (MSA) based on sequences.
        
        Args:
            sequences: List of Sequence objects
        """
        for seq in sequences:
            total_size = 0
            current_position = 0
            last_position = 0
            position_difference = 0
            cumulated_difference = 1
            
            sequence_node_id = seq.get_sequence_node_id()
            
            # Initialize the aligned sequence for this ID if not exists
            if sequence_node_id not in self._aligned_sequence:
                self._aligned_sequence[sequence_node_id] = []
            
            # Handle empty sequence (only gaps)
            if len(seq) == 0:
                self._aligned_sequence[sequence_node_id].append(-self._msa_length)
                continue
            
            # Get first site and handle initial gaps
            previous_site = seq[0]  # Assuming seq supports indexing
            last_position = previous_site()
            print(last_position)
            last_position = previous_site()['absolute_position']
            if last_position > 0:
                self._aligned_sequence[sequence_node_id].append(-last_position)
                total_size += last_position
            
            # Process remaining sites
            for current_site in seq[1:]:  # Skip first element as it's already processed
                
                current_position = current_site()#['absolute_position']
                print(current_position)
                current_position = current_site()['absolute_position']

                position_difference = current_position - last_position - 1
                
                if position_difference == 0:
                    cumulated_difference += 1
                
                if position_difference > 0:
                    self._aligned_sequence[sequence_node_id].append(cumulated_difference)
                    self._aligned_sequence[sequence_node_id].append(-position_difference)
                    total_size += (cumulated_difference + position_difference)
                    cumulated_difference = 1
                
                if total_size > self._msa_length:
                    raise ValueError("Sequence lengths mismatch in fill_msa")
                
                last_position = current_position
            
            # Handle final positions
            if cumulated_difference > 0 and (total_size != self._msa_length):
                self._aligned_sequence[sequence_node_id].append(cumulated_difference)
                total_size += cumulated_difference
            
            if total_size < self._msa_length:
                self._aligned_sequence[sequence_node_id].append(-(self._msa_length - total_size))
            
            # Reset for next sequence
            cumulated_difference = 1
            last_position = 0
            total_size = 0
    
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


