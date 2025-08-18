from pathlib import Path
from indelsim.classes.super_sequence import SuperSequence
from indelsim.classes.sequence import Sequence


class Msa:
    _aligned_sequences: dict[int, str]
    _id_to_name: dict[int,str]
    _substitutions: list[str] # dict[int, list[int]]
    _msa_length: int
    _number_of_sequences: int
    _sequences_to_save: list[Sequence]
    _is_from_naive: bool


    def __init__(self, super_seq: SuperSequence=None):
        if super_seq is not None:
            self._number_of_sequences = super_seq.get_number_of_sequences()
            self._msa_length = super_seq.get_msa_length()
            super_seq.set_absolute_positions()
            self._is_from_naive = False
        else:
            self._is_from_naive = True
        self._aligned_sequences = {}
        self._id_to_name = {}
        
    def compute_msa(self):
        """
        Compute Multiple Sequence Alignment (MSA) based on sequences.
        
        Args:
            sequences: List of Sequence objects
        """
        if self._is_from_naive:
            for idx, seq in self._aligned_sequences.items():
                seq_str = "".join(["X" if (site != -1) else "-" for site in seq])
                self._aligned_sequences[idx] = seq_str
            return

        for idx,seq in enumerate(self._sequences_to_save):            
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
            if (self._msa_length - len(seq_str))  >  0:
                seq_str += "-" * (self._msa_length - len(seq_str))


            self._aligned_sequences[sequence_node_id] = seq_str
            self._sequences_to_save[idx] = 0

    def compute_msa_to_disk(self, output_path: Path) -> Path:
        """
        Compute MSA and write directly to disk to save memory.
        
        Args:
            output_path: Path where the MSA file should be saved
        """
        # Ensure output directory exists
        # output_path.parent.mkdir(parents=True, exist_ok=True)
        if self._is_from_naive:
            with open(output_path, 'a') as f:

                for idx, seq in self._aligned_sequences.items():
                    seq_str = "".join(["X" if (site != -1) else "-" for site in seq])
                    fasta_line_str = f">{self._id_to_name[idx]}\n{seq_str}\n"
                    f.write(fasta_line_str)
            self._msa_length = len(seq_str)
            self._number_of_sequences = len(self._aligned_sequences)
            self._aligned_sequences.clear()
            return
        
        with open(output_path, 'a') as f:
            for idx,seq in enumerate(self._sequences_to_save):            
                sequence_node_id = seq.get_sequence_node_id()
                
                # Write FASTA header
                seq_name = self._id_to_name.get(sequence_node_id, str(sequence_node_id))
                f.write(f">{seq_name}\n")
                
                # Handle empty sequence (only gaps)
                if len(seq) == 0:
                    f.write("-" * self._msa_length + "\n")
                    continue
                
                # Build sequence string (same logic as original compute_msa)
                previous_absolute_position = -1
                seq_str = ""
                
                for site in seq[1:]:
                    current_absolute_position = site()['absolute_position']
                    position_difference = current_absolute_position - previous_absolute_position
                    
                    if position_difference > 1:
                        seq_str += "-" * (position_difference - 1) + "X"
                    else:
                        seq_str += "X"
                        
                    previous_absolute_position = current_absolute_position
                
                # Add trailing gaps if needed
                if (self._msa_length - len(seq_str)) > 0:
                    seq_str += "-" * (self._msa_length - len(seq_str))
                
                # Write sequence directly to file instead of storing in _aligned_sequences
                f.write(seq_str + "\n")
                # seq_str goes out of scope and gets garbage collected
                self._sequences_to_save[idx] = 0

            self._sequences_to_save.clear()
            self._aligned_sequences.clear()
        return output_path

    def get_aligned_sequence(self):
        """Return the aligned sequence dictionary"""
        return self._aligned_sequence
    
    
    def msa_str_rep(self):
        msa_str = ""
        for key,val in self._aligned_sequences.items():
            msa_str += f">{self.id_to_name[key]}\n{val}\n"
        return msa_str
    

    def __repr__(self):
        return self.msa_str_rep()
    
