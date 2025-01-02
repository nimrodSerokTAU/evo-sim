from classes.block import Block
from classes.indel_event import IndelEvent
from enums import EventSubTypes


class SequenceNodeAsList:
    id: int
    blck_list: list[Block]
    my_length: int

    def __init__(self, seq_id: int, original_sequence_length: int):
        self.id = seq_id
        self.my_length = original_sequence_length
        self.blck_list = [Block(index_in_predecessor=0, copy_sites_count=self.my_length,
                                inserted_seq_count=0)]

    def find_block_index_and_sites_count(self, place: int) -> tuple[int, int]:
        agg_seq_length: int = 0
        for i in range(len(self.blck_list)):
            this_block: Block = self.blck_list[i]
            agg_seq_length += this_block.inserted_seq_count + this_block.copy_sites_count
            if place < agg_seq_length:
                return i, agg_seq_length
        return -1, agg_seq_length


    def find_event_sub_type(self, event: IndelEvent) -> tuple[EventSubTypes, int, int]:
        cb_index, seq_length_with_block = self.find_block_index_and_sites_count(event.place)
        # cb_index is the index of the current block
        block_at_inx: Block = self.blck_list[cb_index]
        seq_len_up_to_block: int = (seq_length_with_block - block_at_inx.inserted_seq_count -
                                    block_at_inx.copy_sites_count)
        if event.length < 0 or event.place > seq_length_with_block:
            return EventSubTypes.OUT_OF_SEQUENCE, cb_index, seq_length_with_block
        if event.is_insertion:
            if cb_index == - 1:
                return EventSubTypes.INSERTION_AT_END, cb_index, seq_length_with_block
            if event.place == 0:
                return EventSubTypes.INSERTION_AT_START, cb_index, seq_length_with_block
            if event.place < seq_len_up_to_block + block_at_inx.copy_sites_count:
                return EventSubTypes.INSERTION_INSIDE_COPIED, cb_index, seq_length_with_block
            return EventSubTypes.INSERTION_INSIDE_INSERTED, cb_index, seq_length_with_block
        else:
            if cb_index == -1 and seq_length_with_block == event.place:
                return EventSubTypes.OUT_OF_SEQUENCE, cb_index, seq_length_with_block
            position_in_block: int = event.place - seq_len_up_to_block
            if position_in_block <= block_at_inx.copy_sites_count:
                if position_in_block + event.length < block_at_inx.copy_sites_count:  # contained in copy_sites_count
                    if position_in_block > 0:
                        return EventSubTypes.DELETION_INSIDE_COPIED_CONTAINED, cb_index, seq_length_with_block
                    return EventSubTypes.DELETION_OF_COPIED, cb_index, seq_length_with_block
                return EventSubTypes.DELETION_INSIDE_COPIED_UNCONTAINED, cb_index, seq_length_with_block
            if position_in_block + event.length <= block_at_inx.inserted_seq_count:  # contained in inserted
                if position_in_block > block_at_inx.copy_sites_count:
                    return EventSubTypes.DELETION_INSIDE_INSERTED_CONTAINED, cb_index, seq_length_with_block
                return EventSubTypes.DELETION_INSIDE_INSERTED_UNCONTAINED, cb_index, seq_length_with_block
            return EventSubTypes.DELETION_OF_INSERTED, cb_index, seq_length_with_block

    def calculate_event(self, event: IndelEvent):
        event_type, cb_index, seq_length_with_block = self.find_event_sub_type(event)
        if event_type == EventSubTypes.OUT_OF_SEQUENCE:
            return
        block_at_inx: Block = self.blck_list[cb_index]
        if event.is_insertion:
            self.calculate_insertion_event(event, event_type, block_at_inx, cb_index, seq_length_with_block)
        else:
            self.calculate_deletion_event(event, event_type, block_at_inx, cb_index, seq_length_with_block)

    def calculate_insertion_event(self, event: IndelEvent, event_type: EventSubTypes, block_at_inx: Block,
                                  cb_index: int, seq_length_with_block: int):
        seq_len_up_to_block: int = (seq_length_with_block - block_at_inx.inserted_seq_count -
                                    block_at_inx.copy_sites_count)
        if event_type == EventSubTypes.INSERTION_AT_END:
            block_at_inx.inc_insert_count(event.length)
        elif event_type == EventSubTypes.INSERTION_AT_START:
            block_item = Block(index_in_predecessor=0,
                               copy_sites_count=0,
                               inserted_seq_count=event.length)
            self.blck_list.insert(0, block_item)
        elif event_type == EventSubTypes.INSERTION_INSIDE_COPIED:
            first_block_copy_count = event.place - seq_len_up_to_block
            block_item = Block(index_in_predecessor=block_at_inx.index_in_predecessor + first_block_copy_count,
                               copy_sites_count=block_at_inx.copy_sites_count - first_block_copy_count,
                               inserted_seq_count=block_at_inx.inserted_seq_count)
            block_at_inx.copy_sites_count = first_block_copy_count
            block_at_inx.update_insert_count(event.length)
            self.blck_list.insert(cb_index + 1, block_item)
        elif event_type == EventSubTypes.INSERTION_INSIDE_INSERTED:
            block_at_inx.inc_insert_count(event.length)
        self.my_length += event.length

    def calculate_deletion_event(self, event: IndelEvent, event_type: EventSubTypes, block_at_inx: Block,
                                 cb_index: int, seq_length_with_block: int):
        seq_len_up_to_block: int = (seq_length_with_block - block_at_inx.inserted_seq_count -
                                    block_at_inx.copy_sites_count)
        position_in_block: int = event.place - seq_len_up_to_block
        if event_type == EventSubTypes.DELETION_INSIDE_COPIED_CONTAINED:
            block_item = Block(
                index_in_predecessor=block_at_inx.index_in_predecessor + position_in_block + event.length,
                copy_sites_count=block_at_inx.copy_sites_count - (position_in_block + event.length),
                inserted_seq_count=block_at_inx.inserted_seq_count)
            block_at_inx.update_copy_sites_count(position_in_block)
            block_at_inx.update_insert_count(0)
            self.blck_list.insert(cb_index + 1, block_item)
            self.my_length -= event.length
        elif event_type == EventSubTypes.DELETION_OF_COPIED:
            block_item = Block(
                index_in_predecessor=block_at_inx.index_in_predecessor + position_in_block + event.length,
                copy_sites_count=block_at_inx.copy_sites_count - (position_in_block + event.length),
                inserted_seq_count=block_at_inx.inserted_seq_count)
            self.blck_list[cb_index] = block_item
            self.my_length -= event.length
        elif event_type == EventSubTypes.DELETION_INSIDE_COPIED_UNCONTAINED:
            removed_from_copied: int = block_at_inx.copy_sites_count - position_in_block
            deleted_from_insertion = min((event.length - removed_from_copied), block_at_inx.inserted_seq_count)
            block_at_inx.inc_copy_sites_count(-removed_from_copied)
            deletion_len = event.length - removed_from_copied
            self.my_length -= removed_from_copied
            self.delete_from_insertion_part(block_at_inx, deletion_len, deleted_from_insertion, seq_len_up_to_block)
        elif event_type in [EventSubTypes.DELETION_INSIDE_INSERTED_CONTAINED,
                            EventSubTypes.DELETION_INSIDE_INSERTED_UNCONTAINED, EventSubTypes.DELETION_OF_INSERTED]:  # starts inside insertion part:
            deleted_from_insertion = min(block_at_inx.inserted_seq_count -
                                         (position_in_block - block_at_inx.copy_sites_count), event.length)
            self.delete_from_insertion_part(block_at_inx, event.length, deleted_from_insertion, seq_len_up_to_block)

    def delete_from_insertion_part(self, block: Block, deletion_len: int, deleted_from_insertion: int,
                                   seq_len_up_to_block: int):
        left_to_delete_later = deletion_len - deleted_from_insertion
        block.inc_insert_count(-deleted_from_insertion)
        self.my_length -= deleted_from_insertion
        if block.is_redundant():
            self.blck_list.remove(block)
        if left_to_delete_later > 0:  # continue to next block:
            next_block_start_place = seq_len_up_to_block + block.copy_sites_count + block.inserted_seq_count
            deletion_event = IndelEvent(is_insertion=False, place=next_block_start_place, length=left_to_delete_later)
            self.calculate_event(deletion_event)

    def get_length(self) -> int:
        return self.my_length

    def get_dto(self) -> dict:
        blocks: list[str] = list(map(lambda x: x.get_dto_str(), self.blck_list))
        length: int = self.get_length()
        return {'blocks': blocks, 'length': length}

    def get_blocklist_str(self) -> dict:
        blocks: list[str] = list(map(lambda x: x.get_block_str(), self.blck_list))
        length: int = self.get_length()
        return {'blocks': blocks, 'length': length}

    def blocks_iterator(self):
        return self.blck_list