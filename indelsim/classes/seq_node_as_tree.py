from indelsim.classes.avl_node import AVLNode
from indelsim.classes.avl_tree import AVLTree
from indelsim.classes.block import Block
from indelsim.classes.indel_event import IndelEvent
from indelsim.enums import EventSubTypes


class SequenceNodeAsTree:
    id: int
    block_tree: AVLTree
    my_length: int

    def __init__(self, seq_id: int, original_sequence_length: int):
        self.id = seq_id
        self.my_length = original_sequence_length
        self.block_tree = AVLTree(bl = Block(index_in_predecessor=0, copy_sites_count=self.my_length,
                                             inserted_seq_count=0))

    def find_event_sub_type(self, event: IndelEvent) -> tuple[EventSubTypes, AVLNode | None, int]:
        if event.length < 0 or event.place > self.block_tree.root.length_under_including or \
                (not event.is_insertion and event.place == self.block_tree.root.length_under_including):
            return EventSubTypes.OUT_OF_SEQUENCE, None, -1
        node_at_inx, position_in_block = self.block_tree.search(self.block_tree.root, event.place, event.is_insertion)
        if event.is_insertion:
            if event.place == 0:
                first_node, position_in_block_b = self.block_tree.search(self.block_tree.root, -1, event.is_insertion)
                if first_node.bl.copy_sites_count == 0:
                    return EventSubTypes.INSERTION_AT_START_ADDITION, node_at_inx, position_in_block
                return EventSubTypes.INSERTION_AT_START, node_at_inx, position_in_block
            if position_in_block < node_at_inx.bl.copy_sites_count:
                return EventSubTypes.INSERTION_INSIDE_COPIED, node_at_inx, position_in_block
            return EventSubTypes.INSERTION_INSIDE_INSERTED, node_at_inx, position_in_block
        else:
            if node_at_inx.bl is None and position_in_block == event.place:
                return EventSubTypes.OUT_OF_SEQUENCE, node_at_inx, position_in_block
            if position_in_block < node_at_inx.bl.copy_sites_count:  # <=
                if position_in_block + event.length < node_at_inx.bl.copy_sites_count:  # contained in copy_sites_count
                    if position_in_block > 0:
                        return EventSubTypes.DELETION_INSIDE_COPIED_CONTAINED_AT_MID, node_at_inx, position_in_block
                    return EventSubTypes.DELETION_INSIDE_COPIED_CONTAINED_AT_START, node_at_inx, position_in_block
                if position_in_block + event.length == node_at_inx.bl.copy_sites_count and position_in_block == 0:
                    return EventSubTypes.DELETION_OF_COPIED, node_at_inx, position_in_block
                if position_in_block == 0:
                    if event.place == 0 and node_at_inx.left is None:
                        return EventSubTypes.DELETION_ALL_COPIED_UNCONTAINED_AT_START, node_at_inx, position_in_block
                    return EventSubTypes.DELETION_ALL_COPIED_UNCONTAINED, node_at_inx, position_in_block
                return EventSubTypes.DELETION_INSIDE_COPIED_UNCONTAINED, node_at_inx, position_in_block
            if position_in_block + event.length <= node_at_inx.bl.inserted_seq_count:  # contained in inserted
                if position_in_block > node_at_inx.bl.copy_sites_count:
                    return EventSubTypes.DELETION_INSIDE_INSERTED_CONTAINED, node_at_inx, position_in_block
                return EventSubTypes.DELETION_INSIDE_INSERTED_UNCONTAINED, node_at_inx, position_in_block
            return EventSubTypes.DELETION_OF_INSERTED, node_at_inx, position_in_block

    def calculate_event(self, event: IndelEvent):
        event_type, avl_node, seq_length_with_block = self.find_event_sub_type(event)
        if event_type == EventSubTypes.OUT_OF_SEQUENCE:
            return
        if event.is_insertion:
            self.calculate_insertion_event(event, event_type, avl_node, seq_length_with_block)
        else:
            self.calculate_deletion_event(event, event_type, avl_node, seq_length_with_block)

    def calculate_insertion_event(self, event: IndelEvent, event_type: EventSubTypes, avl_node: AVLNode,
                                  position_in_block: int):
        if event_type == EventSubTypes.INSERTION_AT_END:
            avl_node.inc_on_same_location(0, event.length)
        elif event_type == EventSubTypes.INSERTION_AT_START:
            block_item = Block(index_in_predecessor=-1,
                                       copy_sites_count=0,
                                       inserted_seq_count=event.length)
            self.block_tree.insert_block(block_item)
        elif event_type == EventSubTypes.INSERTION_AT_START_ADDITION:  # this case seems covered
            self.block_tree.inc_on_same_location(avl_node, None, event.length)
        elif event_type == EventSubTypes.INSERTION_INSIDE_COPIED:  # this case seems covered
            first_block_copy_count = position_in_block
            block_item = Block(index_in_predecessor=avl_node.bl.index_in_predecessor + first_block_copy_count,
                               copy_sites_count=avl_node.bl.copy_sites_count - first_block_copy_count,
                               inserted_seq_count=avl_node.bl.inserted_seq_count)
            self.block_tree.update_on_same_location(avl_node, first_block_copy_count, event.length)
            self.block_tree.insert_block(block_item)
        elif event_type == EventSubTypes.INSERTION_INSIDE_INSERTED:
            avl_node.inc_on_same_location(0, event.length)
        self.my_length += event.length

    def calculate_deletion_event(self, event: IndelEvent, event_type: EventSubTypes, avl_node: AVLNode,
                                 position_in_block: int):
        if event_type == EventSubTypes.DELETION_INSIDE_COPIED_CONTAINED_AT_MID:
            block_item = Block(
                index_in_predecessor=avl_node.bl.index_in_predecessor + position_in_block + event.length,
                copy_sites_count=avl_node.bl.copy_sites_count - (position_in_block + event.length),
                inserted_seq_count=avl_node.bl.inserted_seq_count)
            self.block_tree.update_on_same_location(avl_node, position_in_block, 0)
            self.block_tree.insert_block(block_item)
            self.my_length -= event.length
        elif event_type == EventSubTypes.DELETION_INSIDE_COPIED_CONTAINED_AT_START:  # this case seems covered
            block_item = Block(
                index_in_predecessor=avl_node.bl.index_in_predecessor + position_in_block + event.length,
                copy_sites_count=avl_node.bl.copy_sites_count - event.length,
                inserted_seq_count=avl_node.bl.inserted_seq_count)
            self.block_tree.update_to_new_location(avl_node, block_item)
            self.my_length -= event.length
        elif event_type == EventSubTypes.DELETION_OF_COPIED:
            if avl_node.bl.index_in_predecessor > -1:
                inserted_count = avl_node.bl.inserted_seq_count
                self.block_tree.delete_node(avl_node)
                self.my_length -= inserted_count
                if inserted_count > 0:
                    self.calculate_event(IndelEvent(True, event.place, inserted_count))
            else:
                block_item = Block(
                    index_in_predecessor=avl_node.bl.index_in_predecessor + position_in_block + event.length,
                    copy_sites_count=avl_node.bl.copy_sites_count - (position_in_block + event.length),
                    inserted_seq_count=avl_node.bl.inserted_seq_count)
                self.block_tree.update_to_new_location(avl_node, block_item)
            self.my_length -= event.length
        elif event_type == EventSubTypes.DELETION_ALL_COPIED_UNCONTAINED:
            # There are cases of delete all copied with still insertions or all deleted or no insertions at all + is remove all block or not + is going forward to a new block
            deleted_from_insertion = min((event.length - avl_node.bl.copy_sites_count), avl_node.bl.inserted_seq_count)
            deleted_from_copied = avl_node.bl.copy_sites_count
            self.my_length -= (deleted_from_copied + deleted_from_insertion)
            if avl_node.bl.inserted_seq_count - deleted_from_insertion > 0:  # there are still insertions on this block
                to_insert = avl_node.bl.inserted_seq_count - deleted_from_insertion
                self.block_tree.delete_node(avl_node)
                avl_node, position_in_block = self.block_tree.search(self.block_tree.root, event.place, True)
                self.block_tree.inc_on_same_location(avl_node, 0, to_insert)
            else:  # no more insertions on this block
                self.block_tree.delete_node(avl_node)
                if event.length - deleted_from_copied - deleted_from_insertion > 0:
                    deletion_event = IndelEvent(is_insertion=False, place=event.place, length=event.length - deleted_from_copied - deleted_from_insertion)
                    self.calculate_event(deletion_event)
        elif event_type == EventSubTypes.DELETION_ALL_COPIED_UNCONTAINED_AT_START:
            deleted_from_insertion = min((event.length - avl_node.bl.copy_sites_count), avl_node.bl.inserted_seq_count)
            deleted_from_copied = avl_node.bl.copy_sites_count
            self.my_length -= (deleted_from_copied + deleted_from_insertion)
            if avl_node.bl.inserted_seq_count - deleted_from_insertion > 0:  # there are still insertions on this block
                self.block_tree.update_on_same_location(avl_node, 0, avl_node.bl.inserted_seq_count - deleted_from_insertion)
                self.block_tree.update_key_to_insert_only(avl_node)
            else:  # no more insertions on this block
                self.block_tree.delete_node(avl_node)
                if event.length - deleted_from_copied - deleted_from_insertion > 0:
                    deletion_event = IndelEvent(is_insertion=False, place=event.place, length=event.length - deleted_from_copied - deleted_from_insertion)
                    self.calculate_event(deletion_event)
        elif event_type == EventSubTypes.DELETION_INSIDE_COPIED_UNCONTAINED:
            removed_from_copied: int = avl_node.bl.copy_sites_count - position_in_block
            deleted_from_insertion = min((event.length - removed_from_copied), avl_node.bl.inserted_seq_count)
            self.block_tree.inc_on_same_location(avl_node, -removed_from_copied, None)
            deletion_len = event.length - removed_from_copied
            self.my_length -= removed_from_copied
            self.delete_from_insertion_part(avl_node, deletion_len, deleted_from_insertion, event.place)
        elif event_type in [EventSubTypes.DELETION_INSIDE_INSERTED_CONTAINED,
                            EventSubTypes.DELETION_INSIDE_INSERTED_UNCONTAINED, EventSubTypes.DELETION_OF_INSERTED]:  # starts inside insertion part:
            deleted_from_insertion = min(avl_node.bl.inserted_seq_count -
                                         (position_in_block - avl_node.bl.copy_sites_count), event.length)
            self.delete_from_insertion_part(avl_node, event.length, deleted_from_insertion, event.place)

    def delete_from_insertion_part(self, node: AVLNode, deletion_len: int, deleted_from_insertion: int,
                                   event_place: int):
        left_to_delete_later = deletion_len - deleted_from_insertion
        node.inc_on_same_location(None, -deleted_from_insertion)
        self.my_length -= deleted_from_insertion
        if node.is_redundant():
            self.block_tree.delete_node(node)
        if left_to_delete_later > 0:  # continue to next block:
            deletion_event = IndelEvent(is_insertion=False, place=event_place, length=left_to_delete_later)
            self.calculate_event(deletion_event)

    def get_length(self) -> int:
        return self.my_length

    def get_dto(self) -> dict:
        res_list: list[AVLNode] = []
        self.block_tree.inorder_traversal(self.block_tree.root, res_list)
        blocks: list[str] = list(map(lambda x: x.get_dto_str(), res_list))
        length: int = self.get_length()
        return {'blocks': blocks, 'length': length}
    
    def get_blocklist_str(self) -> dict:
        res_list: list[AVLNode] = []
        self.block_tree.inorder_traversal(self.block_tree.root, res_list)
        blocks: list[str] = list(map(lambda x: x.get_block_str(), res_list))
        length: int = self.get_length()
        return {'blocks': blocks, 'length': length}


    def get_clean_dto(self) -> dict:
        res_list: list[AVLNode] = []
        self.block_tree.inorder_traversal(self.block_tree.root, res_list)
        blocks: list[str] = list(map(lambda x: x.get_clean_dto_str(), res_list))
        length: int = self.get_length()
        return {'blocks': blocks, 'length': length}

    def blocks_iterator(self):
        res_list: list[AVLNode] = []
        self.block_tree.inorder_traversal(self.block_tree.root, res_list)
        blocks: list[Block] = list(map(lambda x: x.bl, res_list))

        return blocks