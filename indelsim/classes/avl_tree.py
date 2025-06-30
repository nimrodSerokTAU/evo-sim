from indelsim.classes.avl_node import AVLNode
from indelsim.classes.block import Block
# https://visualgo.net/en/bst


class AVLTree:
    root: AVLNode
    running_nodes_count: int

    def __init__(self, bl: Block):
        self.root = AVLNode(0, bl)
        self.running_nodes_count = 1

    def insert_block(self, bl: Block):
        nodes_to_update_val: set[AVLNode] = set()
        self.root = self.insert(self.root, None, bl, nodes_to_update_val)
        self.running_nodes_count += 1
        self.update_all_needed_nodes(nodes_to_update_val)

    @staticmethod
    def update_on_same_location(avl_node: AVLNode, copy_sites_count: int | None, inserted_seq_count: int | None):
        avl_node.update_on_same_location(copy_sites_count, inserted_seq_count)

    @staticmethod
    def inc_on_same_location(avl_node: AVLNode, delta_copied_count: int | None, delta_inserted_count: int | None):
        avl_node.inc_on_same_location(delta_copied_count, delta_inserted_count)

    @staticmethod
    def update_key_to_insert_only(avl_node: AVLNode):
        avl_node.update_key_to_insert_only()

    def update_to_new_location(self, node_to_update: AVLNode, new_bl: Block):  # TODO: can improve
        nodes_to_update_val: set[AVLNode] = set()
        self.root = self.delete(self.root, node_to_update.bl.index_in_predecessor, nodes_to_update_val)
        self.update_all_needed_nodes(nodes_to_update_val)
        self.insert_block(new_bl)


    def delete_node(self, node_to_delete: AVLNode):
        nodes_to_update_val: set[AVLNode] = set()
        self.root = self.delete(self.root, node_to_delete.bl.index_in_predecessor, nodes_to_update_val)
        self.update_all_needed_nodes(nodes_to_update_val)

    @staticmethod
    def height(node: AVLNode):
        if not node:
            return 0
        return node.height

    def balance(self, node: AVLNode):
        if not node:
            return 0
        return self.height(node.left) - self.height(node.right)

    def insert(self, current_node: AVLNode, father: AVLNode | None, bl: Block, nodes_to_update_val: set[AVLNode]) -> AVLNode:
        if not current_node:
            new_node = AVLNode(self.running_nodes_count, bl)
            new_node.set_a_father(father)
            nodes_to_update_val.add(current_node)
            return new_node
        if bl.index_in_predecessor < current_node.bl.index_in_predecessor:
            current_node.left = self.insert(current_node.left, current_node, bl, nodes_to_update_val)
        else:
            current_node.right = self.insert(current_node.right, current_node, bl, nodes_to_update_val)

        current_node.height = 1 + max(self.height(current_node.left), self.height(current_node.right))
        nodes_to_update_val.add(current_node)
        balance = self.balance(current_node)
        # Left rotation
        if balance > 1 and bl.index_in_predecessor < current_node.left.bl.index_in_predecessor:
            return self.right_rotate(current_node, nodes_to_update_val)

        # Right rotation
        if balance < -1 and bl.index_in_predecessor > current_node.right.bl.index_in_predecessor:
            return self.left_rotate(current_node, nodes_to_update_val)

        # Left-Right rotation
        if balance > 1 and bl.index_in_predecessor > current_node.left.bl.index_in_predecessor:
            current_node.left = self.left_rotate(current_node.left, nodes_to_update_val)
            return self.right_rotate(current_node, nodes_to_update_val)

        # Right-Left rotation
        if balance < -1 and bl.index_in_predecessor < current_node.right.bl.index_in_predecessor:
            current_node.right = self.right_rotate(current_node.right, nodes_to_update_val)
            return self.left_rotate(current_node, nodes_to_update_val)

        # current_node.update_length_under_including()
        return current_node

    def delete(self, current_node: AVLNode, value: int, nodes_to_update_val: set[AVLNode]) -> AVLNode:
        nodes_to_update_val.add(current_node)
        if not current_node:
            stop = True
            return current_node
        if value < current_node.bl.index_in_predecessor:
            if current_node.left is not None:
                current_node.left = self.delete(current_node.left, value, nodes_to_update_val)
            else:
                self.delete(current_node, value, nodes_to_update_val)
        elif value > current_node.bl.index_in_predecessor:
            current_node.right = self.delete(current_node.right, value, nodes_to_update_val)
        elif current_node.left is not None or current_node.right is not None:
            if not current_node.left:
                temp = current_node.right
                temp.set_a_father(current_node.father)
                current_node = None
                nodes_to_update_val.add(temp)
                return temp
            elif not current_node.right:
                temp = current_node.left
                temp.set_a_father(current_node.father)
                current_node = None
                nodes_to_update_val.add(temp)
                return temp

            temp = self.min_value_node(current_node.right)
            current_node.bl = temp.bl
            current_node.right = self.delete(current_node.right, temp.bl.index_in_predecessor, nodes_to_update_val)

        else:
            current_node = None

        if not current_node:
            return current_node

        current_node.height = 1 + max(self.height(current_node.left), self.height(current_node.right))
        balance = self.balance(current_node)

        # Left rotation
        if balance > 1 and self.balance(current_node.left) >= 0:
            return self.right_rotate(current_node, nodes_to_update_val)

        # Right rotation
        if balance < -1 and self.balance(current_node.right) <= 0:
            return self.left_rotate(current_node, nodes_to_update_val)

        # Left-Right rotation
        if balance > 1 and self.balance(current_node.left) < 0:
            current_node.left = self.left_rotate(current_node.left, nodes_to_update_val)
            return self.right_rotate(current_node, nodes_to_update_val)

        # Right-Left rotation
        if balance < -1 and self.balance(current_node.right) > 0:
            current_node.right = self.right_rotate(current_node.right, nodes_to_update_val)
            return self.left_rotate(current_node, nodes_to_update_val)

        return current_node

    def left_rotate(self, grandfather: AVLNode, nodes_to_update_val: set[AVLNode]):
        right_child = grandfather.right
        left_grandchild = right_child.left

        right_child.left = grandfather
        if grandfather is not None:
            right_child.set_a_father(grandfather.father)
            grandfather.set_a_father(right_child)
        grandfather.right = left_grandchild
        if left_grandchild is not None:
            left_grandchild.set_a_father(grandfather)
            left_grandchild.update_length_under_including()
        else:
            grandfather.update_length_under_including()

        grandfather.height = 1 + max(self.height(grandfather.left), self.height(grandfather.right))
        right_child.height = 1 + max(self.height(right_child.left), self.height(right_child.right))
        nodes_to_update_val.update([right_child, left_grandchild, grandfather])
        return right_child

    def right_rotate(self, grandfather: AVLNode, nodes_to_update_val: set[AVLNode]):
        left_child = grandfather.left
        right_grandchild = left_child.right

        left_child.right = grandfather
        if grandfather is not None:
            left_child.set_a_father(grandfather.father)
            grandfather.set_a_father(left_child)
        grandfather.left = right_grandchild
        if right_grandchild is not None:
            right_grandchild.set_a_father(grandfather)
            right_grandchild.update_length_under_including()
        else:
            grandfather.update_length_under_including()

        grandfather.height = 1 + max(self.height(grandfather.left), self.height(grandfather.right))
        left_child.height = 1 + max(self.height(left_child.left), self.height(left_child.right))
        nodes_to_update_val.update([left_child, right_grandchild, grandfather])
        return left_child

    def min_value_node(self, root: AVLNode):
        current = root
        while current.left:
            current = current.left
        return current

    def search(self, root: AVLNode, position_in_block: int, is_insertion: bool) -> tuple[AVLNode, int]:
        if not root:
            return root, 0
        if root.left is not None:
            if ((position_in_block <= root.left.length_under_including and is_insertion) or
                    (position_in_block < root.left.length_under_including and not is_insertion)):
                return self.search(root.left, position_in_block, is_insertion)
            position_in_block -= root.left.length_under_including
        if ((position_in_block < root.get_my_own_length() and not is_insertion) or
                (position_in_block <= root.get_my_own_length() and is_insertion)):
            return root, position_in_block
        elif root.right is not None:
            return self.search(root.right, position_in_block - root.get_my_own_length(), is_insertion)
        else:
            return root, position_in_block

    def inorder_traversal(self, root: AVLNode, res_list: list[AVLNode]):
        if root:
            if root.left is not None:
                self.inorder_traversal(root.left, res_list)
            res_list.append(root),
            if root.right is not None:
                self.inorder_traversal(root.right, res_list)

    @staticmethod
    def update_all_needed_nodes(nodes_to_update_val: set[AVLNode]):
        nodes_list = sorted([n for n in nodes_to_update_val if n is not None], key=lambda x: x.height)
        for n in nodes_list:
            n.update_length_under_including()


    def debug_tree_structure(self):
        if self.root is None:
            return True
        return self.root.debug_node_structure()