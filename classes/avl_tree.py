from classes.avl_node import AVLNode
from classes.block import Block


class AVLTree:
    root: AVLNode
    nodes_count: int

    def __init__(self, bl: Block):
        self.root = AVLNode(0, bl)
        self.nodes_count = 1

    @staticmethod
    def height(node: AVLNode):
        if not node:
            return 0
        return node.height

    def balance(self, node: AVLNode):
        if not node:
            return 0
        return self.height(node.left) - self.height(node.right)

    def insert(self, root: AVLNode, bl: Block):
        if not root:
            return AVLNode(self.nodes_count, bl)
        if bl.index_in_predecessor < root.bl.index_in_predecessor:
            root.left = self.insert(root.left, bl)
        else:
            root.right = self.insert(root.right, bl)

        root.height = 1 + max(self.height(root.left), self.height(root.right))
        balance = self.balance(root)
        self.nodes_count += 1
        # Left rotation
        if balance > 1 and bl.index_in_predecessor < root.left.bl.index_in_predecessor:
            return self.right_rotate(root)

        # Right rotation
        if balance < -1 and bl.index_in_predecessor > root.right.bl.index_in_predecessor:
            return self.left_rotate(root)

        # Left-Right rotation
        if balance > 1 and bl.index_in_predecessor > root.left.bl.index_in_predecessor:
            root.left = self.left_rotate(root.left)
            return self.right_rotate(root)

        # Right-Left rotation
        if balance < -1 and bl.index_in_predecessor < root.right.bl.index_in_predecessor:
            root.right = self.right_rotate(root.right)
            return self.left_rotate(root)

        root.update_val_up()
        return root

    @staticmethod
    def update_on_same_location(avl_node: AVLNode, copy_sites_count: int | None, inserted_seq_count: int | None):
        avl_node.update_on_same_location(copy_sites_count, inserted_seq_count)

    def update_to_new_location(self, root: AVLNode, new_bl: Block, new_value: int, old_value):  # TODO: can improve
        self.delete(root, old_value) # TODO: need to update up
        self.insert(root, new_bl, new_value) # TODO: need to update up

    def delete(self, root: AVLNode, value: int):
        if not root:
            return root
        self.nodes_count -= 1
        if value < root.value:
            if root.left is not None:
                root.left = self.delete(root.left, value)
            else:
                self.delete(root, value)
        elif value > root.value:
            root.right = self.delete(root.right, value)
        else:
            if not root.left:
                temp = root.right
                root = None
                return temp
            elif not root.right:
                temp = root.left
                root = None
                return temp

            temp = self.min_value_node(root.right)
            root.value = temp.value
            root.right = self.delete(root.right, temp.value)

        if not root:
            return root

        root.height = 1 + max(self.height(root.left), self.height(root.right))
        balance = self.balance(root)

        # Left rotation
        if balance > 1 and self.balance(root.left) >= 0:
            return self.right_rotate(root)

        # Right rotation
        if balance < -1 and self.balance(root.right) <= 0:
            return self.left_rotate(root)

        # Left-Right rotation
        if balance > 1 and self.balance(root.left) < 0:
            root.left = self.left_rotate(root.left)
            return self.right_rotate(root)

        # Right-Left rotation
        if balance < -1 and self.balance(root.right) > 0:
            root.right = self.right_rotate(root.right)
            return self.left_rotate(root)

        return root

    def left_rotate(self, z):
        y = z.right
        T2 = y.left

        y.left = z
        z.father = y
        z.right = T2
        T2.father = z
        T2.update_val_up()

        z.height = 1 + max(self.height(z.left), self.height(z.right))
        y.height = 1 + max(self.height(y.left), self.height(y.right))

        return y

    def right_rotate(self, z):
        y = z.left
        T3 = y.right

        y.right = z
        z.father = y
        z.left = T3
        T3.father = z
        T3.update_val_up()

        z.height = 1 + max(self.height(z.left), self.height(z.right))
        y.height = 1 + max(self.height(y.left), self.height(y.right))

        return y

    def min_value_node(self, root: AVLNode):
        current = root
        while current.left:
            current = current.left
        return current

    def search(self, root: AVLNode, value: int, ancestor: AVLNode | None) -> tuple[AVLNode, AVLNode | None]:
        if not root or root.bl.index_in_predecessor == value or \
            (value < root.bl.index_in_predecessor and root.left is None) or \
            (value > root.bl.index_in_predecessor and root.right is None):
            return root, ancestor
        if root.bl.index_in_predecessor < value:
            return self.search(root.right, value, root)
        return self.search(root.left, value, root)

    # def insert_value(self, value: int):
    #     self.root = self.insert(self.root, value)
    #
    # def delete_value(self, value: int):
    #     self.root = self.delete(self.root, value)
    #
    # def search_value(self, value: int):
    #     return self.search(self.root, value)

    def inorder_traversal(self, root: AVLNode, res_list: list[AVLNode]):
        if root:
            if root.left is not None:
                self.inorder_traversal(root.left, res_list)
            res_list.append(root),
            if root.right is not None:
                self.inorder_traversal(root.right, res_list)