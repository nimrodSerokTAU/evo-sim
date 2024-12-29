
from llist import sllist

class LinkedList:
    size: int
    list_container: sllist

    def __init__(self):
        self.list_container = sllist()
    
    def append(self, item) -> None:
        self.list_container.append(item)

    def insert(self, item, index) -> None:
        self.list_container.insert(index ,item)

    def get_item(self, index):
        return self.list_container[index]
    
    def __repr__(self):
        return (self.list_container.__repr__())