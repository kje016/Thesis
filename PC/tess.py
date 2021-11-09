# cd Desktop/Thesis/PySageMath/PC
from sage.all import *


class Node:
    def __init__(self, beliefs, state):
        self.beliefs = beliefs
        self.state = state
        self.left_child = None
        self.right_child = None
        self.values = (self.beliefs, self.state, self.right_child, self.left_child)

    def __str__(self):
        return f'beliefs:{str(self.beliefs)}, state:{self.state}, left_child: {str(self.left_child)}, right_child:{str(self.right_child)}'

    def __repr__(self):
        return repr(self.beliefs + self.state)


channel = "BSC"
r = [0, 0, 0, 0, 1, 1, 1, 1]
root = Node(r, 'l')
done, depth = False, 0
while not done:

breakpoint()
