
import codecs
from string import strip
from random import shuffle
import sys
from time import sleep

# Contains a sentence's parse tree
class Tree:
  def TreebankRecursion(self, children, node_source, current_id, current_node):
    for ii in children.get(current_id, []):
      id, word, role, pos = node_source[ii]
      #print word
      child_node = Node(word, role, pos)
      current_node.children.add(child_node)

      self.TreebankRecursion(children,
                             node_source,
                             ii,
                             child_node)
    return current_node

  def __init__(self, head, children, nodes):
        self.head_id = head
        # print head, nodes, children
        assert self.head_id in nodes

        id, word, role, pos = nodes[self.head_id]
        self.head = self.TreebankRecursion(children, nodes, self.head_id,
                                           Node(word, role, pos))


  def PrintString(self, use_relation, min_nodes = 2):
        trees = []
        # Minipar had sentence fragments
        tree_string = self.head.PrintString(use_relation)
        # print "Returned string", tree_string
        if tree_string.count("-1") > min_nodes:
            return tree_string
        else:
            return ""

class Node:
    def __init__(self, word, relation, pos):
	self.word = word
	self.relation = relation
        self.pos = pos
	self.children = set([])

    def __str__(self):
        return "{%s, %s, %s}" % (self.word, self.relation, self.pos)

    def PrintString(self, use_relation):
        word = str(self.word)
        s = ""
        s += word
        s += " "
        if use_relation:
            s += str(self.relation)
        else:
            s += str(self.pos)
        s += " "
        s += str(len(self.children))
        s += " "

        for i in self.children:
            s += i.PrintString(use_relation)

        s += "-1"
        s += " "
        return s
