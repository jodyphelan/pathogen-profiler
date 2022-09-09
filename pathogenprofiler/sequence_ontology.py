import json

class sequence_ontology:
    def __init__(self,json_file):
        self.data = json.load(open(json_file))
        self._id2label = {}
        self._labels2id = {}
        for t in self.data['graphs'][0]['nodes']:
            if "lbl" in t:
                self._id2label[t['id'].split("/")[-1]] = t['lbl']
                self._labels2id[t['lbl']] = t['id'].split("/")[-1]
        self.edges = []
        for e in self.data['graphs'][0]['edges']:
            self.edges.append((e['obj'].split("/")[-1],e['sub'].split("/")[-1]))
    def get_children(self,node):
        children = set()
        for e in self.edges:
            if e[0]==node:
                children.add(e[1])
        return children
    def get_sucessors_id(self,current_node):
        nodes = set([current_node])
        def recurse(n):
            children = self.get_children(n)
            if len(children)==0:
                return 
            else:
                for c in children:
                    nodes.add(c)
                    recurse(c)
            return 
        recurse(current_node)
        return nodes
    def get_sucessors(self,label):
        if self.is_valid_term(label):
            return [self.id2label(x) for x in  self.get_sucessors_id(self.label2id(label))]
        else:
            return set()
    def label2id(self,label):
        return self._labels2id[label]
    def id2label(self,label):
        return self._id2label[label]
    def is_valid_term(self,term):
        if term in self._labels2id:
            return True
        else:
            return False