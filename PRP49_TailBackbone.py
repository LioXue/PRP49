class TailBackbone:
    def __init__(self, start = None, end = None):
        self.start = start
        self.end = end

    class _N:
        def __init__(self, position):
            self.position = position
            self.next = None
            self.hydrogen = None

    class _CA:
        def __init__(self, position):
            self.position = position
            self.next = None
            self.sidechain = [position.atom_info.sidechain] if type(position.atom_info.sidechain) is not list else position.atom_info.sidechain

    class _C:
        def __init__(self, position):
            self.position = position
            self.next = None
            self.oxygen = []

    class NodePosition:
        def __init__(self, node):
            self.node = node

    def make_NodePosition(self, node):
        if isinstance(node, TailBackbone._N) or isinstance(node, TailBackbone._CA) or isinstance(node, TailBackbone._C):
            return self.NodePosition(node)
        else:
            return None

    def add(self,dict):
        if self.start is None:
            ca = self._CA(dict["ca"])
            c = self._C(dict["c"])
            self.start = ca
            ca.next = c
            c.oxygen.append(dict["o"])
            if dict["oxt"]:
                c.oxygen.append(dict["oxt"])
            self.end = c
        else:
            n = self._N(dict["n"])
            n.hydrogen = dict["h"]
            ca = self._CA(dict["ca"])
            self.end.next = n
            n.next = ca
            if not dict["c"].atom_info.root:
                c = self._C(dict["c"])
                c.oxygen.append(dict["o"])
                if dict["oxt"]:
                    c.oxygen.append(dict["oxt"])
                ca.next = c
                self.end = c
            else:
                self.end = ca
    def get_start(self):
        return self.make_NodePosition(self.start)

    def get_next(self, NodePosition):
        return self.make_NodePosition(NodePosition.node.next)

    def get_next_branch(self, NodePosition = None):
        if NodePosition is None:
            return self.start.position
        elif isinstance(NodePosition.node.next, TailBackbone._CA):
            return NodePosition.node.next.position
        elif isinstance(NodePosition.node.next, TailBackbone._C):
            if NodePosition.node.next.next and NodePosition.node.next.next.hydrogen is None: #为PRO打补丁：
                return NodePosition.node.next.position, NodePosition.node.next.oxygen[0], NodePosition.node.next.next.position,
            elif NodePosition.node.next.next and len(NodePosition.node.next.oxygen) == 1:
                return NodePosition.node.next.position, NodePosition.node.next.oxygen[0], NodePosition.node.next.next.position, NodePosition.node.next.next.hydrogen
            elif not NodePosition.node.next.next and len(NodePosition.node.next.oxygen) > 1:
                return NodePosition.node.next.position, NodePosition.node.next.oxygen[0], NodePosition.node.next.oxygen[1]
            else:
                return NodePosition.node.next.position, NodePosition.node.next.oxygen[0]
        else:
            return None

    def __iter__(self):
        i = self.make_NodePosition(self.start)
        yield i.node.position
        while i:
            yield self.get_next_branch(i)
            i = self.get_next(i)




