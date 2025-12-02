import random
import networkx as nx
from typing import Iterable

class Constants:
    N0 = 10
    M0 = 3
    T_MAX = 10 ** 5

class GraphGenerator():

    name: str
    dir_name: str

    @staticmethod
    def generate(        
        n0: int = Constants.N0,
        m0: int = Constants.M0,
        t_max: int = Constants.T_MAX,
        *,
        seed: int | None = None
        ) -> Iterable[tuple[int, nx.Graph]]:
        """
        Paràmetres
        ----------
        n0 : int
            Nombre inicial de vèrtexs (al temps t = 0).
        m0 : int
            Nombre d'arestes que afegeix cada nou vèrtex (m0 >= 1).
        t_max : int
            Nombre de passos de temps (nombre de vèrtexs nous que afegirem).
            El nombre final de vèrtexs serà n0 + t_max.
        seed : int | None
            Llavor per al generador aleatori (opcional, per reproductibilitat).

        Retorna (t, G) a cada iteració, on t és el temps actual i G és el graf generat fins al moment.
        """
        raise NotImplementedError()
    
    def init(n0: int, m0: int, t_max:int, seed: int | None = None) -> None:
        random.seed(seed)
        if m0 < 1:    raise ValueError("m0 ha de ser >= 1")
        if n0 < m0:   raise ValueError("Cal que n0 >= m0 per poder connectar el nou vèrtex a m0 vèrtexs diferents.")
        if t_max < 0: raise ValueError("t_max ha de ser >= 0")

class BarabasiAlbert(GraphGenerator):
    name = "Barabási-Albert Preferential Attachment"
    dir_name = "barabasi_albert_n" + Constants.N0 + "_m" + Constants.M0
    @staticmethod
    def generate(
        n0: int = Constants.N0,
        m0: int = Constants.M0,
        t_max: int = Constants.T_MAX,
        *,
        seed: int | None = None
    ) -> Iterable[tuple[int, nx.Graph]]:

        GraphGenerator.init(n0, m0, t_max, seed)

        G = nx.Graph()
        G.add_nodes_from(range(-n0 + 1, 1)) # Vèrtexs inicials: -n0+1, ..., 0

        # Creem una anella (cicle)
        stubs = []
        for u in range(-n0 + 1, 0):
            G.add_edge(u, u + 1)
            stubs.append(u)
            stubs.append(u + 1)
        G.add_edge(0, -n0 + 1)
        stubs.append(0)
        stubs.append(-n0 + 1)

        for t in range(1, t_max + 1):
            new_v = t

            if len(stubs) >= m0:
                targets = []
                while len(targets) < m0:
                    u = random.choice(stubs)
                    if u not in targets:
                        targets.append(u)
            else: # Això no hauria de passar
                targets = random.sample(stubs + list(G.nodes), k=m0)

            G.add_node(new_v)

            for u in targets:
                G.add_edge(new_v, u)
                stubs.append(new_v)
                stubs.append(u)

            yield (t, G)

class BA_RandomAttachment(GraphGenerator):
    name = "Barabási-Albert Random Attachment"
    dir_name = "ba_random_attachment_n" + Constants.N0 + "_m" + Constants.M0
    @staticmethod
    def generate(
        n0: int = Constants.N0,
        m0: int = Constants.M0,
        t_max: int = Constants.T_MAX,
        *,
        seed: int | None = None
    ) -> Iterable[tuple[int, nx.Graph]]:

        GraphGenerator.init(n0, m0, t_max, seed)

        G = nx.Graph()
        G.add_nodes_from(range(-n0 + 1, 1)) # Vèrtexs inicials: -n0+1, ..., 0

        # Creem una anella (cicle) per coherència, tot i que aquí no faria falta 
        for u in range(-n0 + 1, 0):
            G.add_edge(u, u + 1)
        G.add_edge(0, -n0 + 1)

        for t in range(1, t_max + 1):
            new_v = t
            targets = random.sample(list(G.nodes), k=m0)
            G.add_node(new_v)
            for u in targets:
                G.add_edge(new_v, u)
            yield (t, G)

class BA_StaticNodes(GraphGenerator):
    name = "Barabási-Albert Static Nodes"
    dir_name = "ba_static_nodes_n" + Constants.N0 + "_m" + Constants.M0

    @staticmethod
    def generate(
        n0: int = 10000,
        m0: int = Constants.M0,
        t_max: int = Constants.T_MAX,
        *,
        seed: int | None = None
    ) -> Iterable[tuple[int, nx.Graph]]:

        GraphGenerator.init(n0, m0, t_max, seed)

        G = nx.Graph()
        G.add_nodes_from(range(1, n0 + 1)) # Vèrtexs: 1, ..., n0

        # Creem una anella (cicle)
        stubs = []
        for u in range(1, n0):
            G.add_edge(u, u + 1)
            stubs.append(u)
            stubs.append(u + 1)
        G.add_edge(n0, 1)
        stubs.append(n0)
        stubs.append(1)

        for t in range(1, t_max + 1):
            v = random.randrange(1, n0 + 1)
            if len(stubs) >= m0:
                targets = []
                while len(targets) < m0:
                    u = random.choice(stubs)
                    if u not in targets:
                        targets.append(u)
            else: # Això no hauria de passar
                targets = random.sample(stubs + list(G.nodes), k=m0)
            for u in targets:
                G.add_edge(v, u)
                stubs.append(v)
                stubs.append(u)
            yield (t, G)