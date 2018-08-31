import networkx as nx
import matplotlib.pyplot as plt

if __name__ == '__main__':
    target='/home/ndh0004/Downloads/qac_analysis_aug30/qac.edges'
    G = nx.read_edgelist(target)
    counter = 0
    plt.hist([len(con) for con in nx.connected_components(G) ], bins=30)
    plt.show()
    for X  in [con for con in nx.connected_components(G) if len(con) >4]:
        k = G.subgraph(X)
        print(X)
        nx.draw(k, with_labels=True,
                font_size=8)
        plt.show()

