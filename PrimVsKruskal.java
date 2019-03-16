/* PrimVsKruskal.java
   CSC 226 - Spring 2019
   Assignment 2 - Prim MST versus Kruskal MST Template
   
   The file includes the "import edu.princeton.cs.algs4.*;" so that yo can use
   any of the code in the algs4.jar file. You should be able to compile your program
   with the command
   
	javac -cp .;algs4.jar PrimVsKruskal.java
	
   To conveniently test the algorithm with a large input, create a text file
   containing a test graphs (in the format described below) and run
   the program with
   
	java -cp .;algs4.jar PrimVsKruskal file.txt
	
   where file.txt is replaced by the name of the text file.
   
   The input consists of a graph (as an adjacency matrix) in the following format:
   
    <number of vertices>
	<adjacency matrix row 1>
	...
	<adjacency matrix row n>
	
   Entry G[i][j] >= 0.0 of the adjacency matrix gives the weight (as type double) of the edge from 
   vertex i to vertex j (if G[i][j] is 0.0, then the edge does not exist).
   Note that since the graph is undirected, it is assumed that G[i][j]
   is always equal to G[j][i].
   R. Little - 03/07/2019
*/

import edu.princeton.cs.algs4.*;
import java.util.Scanner;
import java.io.File;

//Do not change the name of the PrimVsKruskal class
public class PrimVsKruskal {

	/* PrimVsKruskal(G)
		Given an adjacency matrix for connected graph G, with no self-loops or parallel edges,
		determine if the minimum spanning tree of G found by Prim's algorithm is equal to 
		the minimum spanning tree of G found by Kruskal's algorithm.
		
		If G[i][j] == 0.0, there is no edge between vertex i and vertex j
		If G[i][j] > 0.0, there is an edge between vertices i and j, and the
		value of G[i][j] gives the weight of the edge.
		No entries of G will be negative.
    */
    /*for (int a = 0; a < n; a++) {
        System.out.println();
        for (int b = a; b < n; b++) {
            System.out.print((int)prim_matrix[a][b] + " ");
        }
    }

    System.out.println();

    for (int a = 0; a < n; a++) {
        System.out.println();
        for (int b = a; b < n; b++) {
            System.out.print((int)kruskal_matrix[a][b] + " ");
        }
    }*/
    /*
    EdgeWeightedGraph ew_graph = new EdgeWeightedGraph(n);

    //add all of the edges to our graph, except that we are testing at every step as an early failsafe
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (G[i][j] > 0) {
                Edge e = new Edge(i,j, (double)G[i][j]);
                ew_graph.addEdge(e);
            }
        }
    }
    */

    static boolean PrimVsKruskal(double[][] G) {

        Prim_MST prim_mst = new Prim_MST(G);
        Kruskal_MST kruskal_mst = new Kruskal_MST(G);

        Edge k_edge = null;
        Edge p_edge = null;

        while (kruskal_mst.pq.isEmpty() == false && kruskal_mst.mst.size() < G.length - 1 || prim_mst.pq.isEmpty() == false) {
            if (kruskal_mst.pq.isEmpty() == false && kruskal_mst.mst.size() < G.length - 1) {
                try {
                    k_edge = kruskal_mst.greedy_kruskal(G);

                    if (is_a_cycle(k_edge, prim_mst.uf, prim_mst.edges()) == true) {
                        return false;
                    }
                }
                catch (NullPointerException e) {}
            }

            if (prim_mst.pq.isEmpty() == false) {
                try {
                    p_edge = prim_mst.greedy_prim(G);

                    if (is_a_cycle(p_edge, kruskal_mst.uf, kruskal_mst.edges()) == true) {
                        return false;
                    }
                }
                catch (NullPointerException e) {}
            }
        }
        return true;
    }

    public static boolean is_a_cycle(Edge e, UF uf, Iterable<Edge> mst) {
        boolean not_in_mst = true;
        for (Edge edge : mst) {
            if (check_edge(e, edge) == true) {
                not_in_mst = false;
            }
        }
        if (uf.connected(e.either(), e.other(e.either())) == true && not_in_mst == true) {
            return true;
        }
        else {
            return false;
        }
    }
    public static boolean check_edge(Edge edge_1, Edge edge_2) {
        if (edge_1.weight() == edge_2.weight()) {
            if (edge_1.either() == edge_2.either() && edge_1.other(edge_1.either()) == edge_2.other(edge_2.either())) {
                return true;
            }
            if (edge_1.other(edge_1.either()) == edge_2.either() && edge_1.either() == edge_2.other(edge_2.either())) {
                return true;
            }
        }
        return false;
    }
    static class Kruskal_MST {
        public static final double FLOATING_POINT_EPSILON = 1E-12;
        public MinPQ<Edge> pq;
        public Queue<Edge> mst = new Queue<Edge>(); // edges in MST
        public UF uf;

        public Kruskal_MST(double[][] G) {
            // more efficient to build heap by passing array of edges
            this.pq = new MinPQ<Edge>();

            for(int row = 0; row < G.length; row++) {
                for(int col = row + 1; col < G.length; col++) {
                    if (G[row][col] > 0) {
                        Edge temp_edge = new Edge(row, col, G[row][col]);
                        pq.insert(temp_edge);
                    }
                }
            }
            this.uf = new UF(G.length);
        }
        public Edge greedy_kruskal(double[][] G){
            Edge edge = pq.delMin();
            int v = edge.either();
            int w = edge.other(v);

            if (uf.connected(v, w) == false) { // v-w does not create a cycle
                uf.union(v, w); // merge v and w components
                mst.enqueue(edge); // add edge e to mst
                return edge;
            }
            return null;
        }
        public Iterable<Edge> edges() {
            return mst;
        }
    }
    static class Prim_MST {
        public static final double FLOATING_POINT_EPSILON = 1E-12;

        public Edge[] edgeTo;        // edgeTo[v] = shortest edge from tree vertex to non-tree vertex
        public double[] distTo;      // distTo[v] = weight of shortest such edge
        public boolean[] marked;     // marked[v] = true if v on tree, false otherwise
        public IndexMinPQ<Double> pq;
        public UF uf;
        public static int ind = 0;

        public Prim_MST(double G[][]) {
            edgeTo = new Edge[G.length];
            distTo = new double[G.length];
            marked = new boolean[G.length];
            this.pq = new IndexMinPQ<Double>(G.length);

            for (int v = 0; v < G.length; v++) {
                distTo[v] = Double.POSITIVE_INFINITY;
            }

            distTo[ind] = 0.0;
            pq.insert(ind, distTo[ind]);
            this.uf = new UF(G.length);
        }
        // run Prim's algorithm in graph G, starting from vertex s
        public Edge greedy_prim(double G[][]) {
            int v = pq.delMin();
            Edge edge = scan(G, v);
            if (edge != null) {
                uf.union(v, edge.other(v));
            }
            return edge;
        }

        // scan vertex v
        public Edge scan(double G[][], int v) {
            marked[v] = true;
            for (int u = 0; u < G.length; u++) {
                if(G[v][u] > 0) {
                    Edge e = new Edge(v, u, G[v][u]);
                    int w = e.other(v);

                    if (marked[w] == true) {
                        continue;         // v-w is obsolete edge
                    }
                    if (e.weight() < distTo[w]) {
                        distTo[w] = e.weight();
                        edgeTo[w] = e;

                        if (pq.contains(w) == true) {
                            pq.decreaseKey(w, distTo[w]);
                        }
                        else {
                            pq.insert(w, distTo[w]);
                        }
                    }
                }
            }
            return edgeTo[v];
        }
        public Iterable<Edge> edges() {
            Queue<Edge> mst = new Queue<Edge>();
            for (int v = 0; v < edgeTo.length; v++) {
                Edge e = edgeTo[v];
                
                if (e != null) {
                    mst.enqueue(e);
                }
            }
            return mst;
        }
        // check optimality conditions (takes time proportional to E V lg* V)   
    }
	/* main()
	   Contains code to test the PrimVsKruskal function. You may modify the
	   testing code if needed, but nothing in this function will be considered
	   during marking, and the testing process used for marking will not
	   execute any of the code below. 
	*/
	public static void main(String[] args) {
		Scanner s;
		if (args.length > 0){
			try{
				s = new Scanner(new File(args[0]));
			} catch(java.io.FileNotFoundException e){
				System.out.printf("Unable to open %s\n",args[0]);
				return;
			}
			System.out.printf("Reading input values from %s.\n",args[0]);
		}else{
			s = new Scanner(System.in);
			System.out.printf("Reading input values from stdin.\n");
		}
		
		int n = s.nextInt();
		double[][] G = new double[n][n];
		int valuesRead = 0;
		for (int i = 0; i < n && s.hasNextDouble(); i++){
			for (int j = 0; j < n && s.hasNextDouble(); j++){
				G[i][j] = s.nextDouble();
				if (i == j && G[i][j] != 0.0) {
					System.out.printf("Adjacency matrix contains self-loops.\n");
					return;
				}
				if (G[i][j] < 0.0) {
					System.out.printf("Adjacency matrix contains negative values.\n");
					return;
				}
				if (j < i && G[i][j] != G[j][i]) {
					System.out.printf("Adjacency matrix is not symmetric.\n");
					return;
				}
				valuesRead++;
			}
		}
		
		if (valuesRead < n*n){
			System.out.printf("Adjacency matrix for the graph contains too few values.\n");
			return;
		}	
		
        boolean pvk = PrimVsKruskal(G);
        System.out.printf("Does Prim MST = Kruskal MST? %b\n", pvk);
    }
}