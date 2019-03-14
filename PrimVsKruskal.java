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
public class PrimVsKruskal{

	/* PrimVsKruskal(G)
		Given an adjacency matrix for connected graph G, with no self-loops or parallel edges,
		determine if the minimum spanning tree of G found by Prim's algorithm is equal to 
		the minimum spanning tree of G found by Kruskal's algorithm.
		
		If G[i][j] == 0.0, there is no edge between vertex i and vertex j
		If G[i][j] > 0.0, there is an edge between vertices i and j, and the
		value of G[i][j] gives the weight of the edge.
		No entries of G will be negative.
	*/

	static boolean PrimVsKruskal(double[][] G){
		int n = G.length;
		boolean pvk = true;
		
		//Create our edgeweighted graph based off of the length of G
		EdgeWeightedGraph ew_graph = new EdgeWeightedGraph(n);

		//add all of the edges to our graph, except that we are testing at every step as an early failsafe
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (G[i][j] > 0) {
					Edge e = new Edge(i,j, (double)G[i][j]);
					ew_graph.addEdge(e);
				}
				//System.out.println(ew_graph.toString());

				PrimMST mst_prim = new PrimMST(ew_graph);
				KruskalMST mst_kruskal = new KruskalMST(ew_graph);

				//System.out.println(mst_kruskal.edges());

				//I wasn't really sure what the best way to compare the graphs was, so I just made them back into matrices and compared them
				double[][] prim_mat = new double[n][n];
				double[][] kruskal_mat = new double[n][n];

				//ugly unreadable code, but owell
				for (Edge edge : mst_prim.edges()) {
					prim_mat[edge.other(edge.either())][edge.either()] = prim_mat[edge.either()][edge.other(edge.either())] = edge.weight();
				}
				
				for (Edge edge : mst_kruskal.edges()) {
					kruskal_mat[edge.other(edge.either())][edge.either()] = kruskal_mat[edge.either()][edge.other(edge.either())] = edge.weight();
				}

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

				for (int a = 0; a < n; a++) {
					for (int b = a; b < n; b++) {
						if ((int)prim_mat[a][b] != (int)kruskal_mat[a][b]) {
							pvk = false;
							return pvk;
						}
					}
				}
			}
		}
		
		//System.out.println();
		return pvk;
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