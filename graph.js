//Directed graphs
//Directed Graph using an adjacency list:
// Define the Graph class
class Graph {
    constructor() {
      // Create an empty adjacency list using a Map
      this.adjList = new Map();
    }
    
    // Add a vertex to the graph
    addVertex(vertex) {
      // Check if the vertex is already in the graph
      if (!this.adjList.has(vertex)) {
        // If not, add the vertex and initialize an empty array for its neighbors
        this.adjList.set(vertex, []);
      }
    }
    
    // Add a directed edge between two vertices
    addEdge(start, end) {
      // If the starting vertex is not already in the graph, add it with an empty array for its neighbors
      if (!this.adjList.has(start)) {
        this.adjList.set(start, []);
      }
      // Add the ending vertex to the list of neighbors for the starting vertex
      this.adjList.get(start).push(end);
    }
  }
  
  // Create a new instance of the Graph class
  const graph1 = new Graph();
  
  // Add vertices to the graph
  graph1.addVertex('A');
  graph1.addVertex('B');
  graph1.addVertex('C');
  graph1.addVertex('D');
  
  // Add directed edges between the vertices
  graph1.addEdge('A', 'B');
  graph1.addEdge('B', 'C');
  graph1.addEdge('C', 'D');
  graph1.addEdge('D', 'A');
  
/* Notes
In this example, the Graph class has an adjList property that is a Map object.
The addVertex() method adds a vertex to the graph by creating a new entry in the adjList map.
The addEdge() method adds a directed edge from one vertex to another by appending the target vertex to the list of adjacent vertices for the starting vertex.
The code creates a new directed graph with vertices labeled 'A', 'B', 'C', and 'D'. 
The addEdge() method is used to create directed edges between the vertices.
*/

//_______________________________________________________________________________________

//Undirected Graphs

/* Notes
An undirected graph is a type of graph where the edges between the vertices have no orientation.
In other words, the edges are bi-directional, and you can travel from one vertex to another in both directions.
In contrast to directed graphs, there is no concept of in-degree and out-degree in undirected graphs.
Instead, each vertex in an undirected graph has a degree, which represents the number of edges connected to it.
Undirected graphs are commonly used to model pairwise relationships between objects.
Like social networks, road networks, and electrical circuits.
*/

//Undirected Graph using an adjacency list representation:
class Graph {
    constructor() {
      this.adjList = new Map();
    }
    
    addVertex(vertex) {
      if (!this.adjList.has(vertex)) {
        this.adjList.set(vertex, []);
      }
    }
    
    addEdge(v1, v2) {
      // Add v1 as a neighbor of v2
      this.adjList.get(v1).push(v2);
      // Add v2 as a neighbor of v1 (because the graph is undirected)
      this.adjList.get(v2).push(v1);
    }
  }
  
  // Create a new graph
  const graph2 = new Graph();
  
  // Add vertices to the graph
  graph2.addVertex('A');
  graph2.addVertex('B');
  graph2.addVertex('C');
  graph2.addVertex('D');
  graph2.addVertex('E');
  
  // Add edges between vertices
  graph2.addEdge('A', 'B');
  graph2.addEdge('A', 'D');
  graph2.addEdge('B', 'C');
  graph2.addEdge('B', 'D');
  graph2.addEdge('C', 'E');
  
  // Print the adjacency list of the graph
  console.log(graph2.adjList);
/* Notes
In this example, we use a Map data structure to represent the adjacency list of the graph.
The addVertex method adds a new vertex to the graph.
The addEdge method adds an undirected edge between two vertices.
*/

//_______________________________________________________________________________________

//Graph Traversal

//Breadth First Traversal
/* Notes
Otherwise known as Breadth First Search.
Can simply be defined as traversing through the graph horizontally, visiting each node at each row in the graph.
The algorithm for BFS is "vertex-based", meaning it will check each adjacent vertex within the data structure first.
It is considered slower than DFS, but better for searching nodes closer to a given node source position.
*/
//Breadth First Traversal Example:
class Graph {
    constructor() {
      this.adjList = new Map();
    }
    
    addVertex(vertex) {
      if (!this.adjList.has(vertex)) {
        this.adjList.set(vertex, []);
      }
    }
    
    addEdge(start, end) {
      if (!this.adjList.has(start)) {
        this.adjList.set(start, []);
      }
      this.adjList.get(start).push(end);
    }
    
    bfs(start) {
      let visited = new Set();
      let queue = [start];
      visited.add(start);
  
      while (queue.length > 0) {
        let currentVertex = queue.shift();
        console.log(currentVertex);
  
        let neighbors = this.adjList.get(currentVertex);
  
        for (let neighbor of neighbors) {
          if (!visited.has(neighbor)) {
            visited.add(neighbor);
            queue.push(neighbor);
          }
        }
      }
    }
  }
  
  // Create a new graph
  const graph3 = new Graph();
  
  // Add vertices and edges to the graph
  graph3.addVertex('A');
  graph3.addVertex('B');
  graph3.addVertex('C');
  graph3.addVertex('D');
  graph3.addEdge('A', 'B');
  graph3.addEdge('A', 'C');
  graph3.addEdge('B', 'D');
  graph3.addEdge('C', 'D');
  
  // Traverse the graph using BFS
  graph3.bfs('A');
/* Notes
In this example, we define a Graph class with methods to add vertices and edges.
In this example, we simply log each visited vertex to the console.
You could modify the code to perform other operations on the vertices, such as finding a specific vertex or determining the shortest path between two vertices.
*/

//Depth First Traversal
/* Notes
Also known as Depth First Search.
Can be defined as traversing through the graph column by column.
The algorithm for DFS is "edges-based", checking if there are any unvisited vertices connected to the node.
*/
//Depth First Travel example:
class Graph {
    constructor() {
      this.adjList = new Map();
    }
    
    addVertex(vertex) {
      if (!this.adjList.has(vertex)) {
        this.adjList.set(vertex, []);
      }
    }
    
    addEdge(start, end) {
      if (!this.adjList.has(start)) {
        this.adjList.set(start, []);
      }
      this.adjList.get(start).push(end);
    }
  
    // DFS algorithm
    dfs(start) {
      const visited = new Set(); // Set to keep track of visited nodes
      this._dfs(start, visited); // Call recursive DFS helper function
    }
  
    _dfs(node, visited) {
      visited.add(node); // Mark the current node as visited
      console.log(node); // Output the visited node to the console
  
      // Recursively visit all adjacent nodes that haven't been visited yet
      for (const neighbor of this.adjList.get(node)) {
        if (!visited.has(neighbor)) {
          this._dfs(neighbor, visited);
        }
      }
    }
  }
  
  // Create a new graph and add vertices/edges
  const graph4 = new Graph();
  graph4.addVertex('A');
  graph4.addVertex('B');
  graph4.addVertex('C');
  graph4.addVertex('D');
  graph4.addEdge('A', 'B');
  graph4.addEdge('B', 'C');
  graph4.addEdge('C', 'D');
  graph4.addEdge('D', 'A');
  
  // Perform DFS on the graph starting from vertex 'A'
  graph4.dfs('A');  
/* Notes
In this implementation, the DFS algorithm is defined as a method of the Graph class, taking the starting vertex as its only parameter.
The method creates a Set to keep track of visited nodes, and then calls a recursive helper function _dfs with the starting node and the visited set.
The _dfs function takes two parameters: the current node being visited and the set of visited nodes.
It marks the current node as visited and outputs it to the console, then recursively visits all adjacent nodes that haven't already been visited.
This is done by looping over the neighbors of the current node and checking if they're in the visited set.
Finally, we create a new instance of the Graph class and add vertices/edges to it.
We then call the dfs method on the graph starting from vertex 'A', which should output the vertices in depth-first order: 'A', 'B', 'C', 'D'.
*/

//Complete Graphs
/* Notes
A complete graph is an undirected graph in which each vertex is connected to every other vertex.
In other words, there are no isolated vertices, and every vertex has degree n-1, where n is the total number of vertices in the graph.
Complete graphs are denoted by the symbol Kn, where n is the number of vertices.
For example, K3 is a complete graph with three vertices, where each vertex is connected to the other two vertices.
Similarly, K4 is a complete graph with four vertices, where each vertex is connected to the other three vertices.
The number of edges in a complete graph Kn can be calculated using the formula:
E = (n*(n-1))/2
E is the number of edges in the graph.
This formula is derived from the fact that each vertex is connected to every other vertex.
So the total number of edges is the sum of the degrees of all the vertices divided by two.
Complete graphs are also used in various applications, such as network routing, social network analysis, and clustering algorithms.
*/
//Example of a complete graph:
class CompleteGraph {
    constructor(numNodes) {
      // Create an empty adjacency matrix of size numNodes x numNodes
      this.adjMatrix = new Array(numNodes).fill(null).map(() => new Array(numNodes).fill(0));
      
      // Fill the adjacency matrix with edges between every pair of nodes
      for (let i = 0; i < numNodes; i++) {
        for (let j = i+1; j < numNodes; j++) {
          // Add an edge between node i and node j with weight 1
          this.adjMatrix[i][j] = 1;
          this.adjMatrix[j][i] = 1;
        }
      }
    }
  }
  
  // Example usage: create a complete graph with 4 nodes
  const graph5 = new CompleteGraph(4);
  
  // The resulting adjacency matrix looks like this:
  // [
  //   [0, 1, 1, 1],
  //   [1, 0, 1, 1],
  //   [1, 1, 0, 1],
  //   [1, 1, 1, 0]
  // ]
  // This represents a graph with nodes 0, 1, 2, and 3, where every pair of nodes is connected by an edge with weight 1.  
/* Notes
In this example, we define a CompleteGraph class that takes a number of nodes as its argument.
It creates an empty adjMatrix that will hold the edges between nodes.
Fills it with edges between every pair of nodes using two nested loops.
The resulting adjMatrix represents a complete graph with every pair of nodes connected by an edge of weight 1.
We then create an instance of the CompleteGraph class with 4 nodes.
The resulting adjMatrix represents a complete graph with 4 nodes where every pair of nodes is connected by an edge of weight 1.
*/

//Weighted Graphs
/* Notes
A weighted graph is a graph where each edge has a weight or cost associated with it.
This weight can represent any kind of metric, such as distance, time, or cost.
Weighted graphs are used to model a wide variety of problems, such as finding the shortest path between two nodes or minimizing the cost of a network.
A weighted graph can be represented using an adjacency matrix or an adjacency list, just like an unweighted graph.
The difference is that each entry in the matrix or each node in the list includes information about the weight of the edge.
To perform graph algorithms on a weighted graph, you need to take into account the weight of each edge.
For example, to find the shortest path between two nodes in a weighted graph, you can use algorithms such as Dijkstra's algorithm or the Bellman-Ford algorithm.
These algorithms take into account the weight of each edge to find the shortest path between two nodes.
*/
//Example of a weighted graph:
// Define a class to represent a weighted graph
class WeightedGraph {
    constructor() {
      // Initialize an empty adjacency list to store the vertices and their edges
      this.adjList = new Map();
    }
  
    // Method to add a vertex to the graph
    addVertex(vertex) {
      // If the vertex is not already in the graph, add it to the adjacency list with an empty array for its edges
      if (!this.adjList.has(vertex)) {
        this.adjList.set(vertex, []);
      }
    }
  
    // Method to add an edge with a weight between two vertices
    addEdge(start, end, weight) {
      // If either the start or end vertex is not in the graph, add it
      if (!this.adjList.has(start)) {
        this.adjList.set(start, []);
      }
      if (!this.adjList.has(end)) {
        this.adjList.set(end, []);
      }
      // Add the edge to the adjacency list for both the start and end vertices, along with its weight
      this.adjList.get(start).push({ node: end, weight });
      this.adjList.get(end).push({ node: start, weight });
    }
  }
  
  // Create a new weighted graph
  const graph6 = new WeightedGraph();
  
  // Add some vertices to the graph
  graph6.addVertex("A");
  graph6.addVertex("B");
  graph6.addVertex("C");
  graph6.addVertex("D");
  
  // Add some edges with weights to the graph
  graph6.addEdge("A", "B", 5);
  graph6.addEdge("B", "C", 3);
  graph6.addEdge("C", "D", 2);
  graph6.addEdge("D", "A", 1);
/* Notes
In this example, we define a WeightedGraph class that has methods to add vertices and edges with weights.
The addVertex method is similar to the one in the undirected graph example.
The addEdge method takes an additional parameter for the weight of the edge.
*/

//Bipartite Graphs
/* Notes
A bipartite graph is a graph whose vertices can be divided into two independent sets, such that no two vertices within the same set are adjacent.
In other words, it is a graph whose vertices can be partitioned into two sets, and all edges connect vertices from different sets.
A bipartite graph can also be referred to as a 2-partite graph, a bigraph, or a bipartite graph.
Bipartite graphs have many applications, including modeling matching problems, such as matching people with jobs, students with schools, or advertisers with ad slots.
To determine if a graph is bipartite, you can use graph traversal algorithms, such as BFS or DFS, along with a color-coding method.
In this method, you start at any vertex and assign it a color (either red or blue).
Then, you traverse the graph and assign the opposite color to all of its neighbors.
If you encounter a vertex that has already been colored with the same color as its neighbor, then the graph is not bipartite.
Otherwise, continue the traversal until all vertices have been colored.
If the graph is bipartite, you can also determine its two partitions by taking all the vertices of one color and placing them in one set, and all the vertices of the other color in the other set.
Bipartite graphs can be represented in many ways, including adjacency matrix or adjacency list.
*/
//Example of a bipartite graph:
// A bipartite graph is a graph whose vertices can be divided into two independent sets
// such that every edge connects a vertex from one set to a vertex from the other set.

// Here is an example of a bipartite graph with two sets of vertices: {1, 2, 3} and {4, 5, 6}
// Each edge connects a vertex from one set to a vertex from the other set, but there are no edges
// that connect vertices within the same set.

// We can represent the graph using an adjacency list.
const graph7 = new Map([
    [1, [4, 5]],
    [2, [4, 6]],
    [3, [5, 6]],
    [4, [1, 2]],
    [5, [1, 3]],
    [6, [2, 3]]
  ]);
  
  // To check if a graph is bipartite, we can perform a depth-first search (DFS) starting at each vertex.
  // During the DFS, we will color each vertex either red or blue, and we will make sure that no two
  // adjacent vertices have the same color. If we are able to color all vertices without violating this
  // rule, then the graph is bipartite.
  
  function isBipartite(graph7) {
    // Create a color map to keep track of the color of each vertex.
    const colors = new Map();
    
    // Perform a DFS starting at each vertex.
    for (let vertex of graph7.keys()) {
      if (!colors.has(vertex)) {
        if (!dfs(vertex, graph7, colors, 'red')) {
          return false;
        }
      }
    }
    
    // If we have successfully colored all vertices without violating the bipartite property,
    // then the graph is bipartite.
    return true;
  }
  
  function dfs(vertex, graph7, colors, color) {
    // Set the color of the current vertex.
    colors.set(vertex, color);
    
    // Loop through all neighbors of the current vertex.
    for (let neighbor of graph7.get(vertex)) {
      // If the neighbor has not been colored yet, color it the opposite color of the current vertex.
      if (!colors.has(neighbor)) {
        if (!dfs(neighbor, graph7, colors, color === 'red' ? 'blue' : 'red')) {
          return false;
        }
      }
      // If the neighbor has already been colored, make sure its color is different from the current vertex.
      else if (colors.get(neighbor) === color) {
        return false;
      }
    }
    
    // If we have successfully colored all neighbors without violating the bipartite property,
    // then the current vertex can be part of a bipartition.
    return true;
  }
  
  // Now we can test the isBipartite function with our example graph.
  console.log(isBipartite(graph7)); // true
   
//

//Planar Graphs
/* Notes
A planar graph is a graph that can be drawn on a plane without any edges crossing each other.
Planar graphs have many interesting properties and are used in various fields, such as computer science, graph theory, and geography.
One of the key properties of planar graphs is that they have a maximum number of edges for a given number of vertices.
This is known as the Euler's formula, which states that for any planar graph with V vertices, E edges, and F faces (regions bounded by edges):
V - E + F = 2
Another interesting property of planar graphs is that they can be colored using only four colors in such a way that no two adjacent vertices have the same color. 
This is known as the Four Color Theorem and was first proved in 1976 by Kenneth Appel and Wolfgang Haken.
Planar graphs also have applications in various algorithms, such as in finding the shortest path between two points in a 2D space, and in routing algorithms for computer networks.
It's worth noting that not all graphs are planar.
For example, the complete graph K5 and the complete bipartite graph K3,3 are non-planar.
However, many graphs encountered in practice are planar or can be made planar by adding or removing edges.
In summary, planar graphs are a fundamental concept in graph theory and have many interesting properties and applications in various fields.
*/

//Dijkstra’s Shortest Path Algorithm 
/* Notes

*/
//The shortest path between two points in a road network using Dijkstra's algorithm:
function shortestPath(network, start, end) {
  // Initialize distances and visited nodes
  const distances = {};
  const visited = new Set();
  const unvisited = new Set(Object.keys(network));

  for (const node of unvisited) {
    distances[node] = node === start ? 0 : Infinity;
  }

  // Find the shortest path
  while (unvisited.size > 0) {
    // Find the unvisited node with smallest distance
    const current = [...unvisited].reduce((minNode, node) => (
      distances[node] < distances[minNode] ? node : minNode
    ), null);

    // Stop if we've reached the end node
    if (current === end) {
      break;
    }

    // Visit current node and update distances of its neighbors
    visited.add(current);
    unvisited.delete(current);

    for (const neighbor of Object.keys(network[current])) {
      const distance = network[current][neighbor];
      const totalDistance = distances[current] + distance;
      if (totalDistance < distances[neighbor]) {
        distances[neighbor] = totalDistance;
      }
    }
  }

  // Reconstruct the shortest path
  const path = [end];
  let node = end;
  while (node !== start) {
    for (const neighbor of Object.keys(network[node])) {
      if (distances[neighbor] === distances[node] - network[node][neighbor]) {
        path.unshift(neighbor);
        node = neighbor;
        break;
      }
    }
  }

  return path;
}
/* Notes

*/

//Dijkstra's algorithm using a priority queue:
// This function implements Dijkstra's shortest path algorithm to find the shortest path
// between two nodes in a weighted graph represented by an object called "network".

function dijkstraShortestPath(network, start, end) {
  // Initialize distances and visited nodes
  // "distances" is an object that stores the shortest distances from the start node to each node in the network
  // "visited" is a set that stores the nodes that have been visited
  // "unvisited" is a priority queue that stores the unvisited nodes and their distances from the start node
  const distances = {};
  const visited = new Set();
  const unvisited = new PriorityQueue();

  // Set the distance of the start node to 0 and enqueue it to the priority queue
  // Set the distance of all other nodes to Infinity and enqueue them to the priority queue
  for (const node of Object.keys(network)) {
    distances[node] = node === start ? 0 : Infinity;
    unvisited.enqueue(node, distances[node]);
  }

  // Find the shortest path
  while (!unvisited.isEmpty()) {
    // Find the unvisited node with smallest distance
    // Dequeue the node with the smallest distance from the priority queue
    const current = unvisited.dequeue().element;

    // Stop if we've reached the end node
    if (current === end) {
      break;
    }

    // Visit current node and update distances of its neighbors
    visited.add(current);

    // For each neighbor of the current node, calculate the total distance to that neighbor
    // and update its distance if the new distance is shorter than the current distance
    // Enqueue the neighbor to the priority queue if its distance is updated
    for (const neighbor of Object.keys(network[current])) {
      if (visited.has(neighbor)) {
        continue;
      }
      const distance = network[current][neighbor];
      const totalDistance = distances[current] + distance;
      if (totalDistance < distances[neighbor]) {
        distances[neighbor] = totalDistance;
        unvisited.enqueue(neighbor, totalDistance);
      }
    }
  }

  // Reconstruct the shortest path
  // Start from the end node and follow the path by looking at the distances of its neighbors
  const path = [end];
  let node = end;
  while (node !== start) {
    for (const neighbor of Object.keys(network[node])) {
      if (distances[neighbor] === distances[node] - network[node][neighbor]) {
        path.unshift(neighbor);
        node = neighbor;
        break;
      }
    }
  }

  // Return the shortest path
  return path;
}

// This is a helper class for implementing the priority queue used in Dijkstra's algorithm
class PriorityQueue {
  constructor() {
    this.elements = [];
  }

  // Enqueue an element with a priority into the priority queue
  enqueue(element, priority) {
    const item = { element, priority };
    let added = false;
    for (let i = 0; i < this.elements.length; i++) {
      if (item.priority < this.elements[i].priority) {
        this.elements.splice(i, 0, item);
        added = true;
        break;
      }
    }
    if (!added) {
      this.elements.push(item);
    }
  }

  // Dequeue the element with the highest priority from the priority queue
  dequeue() {
    return this.elements.shift();
  }

  // Check if the priority queue is empty
  isEmpty() {
    return this.elements.length === 0;
  }
}
/* Notes

*/

//Bellman-Ford Algorithm
/* Notes

*/
//Bellman-Ford Algorithm example solution:
class Graph {
  constructor(numVertices) {
    this.numVertices = numVertices; // number of vertices in the graph
    this.adjList = new Map(); // adjacency list to store edges and their weights
  }

  addEdge(u, v, weight) {
    // add an edge with weight between vertices u and v
    if (!this.adjList.has(u)) {
      this.adjList.set(u, []);
    }
    this.adjList.get(u).push({ v, weight });
  }

  bellmanFord(start) {
    let dist = new Array(this.numVertices).fill(Number.MAX_SAFE_INTEGER); // initialize all distances to Infinity
    dist[start] = 0; // distance to starting vertex is 0

    for (let i = 0; i < this.numVertices - 1; i++) {
      // loop through all vertices except the last one
      for (let [u, edges] of this.adjList) {
        // loop through each vertex and its edges
        for (let { v, weight } of edges) {
          // loop through each edge of the vertex
          if (dist[u] !== Number.MAX_SAFE_INTEGER && dist[v] > dist[u] + weight) {
            // if the distance to u is not Infinity and the distance to v is greater than the distance to u plus the weight of the edge,
            dist[v] = dist[u] + weight; // update the distance to v
          }
        }
      }
    }

    for (let [u, edges] of this.adjList) {
      // loop through each vertex and its edges
      for (let { v, weight } of edges) {
        // loop through each edge of the vertex
        if (dist[u] !== Number.MAX_SAFE_INTEGER && dist[v] > dist[u] + weight) {
          // if the distance to u is not Infinity and the distance to v is greater than the distance to u plus the weight of the edge,
          console.log("Graph contains negative weight cycle"); // there is a negative weight cycle in the graph
          return;
        }
      }
    }

    console.log(dist); // print the distances from the starting vertex to all other vertices
  }
}

// Example usage
let g = new Graph(5);
g.addEdge(0, 1, 6);
g.addEdge(0, 3, 7);
g.addEdge(1, 2, 5);
g.addEdge(1, 3, 8);
g.addEdge(1, 4, -4);
g.addEdge(2, 1, -2);
g.addEdge(3, 2, -3);
g.addEdge(3, 4, 9);
g.addEdge(4, 0, 2);
g.addEdge(4, 2, 7);
g.bellmanFord(0); // find the shortest paths from vertex 0 to all other vertices
/* Notes

*/

//Floyd-Warshall Algorithm
/* Notes

*/
//Example implementation of the Floyd-Warshall algorithm:
function floydWarshall(graph) {
  // Get the number of vertices in the graph
  const n = graph.length;

  // Create a 2D array of distances initialized to Infinity
  const dist = new Array(n).fill().map(() => new Array(n).fill(Infinity));

  // Initialize the distance array with the edge weights of the graph
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      // If there is an edge between vertices i and j, set the distance to the weight of the edge
      if (graph[i][j]) {
        dist[i][j] = graph[i][j];
      }
    }
    // Set the distance from a vertex to itself to 0
    dist[i][i] = 0;
  }

  // Perform the Floyd-Warshall algorithm
  for (let k = 0; k < n; k++) {
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        // If the path from i to j through vertex k is shorter than the current path from i to j,
        // update the distance between i and j to the new, shorter distance
        if (dist[i][k] + dist[k][j] < dist[i][j]) {
          dist[i][j] = dist[i][k] + dist[k][j];
        }
      }
    }
  }

  // Return the distance array
  return dist;
}
/* Notes

*/