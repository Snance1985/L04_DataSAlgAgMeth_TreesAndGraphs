//Trees

//Binary Trees
/* Notes

*/
//Example implementation of a binary tree:
// This class represents a node in a binary tree
class TreeNode {
    constructor(value, left = null, right = null) {
      // The value of the node
      this.value = value;
      // The left child of the node
      this.left = left;
      // The right child of the node
      this.right = right;
    }
  }
  
  // This class represents a binary tree
  class BinaryTree {
    constructor(root = null) {
      // The root node of the tree
      this.root = root;
    }
  
    // Inserts a value into the tree
    insert(value) {
      // If the tree is empty, create a new root node
      if (!this.root) {
        this.root = new TreeNode(value);
        return;
      }
  
      // Traverse the tree to find the appropriate position for the new node
      let current = this.root;
      while (true) {
        if (value < current.value) {
          // If the value is less than the current node, traverse left
          if (!current.left) {
            // If there is no left child, create a new node
            current.left = new TreeNode(value);
            break;
          }
          current = current.left;
        } else {
          // If the value is greater than or equal to the current node, traverse right
          if (!current.right) {
            // If there is no right child, create a new node
            current.right = new TreeNode(value);
            break;
          }
          current = current.right;
        }
      }
    }
  
    // Searches for a value in the tree
    search(value) {
      let current = this.root;
      while (current) {
        if (current.value === value) {
          // If the value is found, return true
          return true;
        } else if (value < current.value) {
          // If the value is less than the current node, traverse left
          current = current.left;
        } else {
          // If the value is greater than the current node, traverse right
          current = current.right;
        }
      }
      // If the value is not found, return false
      return false;
    }
  }  
/* Notes

*/

//Binary Search Trees
//Example implementation of a binary search tree:
// This class represents a node in a binary search tree
class Node {
    constructor(value) {
      // The value of the node
      this.value = value;
      // The left child of the node
      this.left = null;
      // The right child of the node
      this.right = null;
    }
  }
  
  // This class represents a binary search tree
  class BinarySearchTree {
    constructor() {
      // The root node of the tree
      this.root = null;
    }
  
    // Inserts a value into the tree
    insert(value) {
      // Create a new node with the given value
      let node = new Node(value);
  
      if (this.root === null) {
        // If the tree is empty, set the new node as the root
        this.root = node;
        return this;
      } else {
        // Traverse the tree to find the appropriate position for the new node
        let current = this.root;
        while (true) {
          if (value === current.value) {
            // If the value already exists in the tree, return undefined
            return undefined;
          } else if (value < current.value) {
            // If the value is less than the current node, traverse left
            if (current.left === null) {
              // If there is no left child, set the new node as the left child
              current.left = node;
              return this;
            } else {
              current = current.left;
            }
          } else {
            // If the value is greater than the current node, traverse right
            if (current.right === null) {
              // If there is no right child, set the new node as the right child
              current.right = node;
              return this;
            } else {
              current = current.right;
            }
          }
        }
      }
    }
  
    // Searches for a value in the tree
    find(value) {
      if (this.root === null) {
        // If the tree is empty, return false
        return false;
      }
      let current = this.root;
      let found = false;
      while (current && !found) {
        if (value < current.value) {
          // If the value is less than the current node, traverse left
          current = current.left;
        } else if (value > current.value) {
          // If the value is greater than the current node, traverse right
          current = current.right;
        } else {
          // If the value is equal to the current node, it has been found
          found = true;
        }
      }
      if (!found) {
        // If the value is not found, return false
        return false;
      }
      return current;
    }
  }
/* Notes

*/

//AVL Trees
/* Notes

*/
//Example Implementation of AVL Trees:
class AVLNode {
    constructor(value) {
      this.value = value;
      this.left = null;
      this.right = null;
      this.height = 1; // initialize height to 1
    }
  }
  
  class AVLTree {
    constructor() {
      this.root = null;
    }
  
    // Helper function to get height of a node
    getHeight(node) {
      if (!node) return 0;
      return node.height;
    }
  
    // Helper function to update height of a node
    updateHeight(node) {
      node.height = Math.max(this.getHeight(node.left), this.getHeight(node.right)) + 1;
    }
  
    // Helper function to get balance factor of a node
    getBalanceFactor(node) {
      if (!node) return 0;
      return this.getHeight(node.left) - this.getHeight(node.right);
    }
  
    // Helper function to perform left rotation
    leftRotate(node) {
      let newRoot = node.right;
      node.right = newRoot.left;
      newRoot.left = node;
      this.updateHeight(node);
      this.updateHeight(newRoot);
      return newRoot;
    }
  
    // Helper function to perform right rotation
    rightRotate(node) {
      let newRoot = node.left;
      node.left = newRoot.right;
      newRoot.right = node;
      this.updateHeight(node);
      this.updateHeight(newRoot);
      return newRoot;
    }
  
    // Helper function to balance a node
    balance(node) {
      this.updateHeight(node);
      let balanceFactor = this.getBalanceFactor(node);
      if (balanceFactor > 1) {
        if (this.getBalanceFactor(node.left) < 0) {
          node.left = this.leftRotate(node.left);
        }
        return this.rightRotate(node);
      }
      if (balanceFactor < -1) {
        if (this.getBalanceFactor(node.right) > 0) {
          node.right = this.rightRotate(node.right);
        }
        return this.leftRotate(node);
      }
      return node;
    }
  
    // Function to insert a value into the AVL tree
    insert(value) {
      let node = new AVLNode(value);
      if (!this.root) {
        this.root = node;
        return this;
      }
      this.insertNode(this.root, node);
      return this;
    }
  
    // Helper function to insert a node recursively
    insertNode(current, node) {
      if (!current) {
        return node;
      }
      if (node.value < current.value) {
        current.left = this.insertNode(current.left, node);
      } else {
        current.right = this.insertNode(current.right, node);
      }
      return this.balance(current);
    }
  }  
/* Notes

*/

//Red-Black Trees 
/* Notes

*/
//Example Implementation of Red-Black Trees:


class Node {
    constructor(value, color) {
      this.value = value;
      this.left = null;
      this.right = null;
      this.parent = null;
      this.color = color;
    }
  }
  
  class RedBlackTree {
    constructor() {
      this.root = null;
      this.colors = {
        RED: "red",
        BLACK: "black"
      };
    }
  
    // Helper function to get the grandparent of a node
    getGrandparent(node) {
      if (node && node.parent) {
        return node.parent.parent;
      }
      return null;
    }
  
    // Helper function to get the uncle of a node
    getUncle(node) {
      let grandparent = this.getGrandparent(node);
      if (grandparent) {
        if (node.parent === grandparent.left) {
          return grandparent.right;
        } else {
          return grandparent.left;
        }
      }
      return null;
    }
  
    // Helper function to perform left rotation
    leftRotate(node) {
      let newRoot = node.right;
      node.right = newRoot.left;
      if (newRoot.left !== null) {
        newRoot.left.parent = node;
      }
      newRoot.parent = node.parent;
      if (node.parent === null) {
        this.root = newRoot;
      } else if (node === node.parent.left) {
        node.parent.left = newRoot;
      } else {
        node.parent.right = newRoot;
      }
      newRoot.left = node;
      node.parent = newRoot;
    }
  
    // Helper function to perform right rotation
    rightRotate(node) {
      let newRoot = node.left;
      node.left = newRoot.right;
      if (newRoot.right !== null) {
        newRoot.right.parent = node;
      }
      newRoot.parent = node.parent;
      if (node.parent === null) {
        this.root = newRoot;
      } else if (node === node.parent.right) {
        node.parent.right = newRoot;
      } else {
        node.parent.left = newRoot;
      }
      newRoot.right = node;
      node.parent = newRoot;
    }
  
    // Helper function to insert a node
    insert(value) {
      let newNode = new Node(value, this.colors.RED);
      if (this.root === null) {
        this.root = newNode;
        newNode.color = this.colors.BLACK;
        return;
      }
      let current = this.root;
      while (current) {
        if (value < current.value) {
          if (current.left === null) {
            current.left = newNode;
            newNode.parent = current;
            break;
          }
          current = current.left;
        } else {
          if (current.right === null) {
            current.right = newNode;
            newNode.parent = current;
            break;
          }
          current = current.right;
        }
      }
      this.fixViolation(newNode);
    }
  
    // Helper function to fix any violations after inserting a node
    fixViolation(node) {
      while (node !== this.root && node.parent.color === this.colors.RED) {
        let grandparent = this.getGrandparent(node);
        if (node.parent === grandparent.left) {
          let uncle = this.getUncle(node);
          if (uncle && uncle.color === this.colors.RED) {
            node.parent.color = this.colors.BLACK;
            uncle.color = this.colors.BLACK;
            grandparent.color = this.colors.RED;
            node = grandparent;
          } else {
            if (node === node.parent.right) {
              node = node.parent;
              this.leftRotate(node);
            }
            node.parent.color = this.colors.BLACK;
  grandparent.color = this.colors.RED;
  this.rightRotate(grandparent);
  }
  } else {
  let uncle = this.getUncle(node);
  if (uncle && uncle.color === this.colors.RED) {
  node.parent.color = this.colors.BLACK;
  uncle.color = this.colors.BLACK;
  grandparent.color = this.colors.RED;
  node = grandparent;
  } else {
  if (node === node.parent.left) {
  node = node.parent;
  this.rightRotate(node);
  }
  node.parent.color = this.colors.BLACK;
  grandparent.color = this.colors.RED;
  this.leftRotate(grandparent);
  }
  }
  }
  this.root.color = this.colors.BLACK;
  }
  }
  
  module.exports = {
  Node,
  RedBlackTree
  };
//

//B-trees
/* Notes

*/
//Example Implementation of B-trees:


class BTreeNode {
    constructor(t, leaf = true) {
      this.keys = [];
      this.children = [];
      this.leaf = leaf;
      this.t = t;
    }
  
    // Traverse the subtree rooted with this node
    traverse() {
      let i;
      for (i = 0; i < this.keys.length; i++) {
        if (!this.leaf) {
          this.children[i].traverse();
        }
        console.log(this.keys[i]);
      }
  
      // Print the subtree rooted with last child
      if (!this.leaf) {
        this.children[i].traverse();
      }
    }
  
    // Search for a key in this node
    search(key) {
      let i = 0;
      while (i < this.keys.length && key > this.keys[i]) {
        i++;
      }
      if (this.keys[i] === key) {
        return this;
      }
      if (this.leaf) {
        return null;
      }
      return this.children[i].search(key);
    }
  
    // Insert a key into this node
    insertNonFull(key) {
      let i = this.keys.length - 1;
      if (this.leaf) {
        while (i >= 0 && this.keys[i] > key) {
          this.keys[i + 1] = this.keys[i];
          i--;
        }
        this.keys[i + 1] = key;
      } else {
        while (i >= 0 && this.keys[i] > key) {
          i--;
        }
        if (this.children[i + 1].keys.length === 2 * this.t - 1) {
          this.splitChild(i + 1, this.children[i + 1]);
          if (this.keys[i + 1] < key) {
            i++;
          }
        }
        this.children[i + 1].insertNonFull(key);
      }
    }
  
    // Split a child node of this node
    splitChild(i, y) {
      let z = new BTreeNode(y.t, y.leaf);
      z.keys = y.keys.splice(this.t);
      if (!y.leaf) {
        z.children = y.children.splice(this.t);
      }
      this.keys.splice(i, 0, y.keys[this.t - 1]);
      this.children.splice(i + 1, 0, z);
    }
  }
  
  class BTree {
    constructor(t) {
      this.root = null;
      this.t = t;
    }
  
    // Traverse the tree
    traverse() {
      if (this.root !== null) {
        this.root.traverse();
      }
    }
  
    // Search for a key in the tree
    search(key) {
      return this.root === null ? null : this.root.search(key);
    }
  
    // Insert a key into the tree
    insert(key) {
      if (this.root === null) {
        this.root = new BTreeNode(this.t, true);
        this.root.keys[0] = key;
        return;
      }
      if (this.root.keys.length === 2 * this.t - 1) {
        let s = new BTreeNode(this.t, false);
        s.children[0] = this.root;
        s.splitChild(0, this.root);
        let i = 0;
        if (s.keys[0] < key) {
          i++;
        }
        s.children[i].insertNonFull(key);
        this.root = s;
      } else {
        this.root.insertNonFull(key);
      }
    }
  }
//

//Trie Trees
/* Notes

*/
//Example Implementation of Trie Trees:
class TrieNode {
    constructor() {
      this.children = {};
      this.endOfWord = false;
    }
  }
  
  class Trie {
    constructor() {
      this.root = new TrieNode();
    }
  
    // Inserts a word into the trie
    insert(word) {
      let currentNode = this.root;
      for (let i = 0; i < word.length; i++) {
        let ch = word.charAt(i);
        if (!currentNode.children[ch]) {
          currentNode.children[ch] = new TrieNode();
        }
        currentNode = currentNode.children[ch];
      }
      currentNode.endOfWord = true;
    }
  
    // Searches for a word in the trie
    search(word) {
      let currentNode = this.root;
      for (let i = 0; i < word.length; i++) {
        let ch = word.charAt(i);
        if (!currentNode.children[ch]) {
          return false;
        }
        currentNode = currentNode.children[ch];
      }
      return currentNode.endOfWord;
    }
  
    // Returns true if the trie has any word that starts with the given prefix
    startsWith(prefix) {
      let currentNode = this.root;
      for (let i = 0; i < prefix.length; i++) {
        let ch = prefix.charAt(i);
        if (!currentNode.children[ch]) {
          return false;
        }
        currentNode = currentNode.children[ch];
      }
      return true;
    }
  }
  
  // Example usage:
  let trie = new Trie();
  trie.insert("hello");
  trie.insert("world");
  console.log(trie.search("hello")); // true
  console.log(trie.search("world")); // true
  console.log(trie.search("foo")); // false
  console.log(trie.startsWith("he")); // true
  console.log(trie.startsWith("foo")); // false
//  

//Minimum Spanning Tree in Graph 
/* Notes

*/
//Example of Kruskal's Algorithm that finds the minimum spanning tree of a graph represented as an adjacency list:
// Implementation of the Union-Find data structure to support finding and merging sets
class UnionFind {
    constructor(n) {
      // Initialize the parent array so that each element is in its own set
      this.parent = new Array(n);
      for (let i = 0; i < n; i++) {
        this.parent[i] = i;
      }
    }
  
    // Find the root of the set that x belongs to
    find(x) {
      // Keep traversing up the parent array until we reach the root of the set
      while (x !== this.parent[x]) {
        // Path compression: set the parent of x to the root of its set to flatten the tree
        this.parent[x] = this.parent[this.parent[x]];
        x = this.parent[x];
      }
      // Return the root of the set
      return x;
    }
  
    // Merge the sets that x and y belong to
    union(x, y) {
      let rootX = this.find(x);
      let rootY = this.find(y);
      // If x and y are in different sets, set the parent of the root of x's set to the root of y's set
      if (rootX !== rootY) {
        this.parent[rootX] = rootY;
      }
    }
  }
  
  // Implementation of Kruskal's algorithm to find the minimum spanning tree of an undirected, weighted graph
  function kruskal(graph) {
    // Create a list of all edges in the graph, represented as (u, v, weight) tuples
    let edges = [];
    for (let node in graph) {
      for (let neighbor in graph[node]) {
        edges.push([node, neighbor, graph[node][neighbor]]);
      }
    }
  
    // Sort the edges in ascending order of weight
    edges.sort((a, b) => a[2] - b[2]);
  
    // Initialize the Union-Find data structure with the number of nodes in the graph
    let uf = new UnionFind(Object.keys(graph).length);
  
    // Initialize the minimum spanning tree and its weight
    let mst = [];
    let weight = 0;
  
    // Iterate over the edges in ascending order of weight
    for (let [u, v, w] of edges) {
      // If u and v are not already connected, add the edge to the minimum spanning tree
      if (uf.find(u) !== uf.find(v)) {
        uf.union(u, v);
        mst.push([u, v]);
        weight += w;
      }
    }
  
    // Return the minimum spanning tree and its weight as an object
    return { mst, weight };
  }
  
  // Example usage
  let graph7 = {
    'A': { 'B': 3, 'C': 1 },
    'B': { 'A': 3, 'C': 4, 'D': 2 },
    'C': { 'A': 1, 'B': 4, 'D': 5 },
    'D': { 'B': 2, 'C': 5 },
  };
  
  let result = kruskal(graph7);
  console.log('Minimum spanning tree:', result.mst);
  console.log('Total weight:', result.weight);  
/* Notes

*/

//Example of Primm's Algorithm that finds the minimum spanning tree of a graph represented as an adjacency list:
// Define the graph as an adjacency matrix with edge weights
const graph8 = [
    [0, 2, 3, 0],
    [2, 0, 1, 5],
    [3, 1, 0, 4],
    [0, 5, 4, 0],
  ];
  
  // Find the minimum spanning tree using Prim's algorithm
  function primMST(graph) {
    // Initialize an empty set to store visited nodes
    const visited = new Set();
  
    // Initialize a priority queue to store edges, ordered by their weight
    const queue = new PriorityQueue((a, b) => a.weight - b.weight);
  
    // Initialize a map to store the parent node for each node in the MST
    const parent = new Map();
  
    // Start with an arbitrary node (node 0) and add its edges to the queue
    for (let j = 0; j < graph.length; j++) {
      if (graph[0][j] !== 0) {
        queue.enqueue({ u: 0, v: j, weight: graph[0][j] });
      }
    }
  
    // While there are still nodes to visit and edges in the queue
    while (visited.size < graph.length && !queue.isEmpty()) {
      // Dequeue the edge with the smallest weight
      const edge = queue.dequeue();
  
      // If both nodes of the edge have been visited, skip it
      if (visited.has(edge.u) && visited.has(edge.v)) {
        continue;
      }
  
      // Add the edge to the MST and mark both nodes as visited
      const u = edge.u;
      const v = edge.v;
      const weight = edge.weight;
      parent.set(v, u);
      visited.add(u);
      visited.add(v);
  
      // Add the edges of the newly visited node to the queue
      for (let j = 0; j < graph.length; j++) {
        if (graph[v][j] !== 0 && !visited.has(j)) {
          queue.enqueue({ u: v, v: j, weight: graph[v][j] });
        }
      }
    }
  
    return parent;
  }
  
  // Find the minimum spanning tree of the graph
  const parent = primMST(graph8);
  
  // Print the edges of the minimum spanning tree
  console.log("Edges of the minimum spanning tree:");
  for (let i = 1; i < graph8.length; i++) {
    console.log(`${parent.get(i)} - ${i}`);
  }
/* Notes

*/