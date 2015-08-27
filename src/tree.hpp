#ifndef  __TREE_HPP__
#define __TREE_HPP__
//#include <list>
#include <deque>
#include <vector>
#include "elt.hpp"
#include "node.hpp"

namespace sphereRemap {

using namespace std;

class CBasicTree
{
public:

	NodePtr root; /* The main tree is stored as Nodes which can be reached through traversal starting here */
	NodePtr ref; // FIXME this reference, set by a node is odd, try to remove
	int ri; /** this is set to one by a node in case of reinsertion */
	vector<int> levelSize; /** e.g. levelSize[0] == leafs.size() */
	vector<Node> leafs; /** leafs are stored in vector for easy access and rest of the tree nodes as separate allocations, only reachable through tree traversal */

	CBasicTree(); 
	~CBasicTree(); 
	void build(vector<Node>& nodes);
	void slim(int nbIts = 1);
	virtual void insertNodes(vector<Node>& node) = 0;

	void routeNodes(vector<int>& route, vector<Node>& nodes, int assignLevel);
	void routeIntersections(vector<vector<int> >& route, vector<Node>& nodes);

	void push_back(NodePtr node);
	void push_front(NodePtr node);
	void increaseLevelSize(int level);
	void decreaseLevelSize(int level);
	void newRoot(int level);
	void insertNode(NodePtr node);
  void output(ostream& flux, int level) ;

private:
	deque<NodePtr > pool;
	void emptyPool();

};

class CTree : public CBasicTree
{
public:
	void insertNodes(vector<Node>& nodes);
};

class CSampleTree : public CBasicTree
{
	int keepNodes;
	int assignLevel;
public:
	CSampleTree(int keepNodes, int assignLevel) : keepNodes(keepNodes), assignLevel(assignLevel) {} 

	void insertNodes(vector<Node>& nodes);
};

}
#endif
