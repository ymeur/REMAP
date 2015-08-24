#include <cassert>
#include "tree.hpp"
#include "node.h"
#include "elt.h"
#include <iostream>
#include <mpi.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <vector>
#include <set>

#include "timer.hpp"
#include "mpi_routing.hpp"

static const int MAX_LEVEL_SIZE = 100;

double cputime();
using namespace std;

CBasicTree::CBasicTree() : ri(0), levelSize(MAX_LEVEL_SIZE), root(NULL)
{
}

void CBasicTree::routeNodes(vector<int>& route, vector<Node>& nodes, int assignLevel)
{
	for (int i = 0; i < nodes.size(); i++)
	{
		root->routeNode(&nodes[i], assignLevel);
		route[i] = nodes[i].route;
	}
}

void CBasicTree::routeIntersections(vector<vector<int> >& routes, vector<Node>& nodes)
{
	for (int i = 0; i < nodes.size(); i++)
		root->routeIntersection(routes[i], &nodes[i]);
}

void CBasicTree::build(vector<Node>& nodes)
{
	newRoot(1);
	insertNodes(nodes);
}

void CBasicTree::slim(int nbIts)
{
	//if (!isSampleTree)
	for (int i = 0; i < nbIts; i++)
	{
		for (int level = root->level - 1; level > 0; level--)
		{
			slim2(root, level);
			ri = 0;
			emptyPool();
		}

		for (int level = 2; level < root->level; level++)
		{
			slim2(root, level);
			ri = 0;
			emptyPool();
		}
	}
}

void CBasicTree::insertNode(NodePtr node)
{
	node->tree = this;
	increaseLevelSize(0);
	push_back(node);

	NodePtr q;
	while (pool.size())
	{
		q = pool.front();
		pool.pop_front();
		q = insert(q, root);
		if (ri)
		{
			delete q;
			ri = 0;
		}
	}
}

void CBasicTree::emptyPool(void)
{
	while (pool.size())
	{
		NodePtr q = pool.front();
		pool.pop_front();
		q = insert(q, root);
		if (ri)
		{
			delete q;
			ri = 0;
		}
	}
}

void CBasicTree::push_back(NodePtr node)
{
#ifdef DEBUG
	node->parent = NULL; // will be overridden later, but for mem leak check better trigger leak as soon as possible
#endif
	pool.push_back(node);
}

void CBasicTree::push_front(NodePtr node)
{
#ifdef DEBUG
	node->parent = NULL; // will be overridden later, but for mem leak check better trigger leak as soon as possible
#endif
	pool.push_front(node);
}

void CBasicTree::increaseLevelSize(int level)
{
	levelSize[level]++;
}

void CBasicTree::decreaseLevelSize(int level)
{
	levelSize[level]--;
}

void CBasicTree::newRoot(int level)  // newroot <- root
{
	root = new Node;
	if (level > 1) cout << " " << level << " ";
	root->level = level;
	root->parent = 0;
	root->leafCount = 0;
	root->centre = ORIGIN;
	root->radius = 0.;
	root->reinserted = false;
	root->updateCount = 0;
	root->tree = this;
	levelSize[level]++;
}

CBasicTree::~CBasicTree()
{
	//FIXME uncomment the next line and figure out why it segfault sometimes
	//root->free_descendants(); // recursively deletes all nodes in the tree
	if (root) delete root;
}

void CTree::insertNodes(vector<Node>& nodes)
{
	int stepSlim = MAX_NODE_SZ*MAX_NODE_SZ*2;
	const double ratio = 1.5;
	for (int i = 0; i < nodes.size(); i++)
	{
		insertNode(&nodes[i]);

		if (root->leafCount > stepSlim) // time for slim down
		{
			slim();
			stepSlim = stepSlim * ratio;
		}
	}
}


void CSampleTree::insertNodes(vector<Node>& nodes)
{
	bool first = true;
	int stepSlim = MAX_NODE_SZ * MAX_NODE_SZ * keepNodes / 2;
	double ratio = 1 + (0.5 * (1 - (keepNodes - 1) / (1.0 * keepNodes)));
	for (int i = 0; i < nodes.size(); i++)
	{
		insertNode(&nodes[i]);

		if (root->leafCount > stepSlim) // time for slim down
		{
			slim();
			stepSlim = stepSlim * ratio;
		}
		if (levelSize[assignLevel] >= keepNodes && first)
		{
			slim(3);
			first = false;
		}

		if (levelSize[assignLevel] >= keepNodes)
		{
			cout << "SampleTree build : nb Node inserted : " << i << endl;
			break;
		}
	}
}
