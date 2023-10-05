#pragma once

#include "cgp/cgp.hpp"
#include "octree_node.hpp"

struct particle_element;

/*
** The Morton code length is related to the number of levels in your octree.
** Specifically, the Morton code length should be three times the number of levels in
** your octree since each level corresponds to 3 bits in the Morton code (due to 8 children per node).
*/
template <typename T>
class Octree {
private:
    OctreeNode<T>* _root;
    unsigned int _depth;

    void buildTree(OctreeNode<T>* node, unsigned int mortonCode, unsigned int depth) {
        if (depth == 0) {
            return;
        }

        for (int i = 0; i < 8; ++i) {
            unsigned int childMortonCode = (mortonCode << 3) | i;
            node->_children[i] = new OctreeNode<T>(childMortonCode, T{}); // Initialize value with a default-constructed T
            node->_children[i]->_parent = node;
            buildTree(node->_children[i], childMortonCode, depth - 1);
        }
    }

public:
    Octree(unsigned int depth) : _root(nullptr), _depth(depth) {
        if (_depth > 0) {
            _root = new OctreeNode<T>(0U); // Initialize value with a default-constructed T
            buildTree(_root, 0, depth);
        }
    }

    ~Octree() {
        delete _root;
    }

    unsigned int get_depth() const { return _depth; }
    const OctreeNode<T>* get_root() const { return _root; }

    const OctreeNode<T>* searchValue(unsigned int mortonCode) const {
        return _root->searchValue(mortonCode);
    }
};

void update_grid(cgp::numarray<particle_element>& particles, const Octree<std::set<int>>& octree);
