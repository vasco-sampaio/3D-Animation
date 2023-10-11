#pragma once

#include "cgp/cgp.hpp"
// #include "octree_node.hpp"

struct particle_element;

struct OctantFunctor {
    float _octantSize;
    unsigned int _octantNbPerDim;
    cgp::vec3 _min;

    OctantFunctor(unsigned int depth, cgp::vec3 min, cgp::vec3 max) : _min(min) {
        this->_octantNbPerDim = 1 << depth;
        this->_octantSize = (max.x - min.x) / _octantNbPerDim;
    }

    unsigned int operator() (const cgp::vec3& position) const {
        unsigned int x = static_cast<unsigned int>((position.x - _min.x) / _octantSize);
        unsigned int y = static_cast<unsigned int>((position.y - _min.y) / _octantSize);
        unsigned int z = static_cast<unsigned int>((position.z - _min.z) / _octantSize);

        // assuming the space is a cube
        return x + y * _octantNbPerDim + z * _octantNbPerDim * _octantNbPerDim;
    }
};

template <typename Container>
concept ContainerConcept = requires(Container c) {
    typename Container::value_type;
    typename Container::iterator;
};

/*
** The Morton code length is related to the number of levels in your octree.
** Specifically, the Morton code length should be three times the number of levels in
** your octree since each level corresponds to 3 bits in the Morton code (due to 8 children per node).
*/
template <typename T>
class Octree {
private:
    OctantFunctor _functor;

    unsigned int _depth;
    mutable std::vector<T> _nodes;

    /* builds only the leaf nodes */
    void buildTree(unsigned int depth) {
        if (depth == 0) {
            _nodes.push_back(T{});
            return;
        }

        for (int i = 0; i < 8; ++i)
            buildTree( depth - 1);
    }

public:
    explicit Octree(unsigned int depth, cgp::vec3 min = { -1.0f, -1.0f, -1.0f }, cgp::vec3 max = { 1.0f, 1.0f, 1.0f }) :
            _functor(depth, min, max), _depth(depth), _nodes() {
        if (_depth > 0) {
            _nodes.reserve(std::pow(8, depth));
            buildTree(depth);
        }
    }

    unsigned int get_depth() const { return _depth; }

    unsigned int get_octant(const cgp::vec3& position) const {
        return _functor(position);
    }

    const T& get_node(unsigned int octant) const {
        return _nodes[octant];
    }

    // Enable these methods only if the Container satisfies the ContainerConcept
    template <typename Container = T>
    std::enable_if_t<ContainerConcept<Container>, void> insert_into_node(unsigned int octant, const typename Container::value_type& value) const {
        (void)_nodes[octant].insert(value);
    }

    template <typename Container = T>
    std::enable_if_t<ContainerConcept<Container>, void> erase_from_node(unsigned int octant, const typename Container::value_type& value) const {
        (void)_nodes[octant].erase(value);
    }
};

void update_grid(cgp::numarray<particle_element>& particles, const Octree<std::set<int>>& octree);
