#pragma once

#include "cgp/cgp.hpp"
#include "octree_node.hpp"

struct particle_element;

struct OctreeData {
    const std::unordered_map<unsigned int, unsigned int> octantMortonMap;
    const float octantSizeX;
    const float octantSizeY;
    const float octantSizeZ;

    /* Space dimensions represented by the octree */
    const cgp::vec3 min;
    const cgp::vec3 max;

    OctreeData(std::unordered_map<unsigned int, unsigned int> octantMortonMap, float octantSizeX, float octantSizeY, float octantSizeZ, cgp::vec3 min, cgp::vec3 max) :
        octantMortonMap(std::move(octantMortonMap)), octantSizeX(octantSizeX), octantSizeY(octantSizeY), octantSizeZ(octantSizeZ), min(min), max(max)
    {}

    /* Given a position returns the matching morton code */
    unsigned int mortonCode(const cgp::vec3& position) const {
        static unsigned int height = static_cast<unsigned int>(max.x) - static_cast<unsigned int>(min.x);
        static unsigned int width = static_cast<unsigned int>(max.y) - static_cast<unsigned int>(min.y);

        static unsigned int nbOctantsX = static_cast<unsigned int>(height) / octantSizeX;
        static unsigned int nbOctantsY = static_cast<unsigned int>(width) / octantSizeY;

        unsigned int x = static_cast<unsigned int>((position.x - min.x) / octantSizeX);
        unsigned int y = static_cast<unsigned int>((position.y - min.y) / octantSizeY);
        unsigned int z = static_cast<unsigned int>((position.z - min.z) / octantSizeZ);
        return octantMortonMap.at(x + y * nbOctantsY + z * nbOctantsX * nbOctantsY);
    }
};

/*
** The Morton code length is related to the number of levels in your octree.
** Specifically, the Morton code length should be three times the number of levels in
** your octree since each level corresponds to 3 bits in the Morton code (due to 8 children per node).
*/
template <typename T>
class Octree {
private:
    static std::unordered_map<unsigned int, unsigned int> buildOctantMortonMap(unsigned int depth, cgp::vec3 max, cgp::vec3 min) {
        std::unordered_map<unsigned int, unsigned int> octantMortonMap;

        float xDim = max.x - min.x;
        float yDim = max.y - min.y;
        float zDim = max.z - min.z;

        float octantSizeX = xDim / static_cast<float>(1 << depth);
        float octantSizeY = yDim / static_cast<float>(1 << depth);
        float octantSizeZ = zDim / static_cast<float>(1 << depth);

        unsigned int nbOctantsX = static_cast<unsigned int>(xDim) / octantSizeX;
        unsigned int nbOctantsY = static_cast<unsigned int>(yDim) / octantSizeY;
        unsigned int nbOctantsZ = static_cast<unsigned int>(zDim) / octantSizeZ;

        auto morton_code = [](unsigned int x, unsigned int y, unsigned int z) {
            unsigned int offset = 0;
            for (unsigned int i = 0; i < 10; ++i) {
                offset |= ((x & (1 << i)) << 2 * i) | ((y & (1 << i)) << (2*i + 1)) | ((z & (1 << i)) << (2*i + 2));
            }
            return offset;
        };

        for (unsigned int x = 0; x < nbOctantsX; ++x) {
            for (unsigned int y = 0; y < nbOctantsY; ++y) {
                for (unsigned int z = 0; z < nbOctantsZ; ++z) {
                    unsigned int key = x + nbOctantsX * (y + nbOctantsY * z);
                    octantMortonMap[key] = morton_code(x, y, z);
                }
            }
        }

        return octantMortonMap;
    }

    OctreeNode<T>* _root;
    unsigned int _depth;
    std::unordered_map<unsigned int, OctreeNode<T>*> _nodes;

    void buildTree(OctreeNode<T>* node, unsigned int mortonCode, unsigned int depth) {
        if (depth == 0) {
            _nodes.emplace(mortonCode, node);
            return;
        }

        for (int i = 0; i < 8; ++i) {
            unsigned int childMortonCode = (mortonCode << 3) | i;
            node->_children[i] = new OctreeNode<T>(); // Initialize value with a default-constructed T
            node->_children[i]->_parent = node;
            buildTree(node->_children[i], childMortonCode, depth - 1);
        }
    }

public:
    OctreeData data;

    explicit Octree(unsigned int depth, cgp::vec3 max = { 1.0f, 1.0f, 1.0f },
                    cgp::vec3 min = { -1.0f, -1.0f, -1.0f }) : _root(nullptr), _depth(depth), _nodes(), data(buildOctantMortonMap(depth, max, min),
                                                                                                   (max.x - min.x) / static_cast<float>(1 << depth),
                                                                                                   (max.y - min.y) / static_cast<float>(1 << depth),
                                                                                                   (max.z - min.z) / static_cast<float>(1 << depth),
                                                                                                   min, max) {
        if (_depth > 0) {
            _root = new OctreeNode<T>(); // Initialize value with a default-constructed T
            buildTree(_root, 0, depth);
        }
    }

    ~Octree() {
        delete _root;
    }

    [[nodiscard]] unsigned int get_depth() const { return _depth; }
    const OctreeNode<T>* get_root() const { return _root; }

    const OctreeNode<T>* get_node(unsigned int mortonCode) const {
        return _nodes.at(mortonCode);
    }
};

void update_grid(cgp::numarray<particle_element>& particles, const Octree<std::set<unsigned int>>& octree);
