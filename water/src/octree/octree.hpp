#pragma once

#include "cgp/cgp.hpp"
// #include "octree_node.hpp"

struct particle_element;

struct OctreeData {
    const std::unordered_map<unsigned int, unsigned int> octantMortonMap;
    const float octantSize;
    const float octantNb;

    const cgp::vec3 min;
    const cgp::vec3 max;

    OctreeData(std::unordered_map<unsigned int, unsigned int> octantMortonMap, float octantSize, float octantNb, cgp::vec3 min, cgp::vec3 max) :
        octantMortonMap(std::move(octantMortonMap)), octantSize(octantSize), octantNb(octantNb), min(min), max(max) {}

    /* Given a position returns the matching morton code */
    unsigned int mortonCode(const cgp::vec3& position) const {
        unsigned int x = static_cast<unsigned int>((position.x - min.x) / octantSize);
        unsigned int y = static_cast<unsigned int>((position.y - min.y) / octantSize);
        unsigned int z = static_cast<unsigned int>((position.z - min.z) / octantSize);

        // assuming the space is a cube
        return octantMortonMap.at(x + y * octantNb + z * octantNb * octantNb);
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
    // assuming we are dealing with a cube
    static std::unordered_map<unsigned int, unsigned int> buildOctantMortonMap(float octantNb) {
        std::unordered_map<unsigned int, unsigned int> octantMortonMap;

        auto morton_code = [](unsigned int x, unsigned int y, unsigned int z) {
            unsigned int offset = 0;

            for (unsigned int i = 0; i < 10; ++i) {
                offset |= ((x & (1 << i)) << 2 * i) | ((y & (1 << i)) << (2*i + 1)) | ((z & (1 << i)) << (2*i + 2));
            }
            return offset;
        };

        for (unsigned int x = 0; x < octantNb; ++x) {
            for (unsigned int y = 0; y < octantNb; ++y) {
                for (unsigned int z = 0; z < octantNb; ++z) {
                    unsigned int key = x + y * octantNb + octantNb * octantNb * z;
                    octantMortonMap[key] = morton_code(x, y, z);
                }
            }
        }

        return octantMortonMap;
    }

    unsigned int _depth;
    mutable std::unordered_map<unsigned int, T> _nodes;

    /* builds only the leaf nodes */
    void buildTree(unsigned int mortonCode, unsigned int depth) {
        if (depth == 0) {
            _nodes.emplace(mortonCode, T{});
            return;
        }

        for (int i = 0; i < 8; ++i) {
            unsigned int childMortonCode = (mortonCode << 3) | i;
            buildTree(childMortonCode, depth - 1);
        }
    }


public:
    OctreeData data;


    explicit Octree(unsigned int depth, float dimension = 1.0f) : _depth(depth), _nodes(), data(buildOctantMortonMap(dimension/ static_cast<float>(1 << depth)),
                                                                                                   dimension / static_cast<float>(1 << depth),
                                                                                             dimension / static_cast<float>(1 << depth)) {
        if (_depth > 0)
            buildTree(0, depth);
    }


    [[nodiscard]] unsigned int get_depth() const { return _depth; }


    const T& get_node(unsigned int mortonCode) const {
        return _nodes.at(mortonCode);
    }

    // Enable these methods only if the Container satisfies the ContainerConcept
    template <typename Container = T>
    std::enable_if_t<ContainerConcept<Container>, void> insert_into_node(unsigned int mortonCode, const typename Container::value_type& value) const {
        (void)_nodes.at(mortonCode).insert(value);
    }

    template <typename Container = T>
    std::enable_if_t<ContainerConcept<Container>, void> erase_from_node(unsigned int mortonCode, const typename Container::value_type& value) const {
        (void)_nodes.at(mortonCode).erase(value);
    }
};

void update_grid(cgp::numarray<particle_element>& particles, const Octree<std::set<int>>& octree);
