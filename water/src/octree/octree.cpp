#include "octree.hpp"
#include "../simulation/simulation.hpp"

using namespace cgp;


sph_parameters_structure sph_parameters;

using s_int = std::set<int>;

// 1D offset corresponding to the index (k1, k2, k3) in a 3D grid according to a Morton spatial sorting
unsigned int morton_code(unsigned int x, unsigned int y, unsigned int z) {
    unsigned int offset = 0;
    for (unsigned int i = 0; i < 10; ++i) {
        offset |= ((x & (1 << i)) << 2*i) | ((y & (1 << i)) << (2*i + 1)) | ((z & (1 << i)) << (2*i + 2));
    }
    return offset;
}

unsigned int octant(vec3 pos, unsigned int depth, int3 max = { 1, 1, 1 }, int3 min = { -1, -1, -1 }) {
    float octantSizeX = (max.x - min.x) / (1 << depth);
    float octantSizeY = (min.y - min.y) / (1 << depth);
    float octantSizeZ = (min.z - min.z) / (1 << depth);


    int octantX = (pos.x - min.x) / octantSizeX;
    int octantY = (pos.y - min.y) / octantSizeY;
    int octantZ = (pos.z - min.z) / octantSizeZ;

    return morton_code(octantX, octantY, octantZ);

}

// Assuming the node is a leaf node and there are several particles in each
std::vector<int> get_neighbors(int id, const OctreeNode<s_int>* node) {
    std::vector<int> neighbors;

    // take advantage of the morton spacial sorting
    s_int c_set = node->get_siblings_values();
    for (const auto& p : c_set) {
        if (p != id) { // TODO: calculate norm and remove norm calculation from the simulation
            neighbors.push_back(p);
        }
    }

    return neighbors;
}


void update_grid(numarray<particle_element>& particles, const Octree<s_int>& grid)
{
    for (int i = 0; i < particles.size(); ++i)
    {
        particle_element& prt = particles[i];

        unsigned int code = octant(prt.p, grid.get_depth());

        if (code == prt.morton) {
            continue;
        }

        const OctreeNode<s_int>* node = grid.searchValue(prt.morton);
        node->get_value().erase(prt.morton);

        prt.morton = code;
        const OctreeNode<s_int>* new_node = grid.searchValue(code);
        
        // All the new neighbors won't be totally computed
        prt.neighbors.clear();
        prt.neighbors.reserve(30);
        prt.neighbors = get_neighbors(i, new_node);

    }
}


