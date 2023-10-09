#include "octree.hpp"
#include "../simulation/simulation.hpp"

using namespace cgp;


sph_parameters_structure sph_parameters;

using s_uint = std::set<unsigned int>;

void update_grid(numarray<particle_element>& particles, const Octree<s_uint>& grid)
{
    for (int id = 0; id < particles.size(); ++id)
    {
        particle_element& prt = particles[id];
        unsigned int code = grid.data.mortonCode(prt.p);

        if (code == prt.morton) {
            continue;
        }

        if (code != 0) {
            const OctreeNode<s_uint> *node = grid.get_node(prt.morton);
            if (node != nullptr)
                node->get_value().erase(id);
        }

        prt.morton = code;
        const OctreeNode<s_uint>* new_node = grid.get_node(code);
        if (new_node != nullptr)
            new_node->get_value().insert(id);
    }
}


