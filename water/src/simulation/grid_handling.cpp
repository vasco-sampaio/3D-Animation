#include <algorithm>

#include "grid_handling.hpp"


using namespace cgp;

sph_parameters_structure sph_parameters;

// 1D offset corresponding to the index (k1, k2, k3) in a 3D grid according to a Morton spatial sorting
uint morton_code(uint x, uint y, uint z) {
    uint offset = 0;
    for (uint i = 0; i < 10; ++i) {
        offset |= ((x & (1 << i)) << 2*i) | ((y & (1 << i)) << (2*i + 1)) | ((z & (1 << i)) << (2*i + 2));
    }
    return offset;
}


uint binary_search_lower_bound(const std::vector<uint>& data, uint morton)
{
    auto it = std::lower_bound(data.begin(), data.end(), morton,
        [](uint m1, uint m2) { return m1 < m2; });
    if (it == data.end() || *it != morton)
        return data.size();
    return std::distance(data.begin(), it);
}


uint binary_search_upper_bound(const std::vector<uint>& data, uint morton)
{
    auto it = std::upper_bound(data.begin(), data.end(), morton,
        [](uint m1, uint m2) { return m1 < m2; });
    if (it == data.end() || *it != morton)
        return data.size();
    return std::distance(data.begin(), it);
}


std::pair<uint, uint> binary_search_bounds(
    const std::vector<uint>& data, const vec3& target) 
{
    float h = sph_parameters.h;

    int3 lowerBound;
    lowerBound.x = target.x - h;
    lowerBound.y = target.y - h;
    lowerBound.z = target.z - h;
    uint lowerMorton = morton_code(floor(lowerBound.x * 100), floor(lowerBound.y * 100), floor(lowerBound.z * 100));

    int3 upperBound;
    upperBound.x = target.x + h;
    upperBound.y = target.y + h;
    upperBound.z = target.z + h;
    uint upperMorton = morton_code(floor(upperBound.x * 100), floor(upperBound.y * 100), floor(upperBound.z * 100));

    auto lower = binary_search_lower_bound(data, lowerMorton);
    auto upper = binary_search_upper_bound(data, upperMorton);

    return std::make_pair(lower, upper);
}


void get_neighbors(numarray<particle_element>& particles, const std::vector<uint>& grid)
{
    float h = sph_parameters.h;

    for (particle_element& prt : particles)
    {
        prt.neighbors.clear();
        
        std::pair<uint, uint> bounds = binary_search_bounds(grid, prt.p);

        for (uint i = bounds.first; i <= bounds.second; ++i)
        {
            if (particles[i].morton != prt.morton && norm(prt.p - particles[i].p) <= h) {
                prt.neighbors.push_back(i);
            }
        }
    }
}


void update_grid(numarray<particle_element>& particles, std::vector<uint>& grid)
{
    float h = sph_parameters.h;
    // update grid
    for (particle_element& prt : particles)
    {
        uint new_code = morton_code(floor((prt.p.x + 1)/h), floor((prt.p.y + 1)/h), floor((prt.p.z + 1)/h));
        if (prt.morton != new_code)
        {
            grid.erase(std::find(grid.begin(), grid.end(), prt.morton));
            prt.morton = new_code;
            grid.push_back(new_code);
        }
    }

    std::sort(grid.begin(), grid.end());
}