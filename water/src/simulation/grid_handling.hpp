#pragma once

#include "cgp/cgp.hpp"
#include "simulation.hpp"

uint morton_code(uint x, uint y, uint z);
void update_grid(cgp::numarray<particle_element>& particles, std::vector<uint>& grid);
