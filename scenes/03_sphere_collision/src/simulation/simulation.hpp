#pragma once

#include "cgp/cgp.hpp"


struct particle_structure
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed

    cgp::vec3 c; // Color
    float r;     // Radius
    float m;     // mass
};

void simulate(std::vector<particle_structure>& particles, float dt);
void manage_cube_collisions(std::vector<particle_structure>& particles);
void manage_sphere_collisions(std::vector<particle_structure>& particles);

