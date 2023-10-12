#pragma once

#include "../octree/octree.hpp"
#include "cgp/cgp.hpp"

// SPH Particle
struct ParticleArray {
    std::vector<cgp::vec3> positions;
    std::vector<cgp::vec3> speeds;
    std::vector<cgp::vec3> forces;
    std::vector<float> rho;
    std::vector<float> pressure;
    std::vector<unsigned int> octants;

    ParticleArray(size_t numParticles) :
            positions(numParticles),
            speeds(numParticles),
            forces(numParticles),
            rho(numParticles),
            pressure(numParticles),
            octants(numParticles) {
        // You can add additional initialization here if needed.
    }

    [[nodiscard]] unsigned int size() const {
        return positions.size();
    }

    void push_back(cgp::vec3 const& position, cgp::vec3 const& speed, cgp::vec3 const& force, float rho, float pressure, unsigned int octant) {
        positions.push_back(position);
        speeds.push_back(speed);
        forces.push_back(force);
        this->rho.push_back(rho);
        this->pressure.push_back(pressure);
        octants.push_back(octant);
    }
};

// SPH simulation parameters
struct sph_parameters_structure
{
    // Initial particle spacing (relative to h)
    float const c = 0.7f;

    // Influence distance of a particle (size of the kernel)
    float h = 0.2f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

     // Total mass of a particle (consider rho0 h^2)
    float m = rho0*h*h;

    // viscosity parameter
    float nu = 0.02f;
     
    // Stiffness converting density to pressure
    float stiffness = 8.0f;
    
};


void simulate(float dt, ParticleArray& particles, sph_parameters_structure const& sph_parameters,
              const Octree<std::set<int>>& grid);