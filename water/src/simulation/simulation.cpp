#include <algorithm>

#include "simulation.hpp"

using namespace cgp;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    float d = norm(p_i - p_j);
    assert_cgp_no_msg(d <= h);
    return (45/(3.14159f * std::pow(h, 6))) * (h - d);
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    float d = norm(p_i - p_j);
    assert_cgp_no_msg(d <= h);
    return (-45/(3.14159f * std::pow(h, 6))) * std::pow((h - d), 2) * ((p_i-p_j) / d);
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
	return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
}

void update_density_pressure_force(ParticleArray& particles, sph_parameters_structure const& sph_parameters, const Octree<std::set<int>>& grid)
{
    float const h = sph_parameters.h;
    float const m = sph_parameters.m;
    float const rho0 = sph_parameters.rho0;
    float const nu = sph_parameters.nu;
    float const stiffness = sph_parameters.stiffness;

    for (int i = 0; i < particles.size(); ++i)
    {
        // Update density
        float rho = 0;

        for (int j : grid.get_node(particles.octants[i])) {
            if (norm(particles.positions[i] - particles.positions[j]) <= h)
                rho += m * W_density(particles.positions[i], particles.positions[j], h);
        }

        particles.rho[i] = rho;

        // Compute pressure
        particles.pressure[i] = density_to_pressure(rho, rho0, stiffness);

        // Compute forces
        vec3 F_pressure;
        vec3 F_viscosity;
        vec3 sum_press;
        vec3 sum_visc;
        sum_press = vec3{0,0,0};
        sum_visc = vec3{0,0,0};

        for (int j : grid.get_node(particles.octants[i]))
        {
            if (i != j && norm(particles.positions[i] - particles.positions[j]) <= h) {
                sum_press += m * ((particles.pressure[i] + particles.pressure[j]) / (2.0f * particles.rho[j])) *
                             W_gradient_pressure(particles.positions[i], particles.positions[j], h);

                sum_visc += m * ((particles.speeds[j] - particles.speeds[i]) / particles.rho[j]) *
                            W_laplacian_viscosity(particles.positions[i], particles.positions[j], h);
            }
        }

        F_pressure = float(-m / particles.rho[i]) * sum_press;
        F_viscosity = m * nu * sum_visc;
        particles.forces[i] = m * vec3{0,-9.81f,0} + F_pressure + F_viscosity;
    }
}

void simulate(float dt, ParticleArray& particles, sph_parameters_structure const& sph_parameters, const Octree<std::set<int>>& grid)
{
	// Update values
    update_density_pressure_force(particles, sph_parameters, grid);

    // Numerical integration
	float const damping = 0.005f;
	float const m = sph_parameters.m;

    // Collision
    float const epsilon = 1e-3f;

    for (int i = 0; i < particles.size(); ++i)
    {
        vec3& v = particles.speeds[i];
		vec3& p = particles.positions[i];
		vec3& f = particles.forces[i];

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;

        // small perturbation to avoid alignment
        if( p.y<=-1 ) {p.y = -1+rand_interval(-p.y - 1 + epsilon);  v.y *= -0.5f;}
        if( p.x<=-1 ) {p.x = -1+rand_interval(-p.x - 1 + epsilon);  v.x *= -0.5f;}
        if( p.x>=1 )  {p.x =  1-rand_interval(p.x - 1 + epsilon);  v.x *= -0.5f;}
        if( p.y>=1 ) { p.y = 1-rand_interval(p.y - 1 + epsilon);  v.y *= -0.5f;}
        if( p.z<=-1 ) { p.z = -1+rand_interval(-p.z - 1 + epsilon);  v.z *= -0.5f;}
        if( p.z>=1 )  { p.z =  1-rand_interval(p.z - 1 + epsilon);  v.z *= -0.5f; }
    }
}