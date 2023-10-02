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

int3 get_cell_coordinates(int cell, int3 dimension)
{
    int3 cell_coordinates;
    cell_coordinates.x = cell % dimension.x;
    cell_coordinates.y = (cell / dimension.x) % dimension.y;
    cell_coordinates.z = cell / (dimension.x * dimension.y);
    return cell_coordinates;
}

std::vector<int> get_neighbors(const numarray<particle_element>& particles, const cgp::grid_3D<std::vector<int>>& grid, float h, int idx)
{
    // TO DO: Add neighbors property to particle_element to avoid calling it every time
    std::vector<int> neighbors_indexes;
    int3 dimensions = grid.dimension;
    int3 cell_coords = get_cell_coordinates(particles[idx].cell, dimensions);

    for (int i = std::max(0, cell_coords.x - 1); i <= std::min(dimensions.x - 1, cell_coords.x + 1); ++i)
    {
        for (int j = std::max(0, cell_coords.y - 1); j <= std::min(dimensions.y - 1, cell_coords.y + 1); ++j)
        {
            for (int k = std::max(0, cell_coords.z - 1); k <= std::min(dimensions.z - 1, cell_coords.z + 1); ++k)
            {
                for (auto& p : grid.data[grid.index_to_offset(i, j, k)])
                {
                    if (idx != p && norm(particles[idx].p - particles[p].p) <= h) {
                        neighbors_indexes.push_back(p);
                    }
                }
            }
        }
    }

    return neighbors_indexes;
}

void update_density(numarray<particle_element>& particles, float h, float m, cgp::grid_3D<std::vector<int>>& grid)
{
    //  Compute the density value (particles[i].rho) at each particle position
    //  rho_i = \sum_j m W_density(pi,pj)
    int const N = particles.size();
    for(int i=0; i<N; ++i) {
        // rho = density applied to self + sum of neighbors density applied to self
        particles[i].rho = m * W_density(particles[i].p, particles[i].p, h);
        
        std::vector<int> neighbors_indexes = get_neighbors(particles, grid, h, i);
        for (int j : neighbors_indexes)
        {
            if (norm(particles[i].p - particles[j].p) <= h)
                particles[i].rho += m * W_density(particles[i].p, particles[j].p, h);
        }
    }
}

// Convert the particle density to pressure
void update_pressure(numarray<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();
    for(int i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

// Compute the forces and update the acceleration of the particles
void update_force(numarray<particle_element>& particles, float h, float m, float nu, cgp::grid_3D<std::vector<int>>& grid)
{
    // For all particles i
    //   Compute F_pressure
    //   Compute F_viscosity
    //   particles[i].f = F_gravity + (F_pressure + F_viscosity)
    const int N = particles.size();
    vec3 F_pressure;
    vec3 F_viscosity;
    vec3 sum_press;
    vec3 sum_visc;
    for (int i = 0; i < N; ++i)
    {
        sum_press = vec3{0,0,0};
        sum_visc = vec3{0,0,0};

        std::vector<int> neighbors_indexes = get_neighbors(particles, grid, h, i);
        for (int j : neighbors_indexes)
        {
            if (norm(particles[i].p - particles[j].p) <= h)
            {
                sum_press += m * ((particles[i].pressure + particles[j].pressure) / (2.0f * particles[j].rho)) *
                             W_gradient_pressure(particles[i].p, particles[j].p, h);

                sum_visc += m * ((particles[j].v - particles[i].v) / (particles[j].rho)) *
                        W_laplacian_viscosity(particles[i].p, particles[j].p, h);
            }
        }

        F_pressure = float(-m / particles[i].rho) * sum_press;
        F_viscosity = m * nu * sum_visc;
        particles[i].f = m * vec3{0,-9.81f,0} + F_pressure + F_viscosity;
    }
}

void simulate(float dt, numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters, cgp::grid_3D<std::vector<int>>& grid)
{
	// Update values
	float h = sph_parameters.h;
	update_density(particles, h, sph_parameters.m, grid);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, h, sph_parameters.m, sph_parameters.nu, grid);  // Update forces

    // Numerical integration
	float const damping = 0.005f;
    int const N = particles.size();
	float const m = sph_parameters.m;
	for(int k=0; k<N; ++k)
	{
        vec3& v = particles[k].v;
		vec3& p = particles[k].p;
		vec3& f = particles[k].f;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}

    // Collision and update grid
    float const epsilon = 1e-3f;
    int new_cell;
    for(int k=0; k<N; ++k)
    {
        particle_element& particle = particles[k];

        vec3& p = particle.p;
        vec3& v = particle.v;

        // small perturbation to avoid alignment
        if( p.y<=-1 ) {p.y = -1+rand_interval(-p.y - 1 + epsilon);  v.y *= -0.5f;}
        if( p.x<=-1 ) {p.x = -1+rand_interval(-p.x - 1 + epsilon);  v.x *= -0.5f;}
        if( p.x>=1 )  {p.x =  1-rand_interval(p.x - 1 + epsilon);  v.x *= -0.5f;}
        if( p.y>=1 ) { p.y = 1-rand_interval(p.y - 1 + epsilon);  v.y *= -0.5f;}
        if( p.z<=-1 ) { p.z = -1+rand_interval(-p.z - 1 + epsilon);  v.z *= -0.5f;}
        if( p.z>=1 )  { p.z =  1-rand_interval(p.z - 1 + epsilon);  v.z *= -0.5f; }

        // update grid
        new_cell = grid.index_to_offset(floor((p.x + 1)/h), floor((p.y + 1)/h), floor((p.z + 1)/h));
        if (particle.cell != new_cell)
        {
            grid.data[particle.cell].erase(std::find(grid.data[particle.cell].begin(), grid.data[particle.cell].end(), k), grid.data[particle.cell].end());
            particle.cell = new_cell;
            grid.data[particle.cell].push_back(k);
        }
    }
}