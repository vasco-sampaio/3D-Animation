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

std::vector<int> get_neighbors(vec3 pos, cgp::grid_3D<numarray<int>> grid, float h, int val)
{
    std::vector<int> neighbors_indexes;

    for (int i = -1; i <= 1; ++i)
    {
        for (int j = -1; j <= 1; ++j)
        {
            for (int k = -1; k <= 1; ++k)
            {
                if (i + int(pos.x/h) >= 0 && i + int(pos.x/h) <= grid.dimension.x - 1 &&
                j + int(pos.y/h) >= 0 && j + int(pos.y/h) <= grid.dimension.y - 1 &&
                k + int(pos.z/h) >= 0 && k + int(pos.z/h) <= grid.dimension.z - 1)
                {
                    numarray<int> grid_cell = grid.data[grid.index_to_offset(i + int(pos.x/h), j + int(pos.y/h),
                                                        k + int(pos.z/h))];
                    for (int & o : grid_cell)
                    {
                        if (o != val)
                            neighbors_indexes.push_back(o);
                    }

                }
            }
        }
    }

    return neighbors_indexes;
}

void update_density(numarray<particle_element>& particles, float h, float m, cgp::grid_3D<numarray<int>> grid)
{
    //  Compute the density value (particles[i].rho) at each particle position
    //  rho_i = \sum_j m W_density(pi,pj)
    int const N = particles.size();
    for(int i=0; i<N; ++i) {
        particles[i].rho = 0;
        std::vector<int> neighbors_indexes = get_neighbors(particles[i].p, grid, h, i);
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
void update_force(numarray<particle_element>& particles, float h, float m, float nu, cgp::grid_3D<numarray<int>> grid)
{
	// gravity
    const int N = particles.size();
    for(int i=0; i<N; ++i)
        particles[i].f = m * vec3{0,-9.81f,0};

    //To Do
    // For all particles i
    //   Compute F_pressure
    //   Compute F_viscosity
    //   particles[i].f += (F_pressure + F_viscosity)
    vec3 F_pressure;
    vec3 F_viscosity;
    for (int i = 0; i < N; ++i)
    {
        vec3 sum_press = {0,0,0};
        vec3 sum_visc = {0,0,0};

        particle_element particle = particles[i];
        std::vector<int> neighbors_indexes = get_neighbors(particle.p, grid, h, i);
        for (int j : neighbors_indexes)
        {
            if (norm(particles[i].p - particles[j].p) <= h)
            {
                sum_press += m * ((particles[i].pressure + particles[j].pressure) / (2 * particles[j].rho)) *
                             W_gradient_pressure(particles[i].p, particles[j].p, h);

                sum_visc += m * ((particles[j].v - particles[i].v) / (particles[j].rho)) *
                        W_laplacian_viscosity(particles[i].p, particles[j].p, h);
            }
        }
        F_pressure = (-m / (particles[i].rho)) * sum_press;
        F_viscosity = m * nu * sum_visc;
        particles[i].f += F_pressure + F_viscosity;
    }
}

void simulate(float dt, numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters, cgp::grid_3D<numarray<int>> grid)
{

	// Update values
	float h = sph_parameters.h;
	float c = sph_parameters.c;
	update_density(particles, h, sph_parameters.m, grid);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, h, sph_parameters.m, sph_parameters.nu, grid);  // Update forces

	// Numerical integration
	float const damping = 0.005f;
    int const N = particles.size();
	float const m = sph_parameters.m;
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}


	for (auto& vec : grid.data)
	    vec.clear();

	// Collision and update grid
    float const epsilon = 1e-3f;
    for(int k=0; k<N; ++k)
    {
        particle_element particle = particles[k];

        vec3& p = particle.p;
        vec3& v = particle.v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_interval();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_interval();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_interval();  v.x *= -0.5f;}
        if( p.y>1 ) {p.y = 1-epsilon*rand_interval();  v.y *= -0.5f;}
        if( p.z<-1 ) {p.z = -1+epsilon*rand_interval();  v.z *= -0.5f;}
        if( p.z>1 )  {p.z =  1-epsilon*rand_interval();  v.z *= -0.5f;}

        grid.data[grid.index_to_offset(((particle.p.x + 1)/h), ((particle.p.y + 1)/h), ((particle.p.z + 1)/h))].push_back(k);
    }
}