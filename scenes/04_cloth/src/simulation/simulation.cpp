#include "simulation.hpp"

using namespace cgp;




// Fill value of force applied on each particle
// - Gravity
// - Drag
// - Spring force
// - Wind force
void simulation_compute_force(cloth_structure& cloth, simulation_parameters const& parameters)
{
    // Direct access to the variables
    //  Note: A grid_2D is a structure you can access using its 2d-local index coordinates as grid_2d(k1,k2)
    //   The index corresponding to grid_2d(k1,k2) is k1 + N1*k2, with N1 the first dimension of the grid.
    //   
    grid_2D<vec3>& force = cloth.force;  // Storage for the forces exerted on each vertex

    grid_2D<vec3> const& position = cloth.position;  // Storage for the positions of the vertices
    grid_2D<vec3> const& velocity = cloth.velocity;  // Storage for the normals of the vertices
    grid_2D<vec3> const& normal = cloth.normal;      // Storage for the velocity of the vertices
    

    size_t const N_total = cloth.position.size();       // total number of vertices
    size_t const N = cloth.N_samples();                 // number of vertices in one dimension of the grid

    // Retrieve simulation parameter
    //  The default value of the simulation parameters are defined in simulation.hpp
    float const K = parameters.K;              // spring stifness
    float const m = parameters.mass_total / N_total; // mass of a particle
    float const mu = parameters.mu;            // damping/friction coefficient
    float const	L0 = 1.0f / (N - 1.0f);        // rest length between two direct neighboring particle


    // Gravity
    const vec3 g = { 0,0,-9.81f };
    for (int ku = 0; ku < N; ++ku)
        for (int kv = 0; kv < N; ++kv)
            force(ku, kv) = m * g;

    // Drag (= friction)
    for (int ku = 0; ku < N; ++ku)
        for (int kv = 0; kv < N; ++kv)
            force(ku, kv) += -mu * m * velocity(ku, kv);

    // Spring forces
    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {
            // force(ku,kv) = ... fill here the force exerted by all the springs attached to the vertex at coordinates (ku,kv).
            // Notes:
            //   - The vertex positions can be accessed as position(ku,kv)
            //   - The neighbors are at position(ku+1,kv), position(ku-1,kv), position(ku,kv+1), etc. when ku+offset is still in the grid dimension.
            //   - You may want to loop over all the neighbors of a vertex to add each contributing force to this vertex
            //   - To void repetitions and limit the need of debuging, it may be a good idea to define a generic function that computes the spring force between two positions given the parameters K and L0
            //   - If the simulation is a bit too slow, you can speed it up in adapting the parameter N_step in scene.cpp that loops over several simulation step between two displays.
            for (int i = -1; i <= 1; ++i)
            {
                for (int j = -1; j <= 1; ++j)
                {
                    if (!(ku + i >= N or ku + i < 0 or kv + j >= N or kv + j < 0))
                    {
                        if (! (i == j or i == -j))
                            force(ku, kv) += K * (norm(position(ku + i, kv + j) - position(ku, kv)) - L0) *
                                    ((position(ku + i, kv + j) - position(ku, kv)) /
                                    norm(position(ku + i, kv + j) - position(ku, kv)));
                        else if (! (i == 0 and j == 0))
                            force(ku, kv) += K * (norm(position(ku + i, kv + j) - position(ku, kv)) - sqrt(2 * L0 * L0)) *
                                    ((position(ku + i, kv + j) - position(ku, kv)) /
                                    norm(position(ku + i, kv + j) - position(ku, kv)));
                    }
                }
            }
            for (int i = -2; i <= 2; i+=2)
            {
                for (int j = -2; j <= 2; j+=2)
                {
                    if (! (i == j or i == -j or ku + i >= N or ku + i < 0 or kv + j >= N or kv + j < 0))
                    {
                        force(ku, kv) += K * (norm(position(ku + i, kv + j) - position(ku, kv)) - 2 * L0) *
                                ((position(ku + i, kv + j) - position(ku, kv)) /
                                norm(position(ku + i, kv + j) - position(ku, kv)));
                    }
                }
            }
            force(ku, kv) += parameters.wind.magnitude * dot(parameters.wind.direction, normal(ku, kv)) * normal(ku,kv);
        }
    }
}

void simulation_numerical_integration(cloth_structure& cloth, simulation_parameters const& parameters, float dt)
{
    int const N = cloth.N_samples();
    int const N_total = cloth.position.size();
    float const m = parameters.mass_total/ static_cast<float>(N_total);

    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {
            vec3& v = cloth.velocity(ku, kv);
            vec3& p = cloth.position(ku, kv);
            vec3 const& f = cloth.force(ku, kv);

            // Standard semi-implicit numerical integration
            v = v + dt * f / m;
            p = p + dt * v;
        }
    }
    
}

void simulation_apply_constraints(cloth_structure& cloth, constraint_structure const& constraint)
{
    // Fixed positions of the cloth
    for (auto const& it : constraint.fixed_sample) {
        position_contraint c = it.second;
        cloth.position(c.ku, c.kv) = c.position; // set the position to the fixed one
    }

    grid_2D<vec3>& force = cloth.force;  // Storage for the forces exerted on each vertex

    grid_2D<vec3> & position = cloth.position;  // Storage for the positions of the vertices
    grid_2D<vec3> & velocity = cloth.velocity;  // Storage for the normals of the vertices
    grid_2D<vec3> & normal = cloth.normal;      // Storage for the velocity of the vertices

    size_t const N = cloth.N_samples();                 // number of vertices in one dimension of the grid

    // ToDo: apply external constraints
    // For all vertex:
    //   If vertex is below floor level ...
    //   If vertex is inside collision sphere ...
    for (int ku = 0; ku < N; ++ku) {
        for (int kv = 0; kv < N; ++kv) {
            if (position(ku, kv).z < constraint.ground_z + 0.01f)
            {
                position(ku, kv).z = constraint.ground_z + 0.01f;
            }
            if (norm(position(ku,kv) - constraint.sphere.center) < constraint.sphere.radius + 0.01f)
            {
                vec3 n = (position(ku,kv) - constraint.sphere.center) / norm(position(ku,kv) - constraint.sphere.center);
                position(ku, kv) = constraint.sphere.center + (constraint.sphere.radius + 0.01f) * n;
            }
        }
    }
}



bool simulation_detect_divergence(cloth_structure const& cloth)
{
    bool simulation_diverged = false;
    const size_t N = cloth.position.size();
    for (size_t k = 0; simulation_diverged == false && k < N; ++k)
    {
        const float f = norm(cloth.force.data.at_unsafe(k));
        const vec3& p = cloth.position.data.at_unsafe(k);

        if (std::isnan(f)) // detect NaN in force
        {
            std::cout << "\n **** NaN detected in forces" << std::endl;
            simulation_diverged = true;
        }

        if (f > 600.0f) // detect strong force magnitude
        {
            std::cout << "\n **** Warning : Strong force magnitude detected " << f << " at vertex " << k << " ****" << std::endl;
            simulation_diverged = true;
        }

        if (std::isnan(p.x) || std::isnan(p.y) || std::isnan(p.z)) // detect NaN in position
        {
            std::cout << "\n **** NaN detected in positions" << std::endl;
            simulation_diverged = true;
        }
    }

    return simulation_diverged;
}

