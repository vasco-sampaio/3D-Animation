#include "simulation.hpp"

using namespace cgp;



void simulate(std::vector<particle_structure>& particles, float dt)
{
	vec3 const g = { 0,0,-9.81f };
	size_t const N = particles.size();
	for (size_t k = 0; k < N; ++k)
	{
		particle_structure& particle = particles[k];

		vec3 const f = particle.m * g;

		particle.v = (1 - 0.9f * dt) * particle.v + dt * f;
		particle.p = particle.p + dt * particle.v;
	}

	// To do :

	manage_sphere_collisions(particles);
	manage_cube_collisions(particles);
}


void manage_cube_collisions(std::vector<particle_structure>& particles)
{
    numarray<vec3> cube_normals = {{-1, 0, 0}, {1, 0, 0},
                                   {0, -1, 0}, {0, 1, 0},
                                   {0, 0, -1}, {0, 0, 1}};

    numarray<vec3> cube_faces = {{1, 0, 0}, {-1, 0, 0},
                                 {0, 1, 0}, {0, -1, 0},
                                 {0, 0, 1}, {0, 0, -1}};

    for (int i = 0; i < particles.size(); ++i)
    {
        for (int j = 0; j < cube_faces.size(); ++j)
        {
            float detection = dot(particles[i].p - cube_faces[j], cube_normals[j]);
            if (detection < particles[i].r)
            {
                vec3 v_orth = dot(particles[i].v, cube_normals[j]) * cube_normals[j];
                vec3 v_col = particles[i].v - v_orth;
                float d = particles[i].r - dot((particles[i].p - cube_faces[j]), cube_normals[j]);

                particles[i].p += d * cube_normals[j];
                particles[i].v = 0.8f * v_col - 0.8f * v_orth;
            }
        }
    }
}

void manage_sphere_collisions(std::vector<particle_structure>& particles)
{
    for (int i = 0; i < particles.size() - 1; ++i)
    {
        for (int k = i + 1; k < particles.size(); ++k)
        {
            vec3 u = (particles[i].p - particles[k].p) / norm(particles[i].p - particles[k].p);
            float d = (particles[i].r + particles[k].r) - norm(particles[i].p - particles[k].p);
            if (norm(particles[i].p - particles[k].p) < particles[i].r + particles[k].r)
            {
                if (abs(norm(particles[i].v - particles[k].v)) > 0.1f)
                {
                    float j = dot(2 * ((particles[i].m * particles[k].m) / (particles[i].m + particles[k].m)) *
                            (particles[i].v - particles[k].v), u);

                    particles[i].v = particles[i].v + (j / particles[i].m) * u;
                    particles[k].v = particles[k].v - (j / particles[k].m) * u;

                    /*vec3 v_orth1 = dot(particles[i].v, -u) * -u;
                    vec3 v_col1 = particles[i].v - v_orth1;

                    vec3 v_orth2 = dot(particles[k].v, u) * u;
                    vec3 v_col2 = particles[k].v - v_orth2;

                    particles[i].v = 0.7f * v_col1 + 0.7f * v_orth2;
                    particles[k].v = 0.7f * v_col2 + 0.7f * v_orth1;*/
                }
                else
                {
                    particles[i].v /= 0.2f;
                    particles[k].v /= 0.2f;
                }
            }
            particles[i].p += d / (2 * u);
        }
    }
}

