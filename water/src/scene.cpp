#include <algorithm>

#include "scene.hpp"

using namespace cgp;

void scene_structure::initialize()
{
	camera_projection = camera_projection_orthographic{ -1.1f, 1.1f, -1.1f, 1.1f, -10, 10, window.aspect_ratio() };
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.look_at({ 0.0f, 0.0f, 2.0f }, {0,0,0}, {0,1,0});
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());

	initialize_sph();
	sphere_particle.initialize_data_on_gpu(mesh_primitive_sphere(1.0,{0,0,0},10,10));
    sphere_particle.material.color = { 0,0,1 };
	sphere_particle.model.scaling = 0.01f;
	curve_visual.color = { 1,0,0 };
	curve_visual.initialize_data_on_gpu(curve_primitive_circle());
}

void scene_structure::initialize_sph()
{
    for (int k = 0; k < number_of_particles; k++) {
        float x = rand_interval(-1.0f, 1.0f);
        float y = rand_interval(-1.0f, 1.0f);
        float z = rand_interval(-1.0f, 1.0f);

        cgp::vec3 p = {x, y, z};
        unsigned int octant = grid.get_octant(p);

        particles.push_back(p, {0, 0, 0}, {0, 0, 0}, 0, 0, octant);
        grid.insert_into_node(octant, k);
    }

    std::cout << "Number of particles: " << particles.size() << std::endl;
}


float t = 0.0f;

void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	timer.update(); // update the timer to the current elapsed time
	float const dt = 0.005f * timer.scale;

	if (timer.t - t > 0.5f) {
        update_grid(particles, grid);
		t = timer.t;
	}

    simulate(dt, particles, sph_parameters, grid);

	if (gui.display_particles) {
		for (int k = 0; k < particles.size(); ++k) {
			sphere_particle.model.translation = particles.positions[k];
			draw(sphere_particle, environment);
		}
	}

	if (gui.display_radius) {
		curve_visual.model.scaling = sph_parameters.h;
		for (int k = 0; k < particles.size(); k += 10) {
			curve_visual.model.translation = particles.positions[k];
			draw(curve_visual, environment);
		}
	}

}

void scene_structure::display_gui()
{
	ImGui::SliderFloat("Timer scale", &timer.scale, 0.01f, 4.0f, "%0.2f");

	bool const restart = ImGui::Button("Restart");
	if (restart)
		initialize_sph();

	ImGui::Checkbox("Color", &gui.display_color);
	ImGui::Checkbox("Particles", &gui.display_particles);
	ImGui::Checkbox("Radius", &gui.display_radius);
}

void scene_structure::mouse_move_event()
{
	if (!inputs.keyboard.shift)
		camera_control.action_mouse_move(environment.camera_view);
}
void scene_structure::mouse_click_event()
{
	camera_control.action_mouse_click(environment.camera_view);
}
void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
}

