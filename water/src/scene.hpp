#pragma once


#include "cgp/cgp.hpp"
#include "environment.hpp"

#include "simulation/simulation.hpp"
#include "octree/octree.hpp"

using cgp::mesh_drawable;


struct gui_parameters {
	bool display_color = true;
	bool display_particles = true;
	bool display_radius = false;
};

// The structure of the custom scene
struct scene_structure : cgp::scene_inputs_generic {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	camera_controller_orbit camera_control;
	camera_projection_orthographic camera_projection;
	window_structure window;

	mesh_drawable global_frame;          // The standard global frame
	environment_structure environment;   // Standard environment controler
	input_devices inputs;                // Storage for inputs status (mouse, keyboard, window dimension)
	gui_parameters gui;                  // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_basic timer;

	sph_parameters_structure sph_parameters; // Physical parameter related to SPH
	cgp::mesh_drawable sphere_particle; // Sphere used to display a particle
	cgp::curve_drawable curve_visual;   // Circle used to display the radius h of influence

    unsigned int number_of_particles = 10000;
	Octree<std::set<int>> grid = Octree<std::set<int>>(std::floor(std::log2(static_cast<double>(number_of_particles)) / 3.0)); // grid used to store the particles
    ParticleArray particles = ParticleArray(number_of_particles); // Array of particles

	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();    // Standard initialization to be called before the animation loop
	void display_frame(); // The frame display to be called within the animation loop
	void display_gui();   // The display of the GUI, also called within the animation loop

	void initialize_sph();

	void mouse_move_event();
	void mouse_click_event();
	void keyboard_event();
	void idle_frame();

};





