#pragma once

#include "cgp/cgp.hpp"
#include "environment.hpp"

using cgp::mesh_drawable;
using cgp::vec3;


struct gui_parameters {
	bool display_frame     = true;
	bool is_full_screen = false;
	bool display_wireframe = false;
};



// The structure of the custom scene
struct scene_structure : scene_inputs_generic {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	camera_controller_orbit_euler camera_control;
	camera_projection_perspective camera_projection;
	window_structure window;


	mesh_drawable global_frame;        // The standard global frame
	environment_structure environment; // Standard environment controler
	input_devices inputs; // Storage for inputs status (mouse, keyboard, window dimension)
	gui_parameters gui;                // Standard GUI element storage
	
	cgp::timer_basic timer;
	mesh_drawable cube; 
	mesh_drawable ground;
	mesh_drawable sphere;

	mesh_drawable camel;

	// (In the scene_structure header)
	std::vector<vec3> cone_positions;
	mesh_drawable cone;

	std::vector<vec3> cylinder_positions;
	mesh_drawable cylinder;

	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();  // Standard initialization to be called before the animation loop
	void display_frame();     // The frame display to be called within the animation loop
	void display_gui(); // The display of the GUI, also called within the animation loop

	void mouse_move_event();
	void mouse_click_event();
	void keyboard_event();
	void idle_frame();

};





