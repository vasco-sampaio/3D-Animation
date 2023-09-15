#include "scene.hpp"



using namespace cgp;

void scene_structure::initialize()
{
	// Basic set-up
	// ***************************************** //
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.set_rotation_axis_z();
	camera_control.look_at({ 2.0f,-2.0f,1.0f }, { 0,0,0 });

	global_frame.initialize_data_on_gpu(mesh_primitive_frame());


	// Create the hierarchy
	// ************************************ //

	// Initialize the temporary mesh_drawable that will be inserted in the hierarchy
	// Body Features
	mesh_drawable body;
	mesh_drawable head;
	mesh_drawable left_foot;
	mesh_drawable right_foot;
	mesh_drawable left_hand;
	mesh_drawable right_hand;

	// Face Features
	// Eyes
	mesh_drawable left_eye;
	mesh_drawable right_eye;
	mesh_drawable left_pupil;
	mesh_drawable right_pupil;

	// Ears
	mesh_drawable left_ear;
	mesh_drawable right_ear;

	// Nose
	mesh_drawable snout;
	mesh_drawable nose;

	// Create the geometry of the meshes
	body.initialize_data_on_gpu(mesh_primitive_ellipsoid({1.0f, 0.7f, 0.55f}));
	body.model.scaling = 1.0f;
	head.initialize_data_on_gpu(mesh_primitive_sphere()); head.model.scaling = 0.45f;
	left_foot.initialize_data_on_gpu(mesh_primitive_ellipsoid({0.6f, 0.7f, 1.0f}));
	left_foot.model.scaling = 0.4f;
	right_foot.initialize_data_on_gpu(mesh_primitive_ellipsoid({0.6f, 0.7f, 1.0f}));
	right_foot.model.scaling = 0.4f;
	left_hand.initialize_data_on_gpu(mesh_primitive_ellipsoid({0.6f, 0.7f, 1.0f}));
	left_hand.model.scaling = 0.4f;
	right_hand.initialize_data_on_gpu(mesh_primitive_ellipsoid({0.6f, 0.7f, 1.0f}));
	right_hand.model.scaling = 0.4f;

	left_pupil.initialize_data_on_gpu(mesh_primitive_sphere(0.08f));
	right_pupil.initialize_data_on_gpu(mesh_primitive_sphere(0.08f));
	left_eye.initialize_data_on_gpu(mesh_primitive_sphere(0.12f));
	right_eye.initialize_data_on_gpu(mesh_primitive_sphere(0.12f));

	left_ear.initialize_data_on_gpu(mesh_primitive_ellipsoid({0.7f, 0.7f, 0.3f}));
	left_ear.model.scaling = 0.3f;
	right_ear.initialize_data_on_gpu(mesh_primitive_ellipsoid({0.7f, 0.7f, 0.3f}));
	right_ear.model.scaling = 0.3f;

	// TODO: FINISH SNOUT AND NOSE
	snout.initialize_data_on_gpu(mesh_primitive_ellipsoid({0.7f, 0.7f, 0.3f}));

	// Set the color of some elements
	vec3 brown1 = { 0.7f, 0.4f, 0 };
	body.material.color = brown1;
	head.material.color = brown1;
	left_foot.material.color = brown1;
	right_foot.material.color = brown1;
	left_hand.material.color = brown1;
	right_hand.material.color = brown1;
	left_ear.material.color = brown1;
	right_ear.material.color = brown1;

	vec3 black = {0,0,0};
	left_pupil.material.color = black;
	right_pupil.material.color = black;

	vec3 white = {1.0f,1.0f,1.0f};
	left_eye.material.color = white;
	right_eye.material.color = white;

	// Add the elements in the hierarchy
	//   The syntax is hierarchy.add(mesh_drawable, "name of the parent element", [optional: local translation in the hierarchy])
	//   Notes: 
	//     - An element must necessarily be added after its parent
	//     - The first element (without explicit name of its parent) is assumed to be the root.
	hierarchy.add(body, "Body");
	hierarchy.add(head, "Head", "Body", {1.2f,0,0});
	hierarchy.add(left_foot, "Left foot", "Body", {-0.6f, -0.45f, 0.4f});
	hierarchy.add(right_foot, "Right foot", "Body", {-0.6f, 0.45f, 0.4f});
	hierarchy.add(left_hand, "Left hand", "Body", {0.6f, -0.45f, 0.4f});
	hierarchy.add(right_hand, "Right hand", "Body", {0.6f, 0.45f, 0.4f});

	hierarchy.add(left_eye, "Left eye", "Head", {0.15f, -0.13f, 0.3f});
	hierarchy.add(right_eye, "Right eye", "Head", {0.15f, 0.13f, 0.3f});
	hierarchy.add(left_pupil, "Left pupil", "Left eye", {0, 0, 0.05f});
	hierarchy.add(right_pupil, "Right pupil", "Right eye", {0, 0, 0.05f});

	hierarchy.add(left_ear, "Left ear", "Head", {0.37f, -0.25f, 0.1f});
	hierarchy.add(right_ear, "Right ear", "Head", {0.37f, 0.25f, 0.1f});

}





void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	if (gui.display_frame)
		draw(global_frame, environment);

	// Update the current time
	timer.update();

	// This function must be called before the drawing in order to propagate the deformations through the hierarchy
	hierarchy.update_local_to_global_coordinates();

	// Draw the hierarchy as a single mesh
	draw(hierarchy, environment);
	if (gui.display_wireframe)
		draw_wireframe(hierarchy, environment);

}




void scene_structure::display_gui()
{
	ImGui::Checkbox("Frame", &gui.display_frame);
	ImGui::Checkbox("Wireframe", &gui.display_wireframe);
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
