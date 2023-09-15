#include "scene.hpp"

using namespace cgp;

void scene_structure::initialize()
{
	// Basic set-up
	// ***************************************** //
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.set_rotation_axis_z();
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());


	// Initialize the shapes of the scene
	// ***************************************** //

	// Create a mesh structure (here a cube)
	mesh cube_mesh = mesh_primitive_cube(/*center*/{ 0,0,0 }, /*edge length*/ 1.0f);
	// a mesh is simply a container of vertex data (position,normal,color,uv) and triangle index
	// the mesh data are stored on the CPU memory - they will need to be sent to the GPU memory before being drawn


	// Initialize a "mesh_drawable" from a mesh structure
	//   - mesh : store buffer of data (vertices, indices, etc) on the CPU. The mesh structure is convenient to manipulate in the C++ code but cannot be displayed directly (data is not on GPU).
	//   - mesh_drawable : store VBO associated to elements on the GPU + associated uniform parameters. A mesh_drawable can be displayed using the function draw(mesh_drawable, environment). It only stores the indices of the buffers on the GPU - the buffer of data cannot be directly accessed in the C++ code through a mesh_drawable.
	//   Note: A mesh_drawable can be initialized from a mesh structure using the syntax
	//    [meshDrawableName].initalize_data_on_gpu( [meshName] );
	float L = 3;
	float z_floor = -0.51f;

	float cone_radius = 0.1f;
	float cone_height = 0.2f;
	float cylinder_height = 0.3f;
	cone.initialize_data_on_gpu(mesh_primitive_cone(cone_radius, cone_height));
	cone.material.color = { 0,0.6f,0 }; // dark green

	cylinder.initialize_data_on_gpu(mesh_primitive_cylinder(0.05f,
                                                            { 0, 0, 0},
                                                            {0, 0, cylinder_height}));
	cylinder.material.color = { 0.4f,0,0 }; // brown

	const int N_cone = 70;
	for (int k = 0; k < N_cone; ++k)
	{
	    float x;
	    float y;
	    bool check = false;
	    int o = 0;

	    do {
	        x = rand_interval(-L + cone_radius, L - cone_radius);
	        y = rand_interval(-L + cone_radius, L - cone_radius);

	        for (int i = 0; i < cylinder_positions.size() and !check; ++i)
	        {
	            check = norm({cone_positions[i * 3].x - x, cone_positions[i * 3].y - y, 0}) > 0.1f;
	        }

	    } while (check and 5 > o++);

	    cone_positions.emplace_back(x, y, z_floor + cylinder_height);
	    cone_positions.emplace_back(x, y, z_floor + cylinder_height - 0.1f);
	    cone_positions.emplace_back(x, y, z_floor + cylinder_height - 0.2f);

	    cylinder_positions.emplace_back(x, y, z_floor);
	}


	// Same process for the ground which is a plane

	mesh ground_mesh = mesh_primitive_quadrangle({ -L,-L,z_floor }, { L,-L,z_floor }, { L,L,z_floor }, { -L,L,z_floor });
	ground.initialize_data_on_gpu(ground_mesh);
}



// This function is constantly called at every frame
void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	// the general syntax to display a mesh is:
	//   draw(mesh_drawableName, environment);
	//     Note: scene is used to set the uniform parameters associated to the camera, light, etc. to the shader
	draw(ground, environment);

	for (int k = 0; k < cylinder_positions.size(); ++k)
	{
	    cone.model.translation = cone_positions[3 * k];
	    draw(cone, environment);
	    cone.model.translation = cone_positions[3 * k + 1];
	    draw(cone, environment);
	    cone.model.translation = cone_positions[3 * k + 2];
	    draw(cone, environment);

	    cylinder.model.translation = cylinder_positions[k];
        draw(cylinder, environment);
	}
}


void scene_structure::display_gui()
{
	ImGui::Checkbox("Frame", &gui.display_frame);
}

void scene_structure::mouse_move_event()
{
	if (!inputs.keyboard.shift)
		camera_control.action_mouse_move(environment.camera_view);
}
void scene_structure::mouse_click_event()
{
	//camera_control.action_mouse_click(environment.camera);
}
void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
	
}