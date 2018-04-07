#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"
#include <algorithm>

//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;

// light 1 position and color
extern Point light1;
extern float light1_ambient[3];
extern float light1_diffuse[3];
extern float light1_specular[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int shadow_on;
extern int reflection_on;
extern int refraction_on;
extern int stochastic_on;
extern int supersampling_on;
extern int chess_board_on;
extern int step_max;

RGB_float get_ambient_component(Spheres *sph) {
	// ambient = global ambient + light1_ambient
	RGB_float ambient_component;
	ambient_component.r = (global_ambient[0] + light1_ambient[0]) * sph->reflectance;
	ambient_component.g = (global_ambient[1] + light1_ambient[1]) * sph->reflectance;
	ambient_component.b = (global_ambient[2] + light1_ambient[2]) * sph->reflectance;
	return ambient_component;
}

void create_chess_board() {
	Point board_ctr = {0, -2, -2};
    float board_rad = 1.0;
    float board_ambient[] = {0, 0, 0};
    float board_diffuse[] = {0, 0, 0};
    float board_specular[] = {0, 0, 0};
    float board_shineness = 20;
    float board_reflectance = 0.2;
    float board_transparency = 1.0;
    scene = add_sphere(scene, board_ctr, board_rad, board_ambient,
               board_diffuse, board_specular, board_shineness,
           board_reflectance, board_transparency, -1);
}
/*********************************************************************
* Phong illumination
  Precondition:
  - q, v: a surface point and direction vector
  - surf_norm: normal vector
  - sph: list of spheres
  Postcondition:returns color at q
*********************************************************************/
RGB_float phong(Point q, Vector v, Vector surf_norm, Spheres *sph) {
	RGB_float ambient_component = get_ambient_component(sph);
	Vector point_light 			= get_vec(q, light1);
	float distance_point_light 	= vec_len(point_light);
	float attenuate_coefficient = std::min(1.0, 1.0 / (decay_a + decay_b * distance_point_light + decay_c * pow(distance_point_light, 2)));

	// perfect diffuse component
	RGB_float diffuse_component = {0, 0, 0}; 	// default for reflection angle < 0
	normalize(&point_light);
	normalize(&surf_norm);
	float surf_norm_dot_point_light = vec_dot(surf_norm, point_light);
	if (surf_norm_dot_point_light >= 0) {
		diffuse_component.r = light1_diffuse[0] * sph->mat_diffuse[0] * surf_norm_dot_point_light;
		diffuse_component.g = light1_diffuse[1] * sph->mat_diffuse[1] * surf_norm_dot_point_light;
		diffuse_component.b = light1_diffuse[2] * sph->mat_diffuse[2] * surf_norm_dot_point_light;
	}

	// specular component
	RGB_float specular_component = {0, 0, 0};
	Vector perfect_relfection = vec_minus(vec_scale(surf_norm, surf_norm_dot_point_light * 2), point_light);
    float view_dot_reflection = vec_dot(v, perfect_relfection);
    if (view_dot_reflection >= 0) {
        specular_component.r = light1_specular[0] * sph->mat_specular[0] * pow(view_dot_reflection, sph->mat_shineness);
        specular_component.g = light1_specular[1] * sph->mat_specular[1] * pow(view_dot_reflection, sph->mat_shineness);
        specular_component.b = light1_specular[2] * sph->mat_specular[2] * pow(view_dot_reflection, sph->mat_shineness);
    }

	// sum up
	RGB_float color;
	color = clr_scale(clr_add(diffuse_component, specular_component), attenuate_coefficient);
	color = clr_add(color, ambient_component);
	return color;
}

/************************************************************************
* recursive ray tracer
  Precondition:
  - o, ray: a point o and direction vector
  - level_recursion: number of recursive rays
  - skip: sphere emanating the ray
  Postcondition: returns color value at the interseciton of the ray with a sphere
************************************************************************/
RGB_float recursive_ray_trace(Point o, Vector ray, int level_recursion, Spheres *skip) {
	RGB_float color = background_clr;
	RGB_float reflection_color, refraction_color;
	RGB_float stochastic_color = {0, 0, 0};
	// prepare for Phong illumination
	Point hit;
	Vector point_eye, point_light;
	Vector surf_norm = {0, 1, 0};
	Spheres * first_hit_sphere = intersect_scene(o, ray, scene, &hit, skip);
	if (first_hit_sphere) {
		if (first_hit_sphere->index > 0)
			surf_norm = sphere_normal(hit, first_hit_sphere);
		point_eye = get_vec(hit, o);
		point_light = get_vec(hit, light1);
		normalize(&point_eye);
		normalize(&surf_norm);
		normalize(&point_light);
		// shadow
		if (shadow_on && is_all_light_blocked(hit, point_light, scene, first_hit_sphere))
			color = get_ambient_component(first_hit_sphere);
		else color = phong(hit, point_eye, surf_norm, first_hit_sphere);
		// reflection
		if (reflection_on && level_recursion++ < step_max) {
			float surf_norm_dot_point_light = vec_dot(surf_norm, point_eye);
			Vector reflection = vec_minus(vec_scale(surf_norm, surf_norm_dot_point_light * 2), point_eye);
            normalize(&reflection);
			reflection_color = recursive_ray_trace(hit, reflection, level_recursion, first_hit_sphere);
			color = clr_add(color, clr_scale(reflection_color, first_hit_sphere->reflectance));
		}
        // transparency/ refraction : wikipedia
		if (refraction_on && level_recursion++ < step_max) {
			Vector light    = ray;
			Vector n        = surf_norm;
			float ratio = 1.0 / first_hit_sphere->transparency;
			float costheta1   = -1.0 * vec_dot(surf_norm, light);
			float sintheta2   = ratio * sqrt(1 - pow(costheta1, 2));
			float costheta2   = sqrt(1 - pow(sintheta2, 2));
			float parameter1  = ratio;
			float parameter2  = (ratio * costheta1 - costheta2);
			Vector refraction = vec_plus(vec_scale(light, parameter1), vec_scale(n, parameter2));
			refraction_color  = recursive_ray_trace(hit, refraction, level_recursion, first_hit_sphere);
			color = clr_add(color, clr_scale(refraction_color, 0.2)); // simply refractance = 0.2 for all materials
		}
		// diffuse reflection using stochastic
		if (stochastic_on && level_recursion++ < step_max) {
			// generate random reflection rays
			for (int i = 0; i < RAY; i++) {
				Vector random_ray = {0.0, 0.0, 0.0};
				do {
					// ensure no reflection to sphere itself
					random_ray.x = float(rand() % 20 - 10);
					random_ray.y = float(rand() % 20 - 10);
					random_ray.z = float(rand() % 20 - 10);
				} while (vec_dot(random_ray, surf_norm) <= 0);
				// add up colors from random rays
				RGB_float random_ray_color = recursive_ray_trace(hit, random_ray, level_recursion, first_hit_sphere);
				stochastic_color = clr_add(stochastic_color, random_ray_color);
			}
			stochastic_color = clr_scale(stochastic_color, 1 / RAY);
			color = clr_add(color, clr_scale(stochastic_color, first_hit_sphere->reflectance));
		}
	}
	return color;
}

/*********************************************************************
* This function traverses all the pixels and cast rays. It calls the
* recursive ray tracer and assign return color to frame
*********************************************************************/
void ray_trace() {
	int i, j;
	float x_grid_size = image_width / float(win_width);
	float y_grid_size = image_height / float(win_height);;
	float x_start = -0.5 * image_width;
	float y_start = -0.5 * image_height;
	RGB_float ret_color;
	Point cur_pixel_pos;
	Vector ray;

	// ray is cast through center of pixel
	cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
	cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
	cur_pixel_pos.z = image_plane;
	// chess_board
	if (chess_board_on) create_chess_board();
	for (i=0; i<win_height; i++) {
		for (j=0; j<win_width; j++) {
			ray = get_vec(eye_pos, cur_pixel_pos);
			normalize(&ray);
			ret_color = recursive_ray_trace(eye_pos, ray, 0, NULL);
			// supersampling
			if (supersampling_on) {
				Point pixel;
				pixel.x = cur_pixel_pos.x - 0.25 * x_grid_size;
				pixel.y = cur_pixel_pos.y - 0.25 * y_grid_size;
				pixel.z = image_plane;
				ray = get_vec(eye_pos, pixel);
				normalize(&ray);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, ray, 0, NULL));
				pixel.x = cur_pixel_pos.x + 0.25 * x_grid_size;
				pixel.y = cur_pixel_pos.y - 0.25 * y_grid_size;
				ray = get_vec(eye_pos, pixel);
				normalize(&ray);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, ray, 0, NULL));
				pixel.x = cur_pixel_pos.x + 0.25 * x_grid_size;
				pixel.y = cur_pixel_pos.y + 0.25 * y_grid_size;
				ray = get_vec(eye_pos, pixel);
				normalize(&ray);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, ray, 0, NULL));
				pixel.x = cur_pixel_pos.x - 0.25 * x_grid_size;
				pixel.y = cur_pixel_pos.y + 0.25 * y_grid_size;
				ray = get_vec(eye_pos, pixel);
				normalize(&ray);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, ray, 0, NULL));
			}
			frame[i][j][0] = GLfloat(ret_color.r);
			frame[i][j][1] = GLfloat(ret_color.g);
			frame[i][j][2] = GLfloat(ret_color.b);
			cur_pixel_pos.x += x_grid_size;
		}
		cur_pixel_pos.y += y_grid_size;
		cur_pixel_pos.x = x_start;
	}
}
