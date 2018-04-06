#include "sphere.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cfloat>

/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection,
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/
 float intersect_board(Point o, Vector u, Spheres *sph, Point *hit) {
     // parametric ray equation: r = o + t.u
     Vector chess_board_norm = {0, 1, 0};
     if (vec_dot(chess_board_norm, u) == 0) return -1; // paralell
     float t = vec_dot(get_vec(o, sph->center), chess_board_norm) / (vec_dot(u, chess_board_norm));
     if (t < 0.00001) return -1.0;  //
     Point intersection = get_point(o, vec_scale(u, t));
     float x_length = 10.0;
     float z_length = -1.0;
     if (intersection.x <= x_length && intersection.x >= -x_length && intersection.z <= z_length && intersection.z >= -10) {
         if (intersection.x >= 0) {
             if ((int(intersection.x) + int(intersection.z)) % 2 == 0) {
                 for (int i = 0; i < 3; i++) {
                     sph->mat_ambient[i]  = 0;
                     sph->mat_diffuse[i]  = 0;
                     sph->mat_specular[i] = 0;
                 }
             } else {
                 for (int i = 0; i < 3; i++) {
                     sph->mat_ambient[i]  = 1;
                     sph->mat_diffuse[i]  = 1;
                     sph->mat_specular[i] = 1;
                 }
             }
         } else {
             if ((int(intersection.x) + int(intersection.z)) % 2 == 0) {
                 for (int i = 0; i < 3; i++) {
                     sph->mat_ambient[i]  = 1;
                     sph->mat_diffuse[i]  = 1;
                     sph->mat_specular[i] = 1;
                 }
             } else {
                 for (int i = 0; i < 3; i++) {
                     sph->mat_ambient[i]  = 0;
                     sph->mat_diffuse[i]  = 0;
                     sph->mat_specular[i] = 0;
                 }
             }
         }
         *hit = intersection;
         return t;
     }
     return -1.0;
 }

float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit) {
	if (sph->index > 0) {
		// parametric ray equation: r = o + t.u
		// at^2 + bt + c = 0
		float a = pow(vec_len(u), 2);
		float b = 2 * vec_dot(get_vec(sph->center, o), u);
		float c = pow(vec_len(get_vec(sph->center, o)), 2) - pow(sph->radius, 2);
		float delta = b * b - 4 * a * c;
		if (delta < 0) return -1.0;
		float t1 = (-b + sqrt(delta)) / (2 * a);
		float t2 = (-b - sqrt(delta)) / (2 * a);
		if (t1 < 0 && t2 < 0) return -1.0;

		Point first_intersection  = get_point(o, vec_scale(u, t1));
		Point second_intersection = get_point(o, vec_scale(u, t2));
        if (t2 >= t1 && t1 >= 0) {
            *hit = first_intersection;
            return t1;
        }
        if (t1 >= t2 && t2 >= 0) {
            *hit = second_intersection;
            return t2;
        }
    } else
        return intersect_board(o, u, sph, hit);
}
/*********************************************************************
 * This function returns 1 if a shadow ray from source o in direction u
 * is blocked by a sphere, otherwise returns 0
 **********************************************************************/
 // this one is only correct for the first shadow, need to be fixed
int is_all_light_blocked(Point o, Vector u, Spheres *sph, Spheres *hit) {
	Spheres * current_sphere = sph;
    Point temp;
    float current_hit_depth;
    int inside = false;
    if (vec_len(get_vec(hit->center, o)) < hit->radius) inside = true;
	while (current_sphere) {
        current_hit_depth = intersect_sphere(o, u, sph, &temp);
        // if (current_hit_depth != 1.0) return 1;
        if (current_hit_depth != -1.0) {
            if (current_sphere != hit) return 1;
            else {
                if (inside) {
                    if (current_hit_depth > 0.00001) return 1;
                }
                else
                    if (current_hit_depth < 0.00001) return 1;
                }
        }
	    current_sphere = current_sphere->next;
	}
	return 0;
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For example, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres *intersect_scene(Point o, Vector u, Spheres *sph, Point *hit, Spheres * skip) {
	Spheres * current_sphere  	= sph;
	Point temp;
	Spheres * first_hit_sphere	= NULL;
	float first_hit_depth	  	= FLT_MAX;
	float current_hit_depth;
	while (current_sphere) {
		current_hit_depth  = intersect_sphere(o, u, current_sphere, &temp);
		if (current_hit_depth != -1.0 && current_sphere != skip) {
			if (current_hit_depth < first_hit_depth) {
				first_hit_depth  = current_hit_depth;
				*hit 			 = temp;
				first_hit_sphere = current_sphere;
			}
		}
		current_sphere = current_sphere->next;
	}
	return first_hit_sphere;
}

/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine,
		    float refl, float trans, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->transparency = trans;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}
