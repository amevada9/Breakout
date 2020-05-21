#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "list.h"
#include "vector.h"
#include "list.h"
#include "collision.h"

const double NINETY_DEGREES = 1.5070796327;
const double NEGATE = -1;
const int TYPE_MIN = 1;
const int TYPE_MAX = 2;
const int LARGE_NUMBER = 99999;

bool get_if_collided(collision_info_t collision_info) {
  return collision_info.collided;
}

vector_t get_collision_axis(collision_info_t collision_info) {
    return collision_info.axis;
}

void set_collision_axis(collision_info_t collision_info, vector_t axis) {
  collision_info.axis = axis;
}

void set_collided(collision_info_t collision_info, bool collided) {
  collision_info.collided = collided;
}

//Will return the min of the list or max of the list depending on what
// type is passed in
double find_extrema(list_t *numbers, int type){
  double *extrema = list_get(numbers, 0);
  for (int i = 1; i < list_size(numbers); i++){
    double *temp = list_get(numbers, i);
    if(type == 1 && *temp < *extrema){
      extrema = temp;
    }
    if(type == 2 && *temp > *extrema){
      extrema = temp;
    }
  }
  return *extrema;
}

/**
* finds the projection of each point in the provided shape onto the line given
* it then adds the magnitude of that projection to the list given
**/
void add_mag(vector_t line, list_t *shape, list_t *magnitude){
  for(int i = 0; i < list_size(shape); i++){
    vector_t *point = list_get(shape, i);
    vector_t projection =  vec_projection(*point, line);
    double *mag = malloc(sizeof(double));
    *mag = vec_magnitude(projection);
    //this corrects for projections that are negative
    if(vec_dot(*point, line) < 0){
      *mag *= NEGATE;
    }
    list_add(magnitude, mag);
  }
}

/**
*takes in a list that will contain the perpendicular vectors we are searching
* for and adds the vectors it gets from the given shape.
**/
void perpendicular_lines(list_t *perp_vectors, list_t *shape){
  size_t size = list_size(shape);
  for(size_t i = 0; i < size; i++){
    vector_t *point1 = list_get(shape, i);
    vector_t *point2 = list_get(shape, (i + 1) % size);
    vector_t difference = vec_subtract(*point1, *point2);
    vector_t *to_add = malloc(sizeof(vector_t));
    *to_add = vec_rotate(difference, NINETY_DEGREES);
    list_add(perp_vectors, to_add);
  }
}

double overlap(vector_t line, list_t *shape1, list_t *shape2){
  //lists of doubles that represent the magnitude of the projected vectors
  list_t *mag1 = list_init(list_size(shape1), (free_func_t) free);
  list_t *mag2 = list_init(list_size(shape2), (free_func_t) free);

  add_mag(line, shape1, mag1);
  add_mag(line, shape2, mag2);

  double max1 = find_extrema(mag1, TYPE_MAX);
  double max2 = find_extrema(mag2, TYPE_MAX);
  double min1 = find_extrema(mag1, TYPE_MIN);
  double min2 = find_extrema(mag2, TYPE_MIN);

  list_free(mag1);
  list_free(mag2);

  if(min1 <= min2){
    return(max1 - min2);
  }
  return(max2 - min1);
}


collision_info_t find_collision(body_t *body1, body_t *body2){
  list_t *shape1 = body_get_shape(body1);
  list_t *shape2 = body_get_shape(body2);
  collision_info_t information = {
    .collided = false,
    .axis = {-1, -1},
  };
  int size_of_perp = list_size(shape1) + list_size(shape2);
    //create a list of vectors that are the perpendicular lines and add to it
  list_t *perp_vectors = list_init(size_of_perp, (free_func_t) vec_free);
  perpendicular_lines(perp_vectors, shape1);
  perpendicular_lines(perp_vectors, shape2);
  vector_t overlap_vec = vec_init(0, 0);
  double least_overlap = LARGE_NUMBER;
  for(int i = 0; i < list_size(perp_vectors); i++){
    //if they do not overlap, return false because they are not colliding
    vector_t *line = list_get(perp_vectors, i);
    double temp_overlap = overlap(*line, shape1, shape2);
    if(temp_overlap < 0){
      list_free(perp_vectors);
      body_collided(body1, false);
      body_collided(body2, false);
      return information;
    }
    else if(temp_overlap < least_overlap){
      overlap_vec = *line;
      least_overlap = temp_overlap;
    }
  }
  vector_t overlap_axis = vec_unit(overlap_vec);
  information.collided = true;
  information.axis = overlap_axis;
  list_free(perp_vectors);
  return information;
}
