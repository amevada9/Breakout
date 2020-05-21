#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "body.h"
#include "polygon.h"

const double ACC_MULT = 0.5;
const size_t DEFAULT_SIZE = 10;

typedef struct body {
  list_t *points;
  vector_t velocity;
  vector_t acceleration;
  vector_t centroid;
  vector_t force;
  vector_t impulse;
  double mass;
  rgb_color_t color;
  double angle;
  void *info;
  free_func_t info_freer;
  bool removed;
  bool collided_with;
} body_t;

/**
 * Allocates memory for a body with the given parameters.
 * The body is initially at rest.
 */
body_t *body_init_with_info(list_t *shape, double mass, rgb_color_t color,
  void *info, free_func_t info_freer){
    body_t *body = malloc(sizeof(body_t));
    assert(body != NULL);
    body->points = list_init(DEFAULT_SIZE, (free_func_t) vec_free);
    assert(body->points != NULL);
    body->points = shape;
    body->centroid = polygon_centroid(shape);
    assert(mass >= 0);
    body->mass = mass;
    body->color = color;
    body->velocity = VEC_ZERO;
    body->force = VEC_ZERO;
    body->impulse = VEC_ZERO;
    body->info = info;
    body->info_freer = info_freer;
    body->removed = false;
    body->angle = 0.0;
    body->collided_with = false;
    return body;
}

/**
Calls body_init_with_info, passing info fields as NULL.
 */
body_t *body_init(list_t *shape, double mass, rgb_color_t color) {
  return body_init_with_info(shape, mass, color, NULL, NULL);
}


/**
Releases the memory allocated for a body.
  */
void body_free(body_t *body) {
  list_free(body->points);
  if(body->info_freer != NULL){ // IS THIS THE RIGHT THING TO DO
    body->info_freer(body->info);
  }
  free(body);
}

/**
Gets the current shape of a body. Returns a newly allocated vector list,
which must be list_free()d.
 */
list_t *body_get_shape(body_t *body) {
  list_t *points = body->points;
  size_t n = list_size(points);
  list_t *returnList = list_init(n, (free_func_t) vec_free);
  for(size_t i = 0; i < n; i++){
    vector_t *newVec = malloc(sizeof(vector_t));
    assert(newVec!= NULL);
    vector_t *oldVec = list_get(points, i);
    newVec->x = oldVec->x;
    newVec->y = oldVec->y;
    list_add(returnList, newVec);
  }
  return returnList;
}

/**
Gets the current center of mass of a body.
*/
vector_t body_get_centroid(body_t *body) {
  return body->centroid;
}

/**
Gets the current velocity of a body.
*/
vector_t body_get_velocity(body_t *body) {
  return body->velocity;
}

/**
Gets the mass of a body.
*/
double body_get_mass(body_t *body) {
  return body->mass;
}

/**
Gets the display color of a body.
*/
rgb_color_t body_get_color(body_t *body) {
  return body->color;
}

/**
 * Gets the information associated with a body.
 */
void *body_get_info(body_t *body){
  return body->info;
}

/**
Translates a body to a new position.
The position is specified by the position of the body's center of mass.
*/
void body_set_centroid(body_t *body, vector_t x) {
  vector_t translationVector = vec_negate(body_get_centroid(body));
  polygon_translate(body->points, translationVector);
  body->centroid.x = x.x;
  body->centroid.y = x.y;
  polygon_translate(body->points, x);
}

/**
Sets a body's point list.
*/
void body_set_points(body_t *body, list_t *list) {
  body->points = list;
  body->centroid = polygon_centroid(list);
}

/**
Changes a body's velocity (the time-derivative of its position).
*/
void body_set_velocity(body_t *body, vector_t v) {
  body->velocity = v;
}

/**
Changes a body's orientation in the plane.
The body is rotated about its center of mass.
Note that the angle is *absolute*, not relative to the current orientation.
*/
void body_set_rotation(body_t *body, double new_angle) {
  double angle_to_rotate = new_angle - body->angle;
  polygon_rotate(body->points, angle_to_rotate, body->centroid);
  body->angle = new_angle;
}

void body_add_force(body_t *body, vector_t force) {
  body->force = vec_add(body->force, force);
}


void body_add_impulse(body_t *body, vector_t impulse) {
  body->impulse = vec_add(body->impulse, impulse);
}

/**
Moves a body at its current velocity over a given time interval.
*/
void body_tick(body_t *body, double dt) {
  if(!body_is_removed(body)){
    vector_t oldVelocity = body_get_velocity(body);
    double impulse_x = body->impulse.x + (dt * (body->force.x));
    double impulse_y = body->impulse.y + (dt * (body->force.y));
    vector_t newVelocity = {oldVelocity.x + ((1.0/body->mass) * impulse_x),
      oldVelocity.y + ((1.0/body->mass) * impulse_y)};
    body_set_velocity(body, newVelocity);
    double translateX = (body->velocity.x + oldVelocity.x) * dt * ACC_MULT;
    double translateY = (body->velocity.y + oldVelocity.y) * dt * ACC_MULT;
    vector_t translationVector = {translateX, translateY};
    body_set_centroid(body, vec_add(body_get_centroid(body), translationVector));
    body->force = VEC_ZERO;
    body->impulse = VEC_ZERO;
  }
}

/**
 * Marks a body for removal--future calls to body_is_removed() will return true.
 * Does not free the body.
 * If the body is already marked for removal, does nothing.
 */
void body_remove(body_t *body) {
  body->removed = true;
}

/**
 * Returns whether a body has been marked for removal.
 * This function returns false until body_remove() is called on the body,
 * and returns true afterwards.
 */
bool body_is_removed(body_t *body){
  return body->removed;
}

void body_collided(body_t *body, bool is_collided) {
  body->collided_with = is_collided;
}

bool body_is_collided(body_t *body) {
  return body->collided_with;
}
