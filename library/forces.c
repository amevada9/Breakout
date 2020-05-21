#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "forces.h"
#include "vector.h"
#include "body.h"
#include "scene.h"
#include "collision.h"

const double ELASTICITY_TERM = 1.0;
const double FALSE_CONSTANT = -1.0;

typedef struct auxillary {
  double constant; // Constnat needed
  list_t *bodies; // Bodies needed for force
  collision_handler_t collision;
  void *aux;
} auxillary_t;


auxillary_t *aux_init(double constant) {
  auxillary_t *new_aux = malloc(sizeof(auxillary_t));
  new_aux->bodies = list_init(2, (free_func_t) body_free);
  new_aux->constant = constant;
  new_aux->collision = NULL;
  return new_aux;
}

void aux_set_collision(auxillary_t *aux, collision_handler_t collide) {
  aux->collision = collide;
}

collision_handler_t aux_get_collision(auxillary_t *aux){
  return aux->collision;
}

void aux_set_aux(auxillary_t *old_aux, void *new_aux) {
  old_aux->aux = new_aux;
}

void *aux_get_aux(auxillary_t *new_aux){
  return new_aux->aux;
}

void aux_add_body(auxillary_t *aux, body_t *body) {
  list_add(aux->bodies, body);
}

void aux_set_bodies(auxillary_t *aux, list_t *bodies) {
  aux->bodies = bodies;
}

body_t *aux_get_body(auxillary_t *aux, size_t index) {
  assert(index >= 0);
  assert(index < list_size(aux->bodies));
  return list_get(aux->bodies, index);
}

double aux_get_constant(auxillary_t *aux) {
  return aux->constant;
}

void aux_free(auxillary_t *aux) {
  list_free(aux->bodies);
  free(aux);
}

/**
 * Creates a gravitational force creator on a body.
 */
void create_newtonian_gravity(scene_t *scene, double G, body_t *body1,
  body_t *body2) {
  auxillary_t *aux = aux_init(G);
  aux_add_body(aux, body1);
  aux_add_body(aux, body2);
  list_t *bodies = list_init(2, (free_func_t) body_free);
  list_add(bodies, body1);
  list_add(bodies, body2);

  scene_add_bodies_force_creator(scene, (force_creator_t) calc_gravity_force,
  aux, bodies,(free_func_t) aux_free);
}

/**
 * Calculates the gravitational force. Helper funciton to
 * create_newtonian_gravity.
 */
void calc_gravity_force(void *aux) {
  double G = aux_get_constant(aux);
  body_t *body1 = aux_get_body(aux, 0);
  body_t *body2 = aux_get_body(aux, 1);

  vector_t centroid_body_1 = body_get_centroid(body1);
  vector_t centroid_body_2 = body_get_centroid(body2);

  vector_t difference_centroids = vec_subtract(centroid_body_2, centroid_body_1);
  double distance_bodies = vec_magnitude(difference_centroids);
  // Calculating values needed to find mutual grav force.
  double mass_body_1 = body_get_mass(body1);
  double mass_body_2 = body_get_mass(body2);
  double magnitude_of_force = (G * ((double) mass_body_1 *
    (double) mass_body_2)) /
    (distance_bodies * distance_bodies * distance_bodies);

  vector_t force = vec_multiply(magnitude_of_force, difference_centroids);

  vector_t body_tip1 = *(vector_t *) list_get(body_get_shape(body1), 0);
  vector_t body_tip2 = *(vector_t *) list_get(body_get_shape(body2), 0);

  vector_t radius1 = vec_subtract(centroid_body_1, body_tip1);
  vector_t radius2 = vec_subtract(centroid_body_2, body_tip2);

  double touching = vec_magnitude(radius1) + vec_magnitude(radius2);

  vector_t force_body_1 = vec_init(force.x, force.y);
  vector_t force_body_2 = vec_negate(force_body_1);
  // Not applying gravity if too close to avoid very large numbers.
  if(distance_bodies >= touching) {
    body_add_force(body1, force_body_1);
    body_add_force(body2, force_body_2);
  }
}


/**
 * Creates a spring force creator on a body.
 */
void create_spring(scene_t *scene, double k, body_t *body1, body_t *body2) {
  auxillary_t *aux = aux_init(k);
  aux_add_body(aux, body1);
  aux_add_body(aux, body2);
  list_t *bodies = list_init(2, (free_func_t) body_free);
  list_add(bodies, body1);
  list_add(bodies, body2);

  scene_add_bodies_force_creator(scene, (force_creator_t) calc_spring_force,
  aux, bodies, (free_func_t) aux_free);

}

/**
 * Calculates the spring force. Helper funciton to create_spring.
 */
void calc_spring_force(void *aux) {
  double k = aux_get_constant(aux);
  body_t *body1 = aux_get_body(aux, 0);
  body_t *body2 = aux_get_body(aux, 1);

  vector_t centroid1 = body_get_centroid(body1);
  vector_t centroid2 = body_get_centroid(body2);
  vector_t difference = vec_subtract(centroid1, centroid2);
  vector_t force_on_2 = vec_multiply(k, difference);
  vector_t force_on_1 = vec_negate(force_on_2);

  body_add_force(body1, force_on_1);
  body_add_force(body2, force_on_2);

}

/**
 * Creates a drag force creator on a body.
 */
void create_drag(scene_t *scene, double gamma, body_t *body) {
  auxillary_t *aux = aux_init(gamma);
  aux_add_body(aux, body);
  list_t *bodies = list_init(1, (free_func_t) body_free);
  list_add(bodies, body);
  scene_add_bodies_force_creator(scene, (force_creator_t) calc_drag_force,
  aux, bodies, (free_func_t) aux_free);

}
/**
 * Calculates the drag force. Helper funciton to create_drag.
 */
void calc_drag_force(void *aux) {
  double gamma = aux_get_constant(aux);
  body_t *body = aux_get_body(aux, 0);
  vector_t velo = body_get_velocity(body);
  vector_t drag = vec_multiply(-1.0 * gamma, velo);
  body_add_force(body, drag);
}

void create_collision(
    scene_t *scene,
    body_t *body1,
    body_t *body2,
    collision_handler_t handler,
    void *aux,
    free_func_t freer
)
  {
    auxillary_t *new_aux = aux_init(FALSE_CONSTANT);
    aux_add_body(new_aux, body1);
    aux_add_body(new_aux, body2);
    list_t *bodies = list_init(1, (free_func_t) body_free);
    list_add(bodies, body1);
    list_add(bodies, body2);
    aux_set_collision(new_aux, handler);
    aux_set_aux(new_aux, aux);

    scene_add_bodies_force_creator(scene, (force_creator_t) calc_collison,
    new_aux, bodies, freer);
  }

void calc_collison(void *aux) {
  body_t *body1 = aux_get_body(aux, 0);
  body_t *body2 = aux_get_body(aux, 1);
  collision_info_t info = find_collision(body1, body2);
  void *other_aux = aux_get_aux(aux);
  if (get_if_collided(info)) {
    vector_t axis = get_collision_axis(info);
    (aux_get_collision(aux))(body1, body2, axis, other_aux);
  }
}

/**
 * Adds a force creator to a scene that destroys two bodies when they collide.
 * The bodies should be destroyed by calling body_remove().
 */
void create_destructive_collision(scene_t *scene, body_t *body1, body_t *body2){
    create_collision(scene, body1, body1, calc_destructive_force, NULL, NULL);
}

/**
 * Marks the bodies as removed as passed in from create_destructive_collision.
 */
void calc_destructive_force(body_t *body1, body_t *body2, vector_t axis, void *aux) {
    body_remove(body1);
    body_remove(body2);
}

void create_physics_collision(
    scene_t *scene,
    double elasticity,
    body_t *body1,
    body_t *body2
){
    auxillary_t *new_aux = aux_init(elasticity);
    aux_add_body(new_aux, body1);
    aux_add_body(new_aux, body2);
    list_t *bodies = list_init(1, (free_func_t) body_free);
    list_add(bodies, body1);
    list_add(bodies, body2);

    create_collision(scene, body1, body2, calc_physics_collision, new_aux, NULL);
  }

void calc_physics_collision(body_t *body1, body_t *body2, vector_t axis, void *aux) {
  if (!body_is_collided(body1) && !body_is_collided(body2)) {
    double elasticity = aux_get_constant(aux);
    double body1_mass = body_get_mass(body1);
    double body2_mass = body_get_mass(body2);
    //reduced mass is calculated, taking into account the cases where one body
    //has mass infinity
    double reduced_mass = 0.0;
    if (body1_mass == INFINITY){
      reduced_mass = body2_mass;
    }
    else if(body2_mass == INFINITY){
      reduced_mass = body1_mass;
    }
    else{
        reduced_mass = (body1_mass * body2_mass) / (body1_mass + body2_mass);
    }

    double elasticity_term = ELASTICITY_TERM + elasticity;
    double body1_velocity = vec_dot(body_get_velocity(body1), axis);
    double body2_velocity = vec_dot(body_get_velocity(body2), axis);
    double velocity_term = body2_velocity - body1_velocity;

    double impulse_magnitude = reduced_mass * elasticity_term * velocity_term;
    vector_t impulse = vec_multiply(impulse_magnitude, axis);
    body_add_impulse(body1, impulse);
    body_add_impulse(body2, vec_multiply(-1, impulse));
    body_collided(body1, true);
    body_collided(body2, true);
  }
}
