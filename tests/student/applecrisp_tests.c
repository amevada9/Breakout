#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "forces.h"
#include "test_util.h"
#include "body.h"

// implemented force creators: gravity between two bodies, a spring force
//between two bodies, and drag on a single body.

// examples:
// Tests that a conservative force (gravity) conserves K + U
// Tests that a mass on a spring oscillates like A cos(sqrt(K / M) * t)

// three ideas:
// conservation of energy for mass on a spring
// conservation of energy for drag
// balance spring forces: test if 2 masses of the same size and spring
//	 constant will cancell each other out when connected to the same mass

list_t *make_shape() {
    list_t *shape = list_init(4, free);
    vector_t *v = malloc(sizeof(*v));
    *v = (vector_t){-1, -1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){+1, -1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){+1, +1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){-1, +1};
    list_add(shape, v);
    return shape;
}

// helper method
double gravity_potential(double G, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2),
                    body_get_centroid(body1));
    return -G * body_get_mass(body1) * body_get_mass(body2) /
                  sqrt(vec_dot(r, r));
}

// helper method
double kinetic_energy(body_t *body) {
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

// helper method
double vec_mag(vector_t vec){
    double mag = sqrt(vec_dot(vec, vec));
    return mag;
}

// helper method
double drag_potential(body_t *body,
                        double GAMMA, double DT) {
    // delta_PE = work done by force = -1 * GAMMA * pow(velocity, 2)
    // where x is the distance stretched
    double drag_mag = pow(vec_mag(body_get_velocity(body)), 2) *
                                          GAMMA * -1;
    double work =  DT * drag_mag;
    return work;
}

double spring_potential(body_t *body, double K, double DT) {
    // PE = (1/2) * K * pow(x, 2)
    // where x is the distance stretched
    double mag_vel = vec_mag(body_get_velocity(body));
    double distance =  DT * mag_vel;
    return K * pow(distance, 2) / 2;
}

// Tests two opposing spring forces will cancel each other out
void test_spring_cancellation_force() {
    const double M1 = 2.0;
    const double DT = 1e-6;
    const int STEPS = 1000000;

    scene_t *scene = scene_init();
    body_t *left = body_init(make_shape(), M1, (rgb_color_t)
                              {(float) 0, (float) 0,  (float) 0});
    body_set_centroid(left, (vector_t) {-10, -10});
    body_t *center = body_init(make_shape(), M1, (rgb_color_t)
                              {(float) 0, (float) 0,  (float) 0});
    body_t *right = body_init(make_shape(), M1, (rgb_color_t)
                              {(float) 0, (float) 0,  (float) 0});
    body_set_centroid(right, (vector_t) {10, 10});
    scene_add_body(scene, right);
    scene_add_body(scene, center);
    scene_add_body(scene, left);

    create_spring(scene, 4, scene_get_body(scene, 0),
                      scene_get_body(scene, 1));

    create_spring(scene, 4, scene_get_body(scene, 1),
                      scene_get_body(scene, 2));

    for (int i = 0; i < STEPS; i++) {
        scene_tick(scene, DT);
      	vector_t vec = vec_subtract(body_get_centroid(scene_get_body(scene, 0)),
                        vec_negate(body_get_centroid(scene_get_body(scene, 2))));
      	assert(vec_equal(body_get_centroid(scene_get_body(scene, 1)),
                  (vector_t) {0,0}));
      	assert(within(1e-4, vec.x, 0) && within(1e-4, vec.y, 0));
    }
    scene_free(scene);
}

// Tests that spring force oscillations conserve K + U
void test_spring_cons_energy(){

    const double M1 = 2.0;
    const double K = 2;
    const double A = 3;
    const double DT = 1e-6;
    const int STEPS = 1000000;

    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M1, (rgb_color_t)
                          {(float) 0, (float) 0,  (float) 0});
    body_set_centroid(mass, (vector_t){A, 0});
    scene_add_body(scene, mass);
    body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t)
                          {(float) 0, (float) 0,  (float) 0});
    scene_add_body(scene, anchor);
    create_spring(scene, K, mass, anchor);
    for (int i = 0; i < STEPS; i++){
        double energy_pre =  kinetic_energy(mass);
        scene_tick(scene, DT);
        double energy_post = spring_potential(mass, K, DT)
                                                + kinetic_energy(mass);
	assert(within(1e-4, energy_pre - energy_post, 0));
    }
    scene_free(scene);
}

// Tests that drag force correctly calculates work
void test_drag_work(){
    const double M1 = 2.0;
    const double A = 3;
    const double GAMMA = 0.5;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    const double drag = 0.1;

    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M1, (rgb_color_t)
                      {(float) 0, (float) 0,  (float) 0});
    body_set_centroid(mass, (vector_t){A, 0});
    create_drag(scene, drag, mass);
    scene_add_body(scene, mass);
    for (int i = 0; i < STEPS; i++) {
        double energy_pre = kinetic_energy(mass);
        scene_tick(scene, DT);
        double energy_post = drag_potential(mass, GAMMA,
                                              DT) + kinetic_energy(mass);
        assert(within(1e-4, energy_pre - energy_post, 0));
    }
    scene_free(scene);
}

int main(int argc, char *argv[]) {
  // Run all tests if there are no command-line arguments
  bool all_tests = argc == 1;
  // Read test name from file
  char testname[100];
  if (!all_tests) {
      read_testname(argv[1], testname, sizeof(testname));
  }

  DO_TEST(test_drag_work);
  DO_TEST(test_spring_cons_energy);
  DO_TEST(test_spring_cancellation_force);

  puts("forces_test PASS");
  }
