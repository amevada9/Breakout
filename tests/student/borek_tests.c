// TODO: IMPLEMENT YOUR TESTS IN THIS FILE
#include "body.h"
#include "test_util.h"
#include "forces.h"
#include "color.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

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

double spring_potential(double k, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return 0.5 * k * vec_dot(r, r);
}

double gravity_potential(double G, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return -G * body_get_mass(body1) * body_get_mass(body2) / sqrt(vec_dot(r, r));
}

double kinetic_energy(body_t *body) {
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

void test_gravity(){
  double M1 = 4.5, M2 = 7.3;
  const double G = 1e3;
  const double DT = 1e-6;
  const int STEPS = 1000000;
  scene_t *scene = scene_init();
  body_t *mass1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
  body_set_centroid(mass1, (vector_t){0, 0});
  scene_add_body(scene, mass1);
  body_t *mass2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
  body_set_centroid(mass2, (vector_t){20, 10});
  scene_add_body(scene,mass2);
  create_newtonian_gravity(scene, G, mass1, mass2);
  double initial_energy = gravity_potential(G, mass1, mass2);
  for (int i = 0; i < STEPS; i++) {
    double energy = gravity_potential(G, mass1, mass2) + kinetic_energy(mass1) +
                    kinetic_energy(mass2);
    assert(within(1e-4, energy / initial_energy, 1));
    scene_tick(scene, DT);
  }
  scene_free(scene);
}

void test_spring(){
  const int STEPS = 1000000;
  const double DT = 1e-6;
  const double K = 2;
  const double M1 = 4.5;
  scene_t *scene = scene_init();
  body_t *mass1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
  body_set_centroid(mass1, (vector_t){50, 50});
  scene_add_body(scene, mass1);
  body_t *mass2 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
  body_set_centroid(mass2, (vector_t){10, 20});
  scene_add_body(scene, mass2);
  create_spring(scene, K, mass1, mass2);
  double initial_energy = spring_potential(K, mass1, mass2);
  for (int i = 0; i < STEPS; i++) {
      double energy = spring_potential(K, mass1, mass2) + kinetic_energy(mass1) +
                      kinetic_energy(mass2);
      assert(within(1e-4, energy / initial_energy, 1));
      scene_tick(scene, DT);
  }
  scene_free(scene);
}

void test_drag(){
  const int STEPS = 1000000;
  const double DT = 1e-6;
  const double M1 = 4.5;
  const double GAMMA = .2;
  scene_t *scene = scene_init();
  body_t *mass = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
  vector_t initial_velocity = {10, 0};
  body_set_velocity(mass, initial_velocity);
  scene_add_body(scene, mass);
  create_drag(scene, GAMMA, mass);
  for (int i = 0; i < STEPS; i++){
    assert(within(1e-4, body_get_velocity(mass).x, 10 * exp(-GAMMA * i * DT / M1)));
    scene_tick(scene, DT);
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

 DO_TEST(test_spring)
 DO_TEST(test_gravity)
 DO_TEST(test_drag)

 puts("student_test PASS");
}
