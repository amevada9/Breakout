// TODO: IMPLEMENT YOUR TESTS IN THIS FILE
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "forces.h"
#include "test_util.h"

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

double gravity_potential(double G, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return -G * body_get_mass(body1) * body_get_mass(body2) / sqrt(vec_dot(r, r));
}
double kinetic_energy(body_t *body) {
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

// Returns the center of mass relative to the center
vector_t center_of_mass(body_t *b1, body_t *b2){
  vector_t coords1 = body_get_centroid(b1);
  vector_t coords2 = body_get_centroid(b2);
  double mass1 = body_get_mass(b1);
  double mass2 = body_get_mass(b2);
  double com_x = (mass1*coords1.x + mass2*coords2.x)/(mass1 + mass2);
  double com_y = (mass1*coords1.y + mass2*coords2.y)/(mass1 + mass2);
  return (vector_t){com_x, com_y};
}

void test_newtonian_gravity(){
  puts("------------------------Gravity tests--------------------------------");
  // check that doesn't move within buffer
  const double M1 = 5, M2 = 5;
  const double G = 1e3;
  const double DT = 1e-6;
  int STEPS = 100;
  scene_t *scene = scene_init();
  body_t *body1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
  body_set_centroid(body1, (vector_t){0.0005, 0.0005});
  scene_add_body(scene, body1);
  body_t *body2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
  body_set_centroid(body2, (vector_t){0, 0});
  scene_add_body(scene, body2);
  create_newtonian_gravity(scene, G, body1, body2);
  for (int i = 0; i < STEPS; i++) {
     assert(isclose(body_get_centroid(body1).x, .0005));
     assert(isclose(body_get_centroid(body1).y, .0005));
     assert(isclose(body_get_centroid(body2).x, 0));
     assert(isclose(body_get_centroid(body2).y, 0));
     scene_tick(scene, DT);
  }
  puts("Does not move when objects very close!");
  // make sure center of mass stays in same spot no matter what
  body_set_centroid(body1, (vector_t){10, 10});
  body_set_centroid(body2, (vector_t){-10, -10});
  STEPS = 1000000;
  for (int i = 0; i < STEPS; i++) {
       vector_t cmass = center_of_mass(body1, body2);
       assert(isclose(cmass.x, 0.0));
       assert(isclose(cmass.y, 0.0));
       scene_tick(scene, DT);
  }
  puts("Center of mass is stationary!");
  scene_free(scene);
}

double elastic_potential(double K, body_t *body1, body_t *body2){
  double ans = 0.0;
  vector_t coords1 = body_get_centroid(body1);
  vector_t coords2 = body_get_centroid(body2);
  vector_t distance = vec_subtract(coords1, coords2);
  double dist = pow(pow(distance.x, 2) + pow(distance.y, 2), 0.5);
  ans += K * pow(dist,2)/2;
  return ans;
}

void test_spring(){
  puts("-------------------------Spring tests--------------------------------");
  // test center of mass stationary and energy is conserved
  const double M1 = 5, M2 = 5;
  const double K = 1e3;
  const double DT = 1e-6;
  int STEPS = 100000;
  scene_t *scene = scene_init();
  body_t *body1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
  body_set_centroid(body1, (vector_t){10, 10});
  scene_add_body(scene, body1);
  body_t *body2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
  body_set_centroid(body2, (vector_t){0, 0});
  scene_add_body(scene, body2);
  create_spring(scene, K, body1, body2);
  double initial_energy = elastic_potential(K, body1, body2);
  for (int i = 0; i < STEPS; i++){
    // testthat the center of mass is stationary
    vector_t cmass = center_of_mass(body1, body2);
    assert(isclose(cmass.x, 5.0));
    assert(isclose(cmass.y, 5.0));
    // test that energy is conserved
    double energy = kinetic_energy(body1) + kinetic_energy(body2) + elastic_potential(K, body1, body2);
    assert(within(1e-4, energy / initial_energy, 1));
    scene_tick(scene, DT);
  }
  puts("Center of mass is stationary!");
  puts("Energy is conserved!");
  scene_free(scene);
}



void test_drag(){
  puts("---------------------------Drag tests--------------------------------");
  // 0% drag
  const double M = 5;
  const double GAMMA1 = 0;
  const double DT = 1e-6;
  const double VEL = 10;
  const double DRAG_NUM = 10;
  int STEPS = 1000000;
  scene_t *scene = scene_init();
  body_t *body = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
  body_set_centroid(body, (vector_t){0, 0});
  body_set_velocity(body, (vector_t){VEL, VEL});
  scene_add_body(scene, body);
  create_drag(scene, GAMMA1, body);
  for (int i = 0; i < STEPS; i++){
    assert(isclose(body_get_velocity(body).x, VEL));
    assert(isclose(body_get_velocity(body).y, VEL));
    scene_tick(scene, DT);
  }
  puts("0% drag does not alter velocity of body!");
  scene_free(scene);
  // 100 % drag
  scene = scene_init();
  body = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
  body_set_centroid(body, (vector_t){0, 0});
  body_set_velocity(body, (vector_t){VEL, VEL});
  /*
  to stop body in one tick, set gamma such that dv = v
  dv = dt * F/m = dt * gamma * v / m
  gamma = m/dt
  */
  const double GAMMA2 = M / DT;
  scene_add_body(scene, body);
  create_drag(scene, GAMMA2, body);
  for(int i = 0; i < DRAG_NUM; ++i){
    scene_tick(scene, DT);
    assert(isclose(body_get_velocity(body).x, 0));
    assert(isclose(body_get_velocity(body).y, 0));
  }
  puts("100% drag halts body!");
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

  DO_TEST(test_newtonian_gravity)
  DO_TEST(test_spring)
  DO_TEST(test_drag)

  puts("forces_test PASS");
}
