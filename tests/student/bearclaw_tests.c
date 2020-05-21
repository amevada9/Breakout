#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "forces.h"
#include "test_util.h"
#include "scene.h"
#include "vector.h"
#include <stdio.h>

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

/* test that a body centered between two bodies of the same mass does
  not change position of the center body since the forces balance */
void test_balanced_grav() {
    const double M = 10;
    const double G = 2;
    const double A = 3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *center_mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(center_mass, (vector_t){A, 0});
    scene_add_body(scene, center_mass);
    body_t *right_mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(right_mass, (vector_t){2.0 * A, 0});
    scene_add_body(scene, right_mass);
    body_t *left_mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(left_mass, (vector_t){0, 0});
    scene_add_body(scene, left_mass);

    create_newtonian_gravity(scene, G, center_mass, right_mass);
    create_newtonian_gravity(scene, G, center_mass, left_mass);
    for (int i = 0; i < STEPS; i++){
      assert(vec_isclose(body_get_centroid(center_mass), (vector_t) {A, 0}));
      scene_tick(scene, DT);
    }
    scene_free(scene);
}

/*test that gravity and spring don't move
objects on top of each other  */
void test_close_objects() {
    const double M = 10;
    const double K = 2;
    const double G = 1e3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass2);
    create_spring(scene, K, mass1, mass2);
    create_newtonian_gravity(scene, G, mass1, mass2);
    for (int i = 0; i < STEPS; i++) {
        assert(vec_isclose(body_get_centroid(mass1),body_get_centroid(mass2)));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

/*Test that drag slows it down */
void test_drag_slowing() {
  const double M = 10;
  const double GAMMA = 4;
  const double DT = 1e-6;
  const int STEPS = 1000000;
  scene_t *scene = scene_init();
  body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
  body_set_velocity(mass, (vector_t){1.0, 1.0});
  scene_add_body(scene, mass);
  create_drag(scene, GAMMA, mass);
  vector_t vel = body_get_velocity(scene_get_body(scene, 0));
  double old_velocity = sqrt(vec_dot(vel, vel));
  double new_velocity;
  for (int i = 0; i < STEPS; i++) {
    vel = body_get_velocity(scene_get_body(scene, 0));
    new_velocity = sqrt(vec_dot(vel, vel));
    assert(new_velocity <= old_velocity);
    scene_tick(scene, DT);
    old_velocity = new_velocity;
  }
  scene_free(scene);
}

/* Orbital velocity test */

void test_orbital_velocity(){

  //const double EMASS = 5.972*pow(10, 24);
  // 5.9*10^24 kg
  const double M1 = 1300000;
  // 1.3 million E mass
  const double G = 6.6*pow(10, -20);
  const double M2 = 1;
  // 1 E mass
  const double rad = 60;
  //60 km
  double v = 3.8*pow(10, -8);
  // v = 3.8 *10^-8 km/s

  scene_t *scene = scene_init();
  body_t *mass1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
  body_t *mass2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
  scene_add_body(scene, mass1);

  body_set_centroid(mass2, (vector_t){rad, 0});
  scene_add_body(scene, mass2);

  body_set_velocity( mass2, (vector_t){0, v});
  create_newtonian_gravity(scene, G, mass1, mass2);

  const double DT = 1e-6;
  const int STEPS = 100000;

  for (int i = 0; i < STEPS; i++) {
    vector_t vel = body_get_velocity(scene_get_body(scene, 1));
    assert(isclose(v, sqrt(vec_dot(vel, vel))));
    scene_tick(scene, DT);
  }
  scene_free(scene);
//v = sqrt(G* M /r), where M = M1 + M2
}


int main(int argc, char *argv[]) {
  // Run all tests if there are no command-line arguments
  bool all_tests = argc == 1;
  // Read test name from file
  char testname[100];
  if (!all_tests) {
      read_testname(argv[1], testname, sizeof(testname));
  }

  DO_TEST(test_balanced_grav)
  DO_TEST(test_close_objects)
  DO_TEST(test_drag_slowing)
  DO_TEST(test_orbital_velocity)

  puts("student_forces_test PASS");
}
